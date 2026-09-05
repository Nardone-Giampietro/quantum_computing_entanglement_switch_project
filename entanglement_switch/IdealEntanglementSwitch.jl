module IdealEntanglementSwitch

using QuantumSavory
using QuantumSavory.ProtocolZoo:
    AbstractProtocol,
    EntanglerProt,
    SwapperProt,
    EntanglementTracker,
    EntanglementCounterpart

using Graphs
using Random
using Statistics

const ConcurrentSim = QuantumSavory.ConcurrentSim
const ResumableFunctions = QuantumSavory.ResumableFunctions

import .ConcurrentSim: Process, Simulation, now, run
import .ResumableFunctions: @resumable, @yield
import QuantumSavory: @process

export PHI_PLUS,
    SwitchStatistics,
    OldestFirstSwitchProt,
    consume_entangled_pairs!,
    prepare_ideal_switch,
    steady_state_throughput,
    run_ideal_switch

# Perfect Bell state used for every link in the ideal case.
const PHI_PLUS = (Z1 ⊗ Z1 + Z2 ⊗ Z2) / sqrt(2)

"""
    SwitchStatistics()

Store the samples collected during an ideal-switch simulation.

# Fields
- `waiting_times`: elapsed time between creation of the oldest selected
  link-level pair and completion of the corresponding BSM.
- `delivery_times`: simulation times at which end-to-end pairs are consumed.
- `fidelities`: fidelities of the delivered pairs with respect to `|Phi+>`.
"""
mutable struct SwitchStatistics
    waiting_times::Vector{Float64}
    delivery_times::Vector{Float64}
    fidelities::Vector{Float64}
end

SwitchStatistics() = SwitchStatistics(
    Float64[],
    Float64[],
    Float64[],
)

"""
    OldestFirstSwitchProt(sim, net, switch_node, client_nodes, statistics)

Implement the switch controller using the oldest-first selection policy.

# Fields
- `sim`: discrete-event simulation in which the protocol runs.
- `net`: register network containing the switch and client memories.
- `switch_node`: index of the switch node in `net`.
- `client_nodes`: indices of all client nodes connected to the switch.
- `statistics`: mutable collection in which waiting-time samples are stored.
"""
struct OldestFirstSwitchProt <: AbstractProtocol
    sim::Simulation
    net::RegisterNet
    switch_node::Int
    client_nodes::Vector{Int}
    statistics::SwitchStatistics
end

"""
    (prot::OldestFirstSwitchProt)()

Continuously select the two oldest link-level Bell pairs stored at the switch
and perform an instantaneous ideal BSM on their switch-side qubits. If fewer
than two pairs are available, the process waits for a register update.

# Arguments
- `prot`: configured oldest-first switch protocol.

# Returns
A resumable scheduler process. It normally runs indefinitely and records each
waiting-time sample by mutating `prot.statistics`.
"""
@resumable function (prot::OldestFirstSwitchProt)()
    switch_register = prot.net[prot.switch_node]

    while true
        # Available pairs in the switch slots, one per connected client.
        available = queryall(
            switch_register,
            EntanglementCounterpart,
            in(prot.client_nodes),
            W,
            W;
            locked=false,
            assigned=true,
        )

        # Suspend the protocol until a second pair becomes available.
        if length(available) < 2
            @yield onchange(switch_register, Tag)
            continue
        end

        # `time` is the creation time; the slot index breaks any ties.
        sort!(available; by=result -> (result.time, result.slot.idx))
        first_pair = available[1]
        second_pair = available[2]

        first_client = first_pair.tag[2]
        second_client = second_pair.tag[2]
        oldest_creation_time = min(first_pair.time, second_pair.time)

        # The swap uses the slots selected by the oldest-first policy.
        swapper = SwapperProt(
            prot.sim,
            prot.net,
            prot.switch_node;
            chooseslots=[first_pair.slot.idx, second_pair.slot.idx],
            nodeL=first_client,
            nodeH=second_client,
            local_busy_time=0.0,
            rounds=1,
        )

        swap_process = @process swapper()
        @yield swap_process

        # Waiting time starts when the older pair was created.
        push!(
            prot.statistics.waiting_times,
            now(prot.sim) - oldest_creation_time,
        )
    end
end

"""
    consume_entangled_pairs!(sim, net, node_a, node_b, statistics)

Wait for end-to-end entanglement involving `node_a` with nother node `node_b`, measure the 
pair fidelity with respect to `|Phi+>`, record its delivery time, and immediately consume the
two qubits. Reciprocal tags are checked before collecting each sample.

# Arguments
- `sim`: discrete-event simulation in which the consumer runs.
- `net`: register network containing the client memories.
- `node_a`: client node monitored by this consumer process.
- `node_b`: client node expected at the remote end of the delivered pair.
- `statistics`: mutable collection receiving delivery-time and fidelity samples.

# Returns
A resumable process that normally runs indefinitely. Results are produced by
mutating `statistics`; no ordinary value is returned.
"""
@resumable function consume_entangled_pairs!(
        sim::Simulation,
        net::RegisterNet,
        node_a::Int,
        node_b::Int,
        statistics::SwitchStatistics,
    )
    register_a = net[node_a]
    register_b = net[node_b]

    while true
        # Wait for a pair entangling the first node with the second node.
        first = @yield query_wait(
            register_a,
            EntanglementCounterpart,
            node_b,
            W,
            W;
            locked=false,
            assigned=true,
        )

        pair_id = first.tag[4]

        # Find the reciprocal endpoint using the same pair identifier.
        second = query(
            register_b,
            EntanglementCounterpart,
            node_a,
            first.slot.idx,
            pair_id;
            locked=false,
            assigned=true,
        )

        # retry if metadata of node_b are not available 
        if isnothing(second)
            @yield onchange(register_b, Tag)
            continue
        end

        qubit_a = first.slot
        qubit_b = second.slot
        @yield lock(qubit_a) & lock(qubit_b)

        # After locking, check that both tags are still valid.
        current_a = query(
            qubit_a,
            EntanglementCounterpart,
            node_b,
            qubit_b.idx,
            pair_id;
            locked=true,
            assigned=true,
        )
        current_b = query(
            qubit_b,
            EntanglementCounterpart,
            node_a,
            qubit_a.idx,
            pair_id;
            locked=true,
            assigned=true,
        )

        if isnothing(current_a) || isnothing(current_b)
            unlock(qubit_a)
            unlock(qubit_b)
            continue
        end

        # Fidelity of the delivered pair with respect to the ideal Bell state.
        fidelity = real(
            observable(
                (qubit_a, qubit_b),
                projector(PHI_PLUS);
                time=now(sim),
            ),
        )

        # Remove the tags and clear the quantum memories for new generations.
        untag!(qubit_a, current_a.id)
        untag!(qubit_b, current_b.id)
        traceout!(qubit_a, qubit_b)

        push!(statistics.delivery_times, now(sim))
        push!(statistics.fidelities, fidelity)

        unlock(qubit_a)
        unlock(qubit_b)
    end
end

"""
    prepare_ideal_switch(; number_of_users=4, seed=1234)

Build an ideal star network and start its entanglers, pair consumers, and
oldest-first switch controller. Link generation always succeeds after one unit
of simulation time, memories are noiseless, and every generated pair is in
`|Phi+>`. The simulation clock is not advanced by this function.

# Keyword arguments
- `number_of_users`: number of client nodes connected to the central switch.
- `seed`: seed used by the simulation's random-number generator.

# Returns
A named tuple containing `sim`, `net`, `statistics`, `switch_node`,
`client_nodes`, and `controller`.
"""
function prepare_ideal_switch(; number_of_users::Int=4, seed::Int=1234)
    Random.seed!(seed)

    switch_node = 1
    client_nodes = collect(2:(number_of_users + 1))
    graph = star_graph(number_of_users + 1)

    # The switch has one slot per link; each client has a single slot.
    switch_register = Register(number_of_users)
    client_registers = [Register(1) for _ in client_nodes]
    net = RegisterNet(
        graph,
        [switch_register, client_registers...];
        classical_delay=0.0,
    )
    sim = get_time_tracker(net)
    statistics = SwitchStatistics()

    # Each link deterministically generates a perfect pair in one time unit.
    for (switch_slot, client_node) in enumerate(client_nodes)
        entangler = EntanglerProt(
            sim,
            net,
            switch_node,
            client_node;
            pairstate=PHI_PLUS,
            success_prob=1.0,
            attempt_time=1.0,
            local_busy_time_pre=0.0,
            local_busy_time_post=0.0,
            retry_lock_time=nothing,
            chooseslotA=switch_slot,
            chooseslotB=1,
        )
        @process entangler()
    end

    # Trackers propagate BSM outcomes and update the Pauli corrections.
    for client_node in client_nodes
        @process EntanglementTracker(sim, net, client_node)()
    end

    # One consumer per client pair accepts any end-to-end delivery.
    for first_index in eachindex(client_nodes)
        for second_index in (first_index + 1):length(client_nodes)
            node_a = client_nodes[first_index]
            node_b = client_nodes[second_index]
            @process consume_entangled_pairs!(
                sim,
                net,
                node_a,
                node_b,
                statistics,
            )
        end
    end

    # Central controller responsible for pair selection and swapping.
    controller = OldestFirstSwitchProt(
        sim,
        net,
        switch_node,
        client_nodes,
        statistics,
    )
    @process controller()

    return (; sim, net, statistics, switch_node, client_nodes, controller)
end

"""
    steady_state_throughput(statistics, simulation_time; warmup_time=1.0)

Compute the delivered-pair throughput after excluding an initial warm-up
interval. Deliveries in `[warmup_time, simulation_time)` are divided by the
length of that interval.

# Arguments
- `statistics`: samples collected during a switch simulation.
- `simulation_time`: final simulation time used for the measurement.

# Keyword arguments
- `warmup_time`: beginning of the interval used to compute throughput.

# Returns
The steady-state throughput in delivered pairs per unit of simulation time.
"""
function steady_state_throughput(statistics, simulation_time; warmup_time=1.0)
    simulation_time > warmup_time ||
        throw(ArgumentError("simulation_time must be greater than warmup_time"))

    delivered = count(
        time -> warmup_time <= time < simulation_time,
        statistics.delivery_times,
    )
    delivered / (simulation_time - warmup_time)
end

"""
    run_ideal_switch(; number_of_users=4, simulation_time=101.0, seed=1234)

Construct and run the ideal switch benchmark, then compute the throughput,
mean delivered-pair fidelity, and mean waiting time from the collected samples.

# Keyword arguments
- `number_of_users`: number of clients connected to the central switch.
- `simulation_time`: time at which the discrete-event simulation stops.
- `seed`: seed used by the simulation's random-number generator.

# Returns
A named tuple containing the simulation objects and statistics returned by
`prepare_ideal_switch`, together with `throughput`, `mean_fidelity`, and
`mean_waiting_time`.
"""
function run_ideal_switch(;
        number_of_users::Int=4,
        simulation_time::Float64=101.0,
        seed::Int=1234,
    )
    setup = prepare_ideal_switch(; number_of_users, seed)
    run(setup.sim, simulation_time)

    # Compute averages over all pairs delivered during the simulation.
    statistics = setup.statistics
    throughput = steady_state_throughput(statistics, simulation_time)
    mean_fidelity = isempty(statistics.fidelities) ? NaN : mean(statistics.fidelities)
    mean_waiting_time =
        isempty(statistics.waiting_times) ? NaN : mean(statistics.waiting_times)

    return (;
        setup...,
        number_of_users,
        simulation_time,
        throughput,
        mean_fidelity,
        mean_waiting_time,
    )
end

end # module IdealEntanglementSwitch
