# Load the shared ideal implementation into the current parent module once.
if !isdefined(@__MODULE__, :IdealEntanglementSwitch)
    include(joinpath(@__DIR__, "IdealEntanglementSwitch.jl"))
end

module NonIdealEntanglementSwitch

using QuantumSavory
using QuantumSavory.ProtocolZoo: EntanglerProt, EntanglementTracker
using QuantumSavory.StatesZoo: DepolarizedBellPair

using Graphs
using Random
using Statistics

using ..IdealEntanglementSwitch:
    SwitchStatistics,
    OldestFirstSwitchProt,
    consume_entangled_pairs!,
    steady_state_throughput

const ConcurrentSim = QuantumSavory.ConcurrentSim

import .ConcurrentSim: Process, run
import QuantumSavory: @process

export retention_to_lifetime,
    prepare_nonideal_switch,
    run_nonideal_switch

"""
    retention_to_lifetime(p_w)

Convert the per-time-unit memory retention factor `p_w` into the characteristic
time `τ` expected by `QuantumSavory.Depolarization`.

# Arguments
- `p_w::Float64`: probability that the stored state is preserved during one
  time unit.

# Returns
- `Float64`: the lifetime `τ = -1 / log(p_w)`. For `p_w == 1`, returns `Inf`.
"""
function retention_to_lifetime(p_w::Float64)
    isone(p_w) && return Inf
    return -1 / log(p_w)
end

"""
    memory_register(number_of_slots, p_w)

Create a qubit register whose slots all have the same memory depolarization.

# Arguments
- `number_of_slots::Int`: number of qubit memories in the register.
- `p_w::Float64`: state-retention factor per simulated time unit.

# Returns
- `Register`: a register with one identical background process per slot. When
  `p_w == 1`, the slots have no background noise.
"""
function memory_register(number_of_slots::Int, p_w::Float64)
    background = isone(p_w) ? nothing : Depolarization(retention_to_lifetime(p_w))
    return Register(
        fill(Qubit(), number_of_slots),
        fill(background, number_of_slots),
    )
end

"""
    prepare_nonideal_switch(; number_of_users=4, w=0.9,
                              p_success=0.5, p_w=0.99, seed=1234)

Build the non-ideal star network and start all discrete-event protocols. The
oldest-first controller, entanglement consumer, ideal BSM, zero-delay classical
communication, and statistics are reused from the ideal implementation.

# Keyword arguments
- `number_of_users::Int`: number `k` of client nodes and physical links.
- `w::Float64`: Werner parameter of every newly generated link-level pair,
  `ρ = w|Φ⁺><Φ⁺| + (1-w)I/4`.
- `p_success::Float64`: success probability of each generation attempt.
- `p_w::Float64`: memory-retention factor per time unit, shared by every client
  and switch memory.
- `seed::Int`: random seed used for generation attempts and measurements.

# Returns
A named tuple containing the simulation, network, statistics, switch index,
client indices, and oldest-first controller. The simulation processes have been
started but simulated time has not yet advanced.
"""
function prepare_nonideal_switch(;
        number_of_users::Int=4,
        w::Float64=0.9,
        p_success::Float64=0.5,
        p_w::Float64=0.99,
        seed::Int=1234,
    )
    Random.seed!(seed)

    switch_node = 1
    client_nodes = collect(2:(number_of_users + 1))
    graph = star_graph(number_of_users + 1)

    # Apply the same memory background to every physical memory.
    switch_register = memory_register(number_of_users, p_w)
    client_registers = [memory_register(1, p_w) for _ in client_nodes]
    net = RegisterNet(
        graph,
        [switch_register, client_registers...];
        classical_delay=0.0,
    )
    sim = get_time_tracker(net)
    statistics = SwitchStatistics()

    # Every successful link attempt produces the requested Werner state.
    link_state = DepolarizedBellPair(w)
    for (switch_slot, client_node) in enumerate(client_nodes)
        entangler = EntanglerProt(
            sim,
            net,
            switch_node,
            client_node;
            pairstate=link_state,
            success_prob=p_success,
            attempt_time=1.0,
            local_busy_time_pre=0.0,
            local_busy_time_post=0.0,
            retry_lock_time=nothing,
            chooseslotA=switch_slot,
            chooseslotB=1,
        )
        @process entangler()
    end

    # Track BSM outcomes and apply the corresponding Pauli corrections.
    for client_node in client_nodes
        @process EntanglementTracker(sim, net, client_node)()
    end

    # Consume any client pair immediately after measuring its fidelity.
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

    # Reuse the same oldest-first swapping policy as in the ideal model.
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
    run_nonideal_switch(; number_of_users=4, simulation_time=10_001.0,
                         warmup_time=1_000.0, w=0.9, p_success=0.5,
                         p_w=0.99, seed=1234)

Prepare and run one non-ideal entanglement-switch simulation, then compute the
three aggregate performance metrics required by the assignment.

# Keyword arguments
- `number_of_users::Int`: number `k` of clients.
- `simulation_time::Float64`: final simulated time.
- `warmup_time::Float64`: initial interval excluded from throughput.
- `w::Float64`: Werner parameter of newly generated link pairs.
- `p_success::Float64`: success probability per one-unit attempt.
- `p_w::Float64`: memory-retention factor per simulated time unit.
- `seed::Int`: random seed for this run.

# Returns
A named tuple containing the prepared network and processes, all input
parameters, and `throughput`, `mean_fidelity`, and `mean_waiting_time`.
"""
function run_nonideal_switch(;
        number_of_users::Int=4,
        simulation_time::Float64=10_001.0,
        warmup_time::Float64=1_000.0,
        w::Float64=0.9,
        p_success::Float64=0.5,
        p_w::Float64=0.99,
        seed::Int=1234,
    )
    setup = prepare_nonideal_switch(;
        number_of_users,
        w,
        p_success,
        p_w,
        seed,
    )
    run(setup.sim, simulation_time)

    statistics = setup.statistics
    throughput = steady_state_throughput(
        statistics,
        simulation_time;
        warmup_time,
    )
    mean_fidelity = isempty(statistics.fidelities) ? NaN : mean(statistics.fidelities)
    mean_waiting_time =
        isempty(statistics.waiting_times) ? NaN : mean(statistics.waiting_times)

    return (;
        setup...,
        number_of_users,
        simulation_time,
        warmup_time,
        w,
        p_success,
        p_w,
        throughput,
        mean_fidelity,
        mean_waiting_time,
    )
end

end # module NonIdealEntanglementSwitch
