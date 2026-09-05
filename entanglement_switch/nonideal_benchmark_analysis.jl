include(joinpath(@__DIR__, "NonIdealEntanglementSwitch.jl"))
using .NonIdealEntanglementSwitch

using CairoMakie
using Distributions: TDist, quantile
using Statistics

const ConcurrentSim = NonIdealEntanglementSwitch.ConcurrentSim

"""Derive a reproducible seed for one point of a parameter sweep."""
parameter_seed(base, row, column, repetition=0) =
    base + 10_000 * row + 100 * column + repetition

"""Return the sample mean and the half-width of its pointwise 95% t interval."""
function mean_and_ci95(values)
    mean_value = mean(values)
    half_width = quantile(TDist(length(values) - 1), 0.975) *
        std(values) / sqrt(length(values))
    return mean_value, half_width
end

"""Return true when a curve has a visually meaningful confidence interval."""
ci_is_visible(half_widths, y_span) =
    maximum(abs, half_widths) > 0.005 * y_span

"""Save a figure in `output_directory` and return its complete path."""
function save_figure(figure, output_directory, filename)
    mkpath(output_directory)
    output_path = joinpath(output_directory, filename)
    save(output_path, figure)
    return output_path
end

"""
    collect_sample_mean(sample_field; number_of_users, w, p_success, p_w,
                        target_samples, burn_in_samples, simulation_chunk, seed)

Run one simulation until enough values of `sample_field` have been collected,
discard the initial samples, and return their mean and collection metadata.
`sample_field` can be `:fidelities` or `:waiting_times`.
"""
function collect_sample_mean(
        sample_field::Symbol;
        number_of_users,
        w,
        p_success,
        p_w,
        target_samples,
        burn_in_samples,
        simulation_chunk,
        seed,
    )
    setup = prepare_nonideal_switch(;
        number_of_users,
        w,
        p_success,
        p_w,
        seed,
    )
    samples = getproperty(setup.statistics, sample_field)
    required_samples = burn_in_samples + target_samples

    while length(samples) < required_samples
        next_stop = ConcurrentSim.now(setup.sim) + simulation_chunk
        ConcurrentSim.run(setup.sim, next_stop)
    end

    retained_samples = @view samples[
        (burn_in_samples + 1):(burn_in_samples + target_samples)
    ]

    return (;
        mean_value=mean(retained_samples),
        sample_count=length(retained_samples),
        simulation_time=ConcurrentSim.now(setup.sim),
    )
end

"""
    run_fidelity_vs_w_sweep(; w_values=0.0:0.05:1.0,
                              simulation_time=21.0, repetitions=5, seed=1234)

Run a two-client switch simulation for every value of `w`. Memory noise is
disabled with `p_w=1`, while deterministic link generation (`p_success=1`)
isolates the effect of the Werner parameter on the swapping result.

# Keyword arguments
- `w_values`: Werner parameters used in the sweep.
- `simulation_time::Float64`: final time of every simulation run.
- `repetitions::Int`: number of independent runs at each point.
- `seed::Int`: base random seed; a different deterministic seed is used for
  each run.

# Returns
A named tuple containing the Werner parameters, simulated and expected mean
fidelities, 95% confidence half-widths, absolute errors, sample counts, and
fixed simulation parameters.
"""
function run_fidelity_vs_w_sweep(;
        w_values=0.0:0.05:1.0,
        simulation_time::Float64=21.0,
        repetitions::Int=5,
        seed::Int=1234,
    )
    werner_parameters = Float64.(collect(w_values))
    simulated_fidelities = Float64[]
    fidelity_ci95 = Float64[]
    expected_fidelities = Float64[]
    absolute_errors = Float64[]
    sample_counts = Int[]

    for (index, w) in enumerate(werner_parameters)
        repetition_fidelities = Float64[]
        point_sample_count = 0

        for repetition in 1:repetitions
            result = run_nonideal_switch(;
                number_of_users=2,
                simulation_time,
                warmup_time=1.0,
                w,
                p_success=1.0,
                p_w=1.0,
                seed=parameter_seed(seed, 1, index, repetition),
            )
            push!(repetition_fidelities, result.mean_fidelity)
            point_sample_count += length(result.statistics.fidelities)
        end

        simulated_fidelity, fidelity_half_width =
            mean_and_ci95(repetition_fidelities)

        expected_fidelity = (1 + 3w^2) / 4
        push!(simulated_fidelities, simulated_fidelity)
        push!(fidelity_ci95, fidelity_half_width)
        push!(expected_fidelities, expected_fidelity)
        push!(absolute_errors, abs(simulated_fidelity - expected_fidelity))
        push!(sample_counts, point_sample_count)

        println(
            "w=$w: simulated=$simulated_fidelity ± $fidelity_half_width, " *
            "expected=$expected_fidelity (95% CI)",
        )
    end

    return (;
        werner_parameters,
        simulated_fidelities,
        fidelity_ci95,
        expected_fidelities,
        absolute_errors,
        sample_counts,
        simulation_time,
        repetitions,
        number_of_users=2,
        p_success=1.0,
        p_w=1.0,
    )
end

"""
    test_fidelity_vs_w_sweep(results; atol=1.0e-12)

Check that every simulated mean fidelity agrees with the analytical swapping
formula within the requested absolute tolerance.

# Arguments
- `results`: output returned by `run_fidelity_vs_w_sweep`.

# Keyword arguments
- `atol`: maximum accepted absolute numerical discrepancy.

# Returns
`nothing` if all points pass; otherwise an assertion error is raised.
"""
function test_fidelity_vs_w_sweep(results; atol=1.0e-12)
    @assert all(
        isapprox.(
            results.simulated_fidelities,
            results.expected_fidelities;
            atol,
            rtol=0.0,
        ),
    )
    return nothing
end

"""
    plot_fidelity_vs_w_sweep(results; output_directory)

Plot the simulated mean fidelities together with the analytical curve used by
the validation test.

# Arguments
- `results`: output returned by `run_fidelity_vs_w_sweep`.

# Keyword arguments
- `output_directory`: directory in which the PNG figure is saved.

# Returns
The path of the generated PNG file.
"""
function plot_fidelity_vs_w_sweep(results; output_directory)
    figure = Figure(size=(800, 500))
    axis = Axis(
        figure[1, 1];
        title="No memory noise: mean fidelity versus Werner parameter\n" *
            "k=$(results.number_of_users), p_success=$(results.p_success), " *
            "p_w=$(results.p_w). 95% pointwise CI: $(results.repetitions) runs",
        xlabel="Link-level Werner parameter w",
        ylabel="Mean fidelity with respect to |Phi+>",
        limits=((0.0, 1.0), (0.2, 1.02)),
    )

    lines!(
        axis,
        results.werner_parameters,
        results.expected_fidelities;
        color=:black,
        linewidth=2,
        label="Analytical: (1 + 3w²) / 4",
    )
    if ci_is_visible(results.fidelity_ci95, 1.02 - 0.2)
        errorbars!(
            axis,
            results.werner_parameters,
            results.simulated_fidelities,
            results.fidelity_ci95;
            color=(:dodgerblue, 0.65),
            whiskerwidth=8,
        )
    end
    scatter!(
        axis,
        results.werner_parameters,
        results.simulated_fidelities;
        color=:dodgerblue,
        marker=:circle,
        markersize=10,
        label="QuantumSavory simulation",
    )
    axislegend(axis; position=:lt)

    return save_figure(figure, output_directory, "mean_fidelity_vs_w.png")
end

"""
    run_throughput_vs_success_probability_sweep(;
        p_success_values=0.05:0.05:1.0,
        client_counts=(2, 4, 8),
        simulation_time=2_001.0,
        warmup_time=200.0,
        repetitions=5,
        w=0.9,
        p_w=0.99,
        seed=4321,
    )

Measure throughput as a function of link-generation success probability for
several numbers of clients. Every plotted point is the mean of independent
simulation runs, which reduces fluctuations caused by random link attempts.

# Keyword arguments
- `p_success_values`: success probabilities used on the horizontal axis.
- `client_counts`: numbers of clients for which separate curves are generated.
- `simulation_time::Float64`: final time of every simulation run.
- `warmup_time::Float64`: initial interval excluded from throughput.
- `repetitions::Int`: number of independent runs averaged at each point.
- `w::Float64`: Werner parameter of every generated link-level pair.
- `p_w::Float64`: memory-retention factor per simulated time unit.
- `seed::Int`: base seed used to derive a distinct seed for every run.

# Returns
A named tuple containing the sampled probabilities, client counts, mean
throughputs, 95% confidence half-widths, saturation values, and fixed benchmark
parameters.
"""
function run_throughput_vs_success_probability_sweep(;
        p_success_values=0.05:0.05:1.0,
        client_counts=(2, 4, 8),
        simulation_time::Float64=2_001.0,
        warmup_time::Float64=200.0,
        repetitions::Int=5,
        w::Float64=0.9,
        p_w::Float64=0.99,
        seed::Int=4321,
    )
    success_probabilities = Float64.(collect(p_success_values))
    numbers_of_clients = Int.(collect(client_counts))
    mean_throughputs = Matrix{Float64}(
        undef,
        length(numbers_of_clients),
        length(success_probabilities),
    )
    throughput_ci95 = similar(mean_throughputs)

    for (client_index, number_of_users) in enumerate(numbers_of_clients)
        for (probability_index, p_success) in enumerate(success_probabilities)
            repetition_throughputs = Float64[]

            for repetition in 1:repetitions
                run_seed = parameter_seed(
                    seed,
                    client_index,
                    probability_index,
                    repetition,
                )
                result = run_nonideal_switch(;
                    number_of_users,
                    simulation_time,
                    warmup_time,
                    w,
                    p_success,
                    p_w,
                    seed=run_seed,
                )
                push!(repetition_throughputs, result.throughput)
            end

            mean_throughput, throughput_half_width =
                mean_and_ci95(repetition_throughputs)
            mean_throughputs[client_index, probability_index] = mean_throughput
            throughput_ci95[client_index, probability_index] =
                throughput_half_width

            println(
                "k=$number_of_users, p_success=$p_success: " *
                "mean throughput=$mean_throughput ± $throughput_half_width " *
                "(95% CI)",
            )
        end
    end

    saturation_throughputs = numbers_of_clients ./ 2

    return (;
        success_probabilities,
        client_counts=numbers_of_clients,
        mean_throughputs,
        throughput_ci95,
        saturation_throughputs,
        simulation_time,
        warmup_time,
        repetitions,
        w,
        p_w,
    )
end

"""
    plot_throughput_vs_success_probability(results; output_directory)

Plot one throughput curve for each number of clients, together with the
corresponding horizontal dashed saturation line.

# Arguments
- `results`: output of `run_throughput_vs_success_probability_sweep`.

# Keyword arguments
- `output_directory`: directory in which the PNG figure is saved.

# Returns
The path of the generated PNG file.
"""
function plot_throughput_vs_success_probability(results; output_directory)
    figure = Figure(size=(900, 560))
    axis = Axis(
        figure[1, 1];
        title="Throughput versus link success probability\n" *
            "w=$(results.w), p_w=$(results.p_w), " *
            "T=$(results.simulation_time), warm-up=$(results.warmup_time). " *
            "95% pointwise CI: $(results.repetitions) runs",
        xlabel="Success probability per attempt",
        ylabel="Throughput [end-to-end pairs / time unit]",
        limits=((0.05, 1.0), (0.0, 4.2)),
    )

    colors = (:dodgerblue, :darkorange, :seagreen)
    markers = (:circle, :rect, :utriangle)
    throughput_span = 4.2

    for (index, number_of_users) in enumerate(results.client_counts)
        color = colors[index]
        if ci_is_visible(results.throughput_ci95[index, :], throughput_span)
            errorbars!(
                axis,
                results.success_probabilities,
                results.mean_throughputs[index, :],
                results.throughput_ci95[index, :];
                color=(color, 0.65),
                whiskerwidth=7,
            )
        end
        scatterlines!(
            axis,
            results.success_probabilities,
            results.mean_throughputs[index, :];
            color,
            marker=markers[index],
            markersize=9,
            linewidth=2,
            label="k = $number_of_users",
        )
        hlines!(
            axis,
            [results.saturation_throughputs[index]];
            color,
            linestyle=:dash,
            linewidth=1.5,
            label="k = $number_of_users saturation",
        )
    end

    axislegend(axis; position=:lt)

    return save_figure(
        figure,
        output_directory,
        "throughput_vs_success_probability.png",
    )
end

"""
    run_fidelity_vs_success_probability_sweep(;
        p_success_values=0.05:0.05:1.0,
        p_w_values=(0.1, 0.9, 1.0),
        number_of_users=4,
        w=0.9,
        target_samples=100,
        burn_in_samples=100,
        simulation_chunk=100.0,
        repetitions=5,
        seed=6543,
    )

Measure mean end-to-end fidelity as a function of link-generation success
probability for several memory-retention factors. Every point is computed from
the same number of delivered-pair samples.

# Keyword arguments
- `p_success_values`: success probabilities shown on the horizontal axis.
- `p_w_values`: memory-retention factors represented by separate curves.
- `number_of_users::Int`: fixed number of clients connected to the switch.
- `w::Float64`: fixed Werner parameter of each generated link-level pair.
- `target_samples::Int`: number of fidelity samples retained in each run.
- `burn_in_samples::Int`: initial delivered pairs excluded from each run.
- `simulation_chunk::Float64`: amount of simulated time added while collecting
  the requested number of samples.
- `repetitions::Int`: number of independent runs at each point.
- `seed::Int`: base seed used to derive a distinct seed for every run.

# Returns
A named tuple containing the success probabilities, memory-retention factors,
simulated mean fidelities, 95% confidence half-widths, sample counts, mean final
simulation times, and fixed benchmark parameters.
"""
function run_fidelity_vs_success_probability_sweep(;
        p_success_values=0.05:0.05:1.0,
        p_w_values=(0.1, 0.9, 1.0),
        number_of_users::Int=4,
        w::Float64=0.9,
        target_samples::Int=100,
        burn_in_samples::Int=100,
        simulation_chunk::Float64=100.0,
        repetitions::Int=5,
        seed::Int=6543,
    )
    success_probabilities = Float64.(collect(p_success_values))
    memory_retention_factors = Float64.(collect(p_w_values))
    dimensions = (
        length(memory_retention_factors),
        length(success_probabilities),
    )

    mean_fidelities = Matrix{Float64}(undef, dimensions)
    fidelity_ci95 = similar(mean_fidelities)
    sample_counts = Matrix{Int}(undef, dimensions)
    simulation_times = Matrix{Float64}(undef, dimensions)

    for (memory_index, p_w) in enumerate(memory_retention_factors)
        for (probability_index, p_success) in enumerate(success_probabilities)
            repetition_means = Float64[]
            repetition_times = Float64[]
            point_sample_count = 0

            for repetition in 1:repetitions
                sample = collect_sample_mean(
                    :fidelities;
                    number_of_users,
                    w,
                    p_success,
                    p_w,
                    target_samples,
                    burn_in_samples,
                    simulation_chunk,
                    seed=parameter_seed(
                        seed,
                        memory_index,
                        probability_index,
                        repetition,
                    ),
                )
                push!(repetition_means, sample.mean_value)
                push!(repetition_times, sample.simulation_time)
                point_sample_count += sample.sample_count
            end

            mean_fidelity, fidelity_half_width =
                mean_and_ci95(repetition_means)
            mean_fidelities[memory_index, probability_index] = mean_fidelity
            fidelity_ci95[memory_index, probability_index] =
                fidelity_half_width
            sample_counts[memory_index, probability_index] = point_sample_count
            simulation_times[memory_index, probability_index] =
                mean(repetition_times)
        end

        println("Completed fidelity curve for p_w=$p_w")
    end

    return (;
        success_probabilities,
        memory_retention_factors,
        mean_fidelities,
        fidelity_ci95,
        sample_counts,
        simulation_times,
        number_of_users,
        w,
        target_samples,
        burn_in_samples,
        simulation_chunk,
        repetitions,
    )
end

"""
    plot_fidelity_vs_success_probability(results; output_directory)

Plot the simulated mean end-to-end fidelity against link-generation success
probability, using a separate curve for every memory-retention factor. Dashed
horizontal lines mark the completely depolarized lower limit and the ideal
memory upper limit for the selected Werner parameter.

# Arguments
- `results`: output returned by `run_fidelity_vs_success_probability_sweep`.

# Keyword arguments
- `output_directory`: directory in which the PNG figure is saved.

# Returns
The path of the generated PNG file.
"""
function plot_fidelity_vs_success_probability(results; output_directory)
    figure = Figure(size=(900, 560))
    axis = Axis(
        figure[1, 1];
        title="Mean fidelity versus link success probability " *
            "(pointwise 95% CI)\n" *
            "Fixed: k=$(results.number_of_users), w=$(results.w), " *
            "$(results.repetitions) runs, " *
            "$(results.target_samples) samples/run, " *
            "burn-in=$(results.burn_in_samples)",
        xlabel="Success probability per attempt",
        ylabel="Mean fidelity with respect to |Phi+>",
        limits=((0.0, 1.0), (0.0, 0.9)),
        xticks=0.0:0.2:1.0,
        yticks=0.0:0.1:0.9,
    )

    colors = (:firebrick, :dodgerblue, :black)
    markers = (:diamond, :circle, :cross)

    fidelity_min = 1 / 4
    fidelity_max = (1 + 3 * results.w^2) / 4

    hlines!(
        axis,
        [fidelity_min];
        color=:gray45,
        linestyle=:dash,
        linewidth=2,
        label="F_min = 1/4 = $(round(fidelity_min; digits=4))",
    )
    hlines!(
        axis,
        [fidelity_max];
        color=:gray45,
        linestyle=:dot,
        linewidth=2,
        label="F_max = (1+3w^2)/4 = $(round(fidelity_max; digits=4))",
    )

    for (index, p_w) in enumerate(results.memory_retention_factors)
        if ci_is_visible(results.fidelity_ci95[index, :], 0.9)
            errorbars!(
                axis,
                results.success_probabilities,
                results.mean_fidelities[index, :],
                results.fidelity_ci95[index, :];
                color=(colors[index], 0.65),
                whiskerwidth=7,
            )
        end
        scatterlines!(
            axis,
            results.success_probabilities,
            results.mean_fidelities[index, :];
            color=colors[index],
            marker=markers[index],
            markersize=9,
            linewidth=2,
            label="p_w = $p_w",
        )
    end

    axislegend(axis; position=:rb, nbanks=2)

    return save_figure(
        figure,
        output_directory,
        "mean_fidelity_vs_success_probability.png",
    )
end

"""
    run_waiting_time_vs_users_sweep(;
        k_values=2:20,
        p_success_values=(0.1, 0.3, 0.9),
        target_samples=200,
        burn_in_samples=100,
        simulation_chunk=100.0,
        repetitions=5,
        w=0.9,
        p_w=0.99,
        seed=9876,
    )

Measure mean waiting time for every combination of client count and success
probability. Each simulation advances in fixed time chunks until it has
collected the requested number of samples.

# Keyword arguments
- `k_values`: numbers of clients shown on the horizontal axis.
- `p_success_values`: success probabilities represented by separate curves.
- `target_samples::Int`: number of retained BSM samples in each run.
- `burn_in_samples::Int`: initial BSM samples excluded from each run.
- `simulation_chunk::Float64`: amount of simulated time added at each step
  while collecting samples.
- `repetitions::Int`: number of independent runs at each point.
- `w::Float64`: Werner parameter of every generated link-level pair.
- `p_w::Float64`: memory-retention factor per simulated time unit.
- `seed::Int`: base seed used to derive one deterministic seed per run.

# Returns
A named tuple containing the simulated waiting times, 95% confidence
half-widths, sample counts, mean final simulation times, and fixed benchmark
parameters.
"""
function run_waiting_time_vs_users_sweep(;
        k_values=2:20,
        p_success_values=(0.1, 0.3, 0.9),
        target_samples::Int=200,
        burn_in_samples::Int=100,
        simulation_chunk::Float64=100.0,
        repetitions::Int=5,
        w::Float64=0.9,
        p_w::Float64=0.99,
        seed::Int=9876,
    )
    client_counts = Int.(collect(k_values))
    success_probabilities = Float64.(collect(p_success_values))
    dimensions = (length(success_probabilities), length(client_counts))

    simulated_waiting_times = Matrix{Float64}(undef, dimensions)
    waiting_time_ci95 = similar(simulated_waiting_times)
    sample_counts = Matrix{Int}(undef, dimensions)
    simulation_times = Matrix{Float64}(undef, dimensions)

    for (probability_index, p_success) in enumerate(success_probabilities)
        for (client_index, number_of_users) in enumerate(client_counts)
            repetition_means = Float64[]
            repetition_times = Float64[]
            point_sample_count = 0

            for repetition in 1:repetitions
                sample = collect_sample_mean(
                    :waiting_times;
                    number_of_users,
                    w,
                    p_success,
                    p_w,
                    target_samples,
                    burn_in_samples,
                    simulation_chunk,
                    seed=parameter_seed(
                        seed,
                        probability_index,
                        client_index,
                        repetition,
                    ),
                )
                push!(repetition_means, sample.mean_value)
                push!(repetition_times, sample.simulation_time)
                point_sample_count += sample.sample_count
            end

            mean_waiting_time, waiting_time_half_width =
                mean_and_ci95(repetition_means)
            simulated_waiting_times[probability_index, client_index] =
                mean_waiting_time
            waiting_time_ci95[probability_index, client_index] =
                waiting_time_half_width
            sample_counts[probability_index, client_index] = point_sample_count
            simulation_times[probability_index, client_index] =
                mean(repetition_times)
        end

        println("Completed waiting-time curve for p_success=$p_success")
    end

    return (;
        client_counts,
        success_probabilities,
        simulated_waiting_times,
        waiting_time_ci95,
        sample_counts,
        simulation_times,
        target_samples,
        burn_in_samples,
        simulation_chunk,
        repetitions,
        w,
        p_w,
    )
end

"""
    plot_waiting_time_vs_users(results; output_directory)

Plot the simulated mean waiting time against the number of clients, using a
separate curve for every link-generation success probability.

# Arguments
- `results`: output returned by `run_waiting_time_vs_users_sweep`.

# Keyword arguments
- `output_directory`: directory in which the PNG figure is saved.

# Returns
The path of the generated PNG file.
"""
function plot_waiting_time_vs_users(results; output_directory)
    figure = Figure(size=(980, 600))
    axis = Axis(
        figure[1, 1];
        title="Waiting time versus number of clients\n" *
            "w=$(results.w), p_w=$(results.p_w). 95% pointwise CI: " *
            "$(results.repetitions) runs × $(results.target_samples) samples, " *
            "burn-in $(results.burn_in_samples)/run",
        xlabel="Number of clients k",
        ylabel="Mean waiting time [time units]",
        xticks=2:2:20,
    )

    colors = cgrad(
        :viridis,
        length(results.success_probabilities);
        categorical=true,
    )
    waiting_time_span = maximum(results.simulated_waiting_times) -
        minimum(results.simulated_waiting_times)

    for (index, p_success) in enumerate(results.success_probabilities)
        color = colors[index]
        if ci_is_visible(results.waiting_time_ci95[index, :], waiting_time_span)
            errorbars!(
                axis,
                results.client_counts,
                results.simulated_waiting_times[index, :],
                results.waiting_time_ci95[index, :];
                color=(color, 0.65),
                whiskerwidth=6,
            )
        end
        scatterlines!(
            axis,
            results.client_counts,
            results.simulated_waiting_times[index, :];
            color,
            marker=:circle,
            markersize=7,
            linewidth=2,
            label="p = $p_success",
        )
    end

    axislegend(axis; position=:rt, nbanks=2)

    return save_figure(
        figure,
        output_directory,
        "waiting_time_vs_number_of_users.png",
    )
end


function print_benchmark_analysis()
    println("\n================ BENCHMARK ANALYSIS ================")
    println(
        "\nVisible error bars are pointwise 95% Student-t confidence " *
        "intervals across five independent runs. Error bars are omitted for " *
        "a curve when all its intervals are smaller than 0.5% of the " *
        "vertical scale used for that plot.",
    )

    println("\n1. Mean fidelity versus link success probability")
    println(
        "Increasing the link success probability makes the mean fidelity " *
        "approach the value expected without memory degradation, independently " *
        "of p_w. A successfully generated Bell pair then waits less time for a " *
        "second pair, so memory depolarization has less time to act before " *
        "entanglement swapping.",
    )
    println(
        "At p_success=1, all simulated curves converge to 0.8575, which " *
        "matches the theoretical value (1+3w^2)/4 for w=0.9.",
    )
    println(
        "For the strongest simulated memory noise, p_w=0.1, the fidelity at " *
        "p_success=0.05 is 0.2981 ± 0.0354, close to the completely " *
        "depolarized limit " *
        "1/4=0.25" *
        ". As p_w approaches one, depolarization becomes negligible and the " *
        "dependence on p_success weakens. For p_w=1, memory depolarization is " *
        "absent and the fidelity remains equal to 0.8575 throughout the sweep.",
    )

    println("\n2. Mean fidelity versus Werner parameter")
    println(
        "The simulated fidelity follows the theoretical quadratic law " *
        "F=(1+3w^2)/4. Across the 21 simulated w values, the maximum absolute " *
        "discrepancy is 3.3306690738754696e-16, which is only numerical " *
        "round-off. This comparison isolates link-state quality by setting " *
        "p_success=1 and p_w=1, so link generation is deterministic and memory " *
        "noise is absent.",
    )

    println("\n3. Throughput versus link success probability")
    println(
        "The rate of end-to-end entanglement generation increases with the " *
        "link success probability and approaches the ideal saturation value " *
        "k/2 pairs per time unit. At p_success=1, the measured endpoint values " *
        "are 1.0 for k=2, 2.0 for k=4, and 4.0 for k=8, exactly matching the " *
        "corresponding saturation limits.",
    )
    println(
        "Larger client populations therefore support a higher total rate, and " *
        "the curves show an approximately linear increase over much of the " *
        "tested probability range before reaching their saturation limits.",
    )

    println("\n4. Waiting time versus number of clients")
    println(
        "Increasing the number of clients raises the probability that a stored " *
        "link-level Bell pair quickly finds another pair for entanglement " *
        "swapping. At p_success=0.9, increasing k from 2 to 20 reduces the " *
        "measured mean waiting time from 0.208 ± 0.010 to 0.048 ± 0.014 " *
        "time units.",
    )
    println(
        "Increasing p_success also shortens the wait because link-level pairs " *
        "are generated more rapidly. For k=2, the measured mean waiting time " *
        "decreases from 9.829 ± 0.680 at p_success=0.1 to 0.208 ± 0.010 at " *
        "p_success=0.9.",
    )
    println(
        "The graph also suggests an odd-even pattern, particularly at high " *
        "success probabilities: waiting times appear higher for odd values of " *
        "k than for nearby even values. What happen is that for odd values of k " *
        "there are more chances for a link-level pair to not find another pair " *
        "for entanglement swapping, needing to wait more time. This is much more " *
        "evident for lower values of k."
    )

    return nothing
end

"""Run every non-ideal benchmark and save all figures."""
function main()
    output_directory = joinpath(
        @__DIR__,
        "..",
        "results",
        "nonideal_switch",
    )

    println("\n[1/4] Fidelity versus Werner parameter")
    fidelity_results = run_fidelity_vs_w_sweep()
    test_fidelity_vs_w_sweep(fidelity_results)
    fidelity_plot_path =
        plot_fidelity_vs_w_sweep(fidelity_results; output_directory)

    println("\n[2/4] Throughput versus success probability")
    throughput_results = run_throughput_vs_success_probability_sweep()
    throughput_plot_path = plot_throughput_vs_success_probability(
        throughput_results;
        output_directory,
    )

    println("\n[3/4] Fidelity versus success probability")
    fidelity_success_results = run_fidelity_vs_success_probability_sweep()
    fidelity_success_plot_path = plot_fidelity_vs_success_probability(
        fidelity_success_results;
        output_directory,
    )

    println("\n[4/4] Waiting time versus number of users")
    waiting_time_results = run_waiting_time_vs_users_sweep()
    waiting_time_plot_path = plot_waiting_time_vs_users(
        waiting_time_results;
        output_directory,
    )

    print_benchmark_analysis()

    return (;
        fidelity_results,
        fidelity_plot_path,
        throughput_results,
        throughput_plot_path,
        fidelity_success_results,
        fidelity_success_plot_path,
        waiting_time_results,
        waiting_time_plot_path,
    )
end


if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
