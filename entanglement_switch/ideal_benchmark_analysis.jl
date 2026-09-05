include(joinpath(@__DIR__, "IdealEntanglementSwitch.jl"))
using .IdealEntanglementSwitch

using CairoMakie
using Distributions: TDist, quantile
using Statistics

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

"""Run independent ideal simulations for every client count in `k_values`."""
function run_ideal_sweep(;
        k_values=2:10,
        simulation_time::Float64=101.0,
        repetitions::Int=5,
        seed::Int=1234,
    )
    client_counts = collect(k_values)

    throughputs = Float64[]
    throughput_ci95 = Float64[]
    mean_fidelities = Float64[]
    mean_fidelity_ci95 = Float64[]
    expected_throughputs = Int[]

    for (client_index, k) in enumerate(client_counts)
        repetition_throughputs = Float64[]
        repetition_fidelities = Float64[]

        for repetition in 1:repetitions
            result = run_ideal_switch(;
                number_of_users=k,
                simulation_time,
                seed=seed + 100 * client_index + repetition,
            )
            push!(repetition_throughputs, result.throughput)
            push!(repetition_fidelities, result.mean_fidelity)
        end

        throughput, throughput_half_width =
            mean_and_ci95(repetition_throughputs)
        mean_fidelity, fidelity_half_width =
            mean_and_ci95(repetition_fidelities)

        push!(throughputs, throughput)
        push!(throughput_ci95, throughput_half_width)
        push!(mean_fidelities, mean_fidelity)
        push!(mean_fidelity_ci95, fidelity_half_width)
        push!(expected_throughputs, fld(k, 2))

        println(
            "k=$k: throughput=$throughput ± $throughput_half_width, " *
            "mean fidelity=$mean_fidelity ± $fidelity_half_width " *
            "(95% CI)",
        )
    end

    return (;
        client_counts,
        throughputs,
        throughput_ci95,
        mean_fidelities,
        mean_fidelity_ci95,
        expected_throughputs,
        simulation_time,
        repetitions,
    )
end

"""Validate the two theoretical predictions for the ideal benchmark."""
function test_ideal_sweep(results; atol=1.0e-12)
    @assert all(
        isapprox.(
            results.throughputs,
            results.expected_throughputs;
            atol,
        ),
    )
    @assert all(fidelity -> isapprox(fidelity, 1.0; atol), results.mean_fidelities)
    @assert all(width -> isapprox(width, 0.0; atol), results.throughput_ci95)
    @assert all(width -> isapprox(width, 0.0; atol), results.mean_fidelity_ci95)
    return nothing
end

"""Generate the throughput and mean-fidelity plots as two separate PNG files."""
function plot_ideal_sweep(results; output_directory)
    mkpath(output_directory)

    throughput_figure = Figure(size=(800, 500))
    throughput_axis = Axis(
        throughput_figure[1, 1];
        title="Ideal switch: throughput versus number of clients\n" *
            "Pointwise 95% CI from $(results.repetitions) independent runs",
        xlabel="Number of clients k",
        ylabel="Throughput [end-to-end pairs / time unit]",
        xticks=results.client_counts,
    )
    throughput_span = max(
        maximum(results.throughputs) - minimum(results.throughputs),
        1.0,
    )
    if ci_is_visible(results.throughput_ci95, throughput_span)
        errorbars!(
            throughput_axis,
            results.client_counts,
            results.throughputs,
            results.throughput_ci95;
            color=(:dodgerblue, 0.65),
            whiskerwidth=8,
        )
    end
    scatterlines!(
        throughput_axis,
        results.client_counts,
        results.throughputs;
        color=:dodgerblue,
        marker=:circle,
        markersize=12,
        linewidth=2,
    )

    throughput_path = joinpath(output_directory, "throughput_vs_clients.png")
    save(throughput_path, throughput_figure)

    fidelity_figure = Figure(size=(800, 500))
    fidelity_axis = Axis(
        fidelity_figure[1, 1];
        title="Ideal switch: mean fidelity versus number of clients\n" *
            "Pointwise 95% CI from $(results.repetitions) independent runs",
        xlabel="Number of clients k",
        ylabel="Mean fidelity with respect to |Phi+>",
        xticks=results.client_counts,
        limits=(nothing, (0.98, 1.005)),
    )
    if ci_is_visible(results.mean_fidelity_ci95, 1.005 - 0.98)
        errorbars!(
            fidelity_axis,
            results.client_counts,
            results.mean_fidelities,
            results.mean_fidelity_ci95;
            color=(:darkorange, 0.65),
            whiskerwidth=8,
        )
    end
    scatterlines!(
        fidelity_axis,
        results.client_counts,
        results.mean_fidelities;
        color=:darkorange,
        marker=:circle,
        markersize=12,
        linewidth=2,
    )

    fidelity_path = joinpath(output_directory, "mean_fidelity_vs_clients.png")
    save(fidelity_path, fidelity_figure)

    return (; throughput_path, fidelity_path)
end

function main()
    results = run_ideal_sweep()
    test_ideal_sweep(results)

    output_directory = joinpath(@__DIR__, "..", "results", "ideal_switch")
    plot_paths = plot_ideal_sweep(results; output_directory)

    println("Ideal benchmark passed for k = $(results.client_counts).")
    println("Throughput plot: $(abspath(plot_paths.throughput_path))")
    println("Mean-fidelity plot: $(abspath(plot_paths.fidelity_path))")

    return (; results, plot_paths)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
