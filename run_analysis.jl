#!/usr/bin/env julia
#
# Reproduce the ZB23.N64 cycle-duration analysis end-to-end.
#
#     julia --project=. run_analysis.jl
#
# Pure Julia: no R, no private packages. Writes result tables to out/ and
# figures to imgs/. The ZB23.RXX astronomical solutions are downloaded from
# the University of Hawaii server on first use and cached as out/*.arrow;
# subsequent runs reuse the cache. A cold run downloads ~63 solutions.
#
# The optional astrochron cross-check lives in crosscheck_astrochron.jl.

using Pkg
Pkg.activate(@__DIR__)

using Downloads, ZipFile, CSV, Arrow, DataFrames, SHA
using Statistics, StatsBase, Distributions, Bootstrap
using DSP      # bandpass filtering
using Peaks    # peak detection
using CairoMakie, AlgebraOfGraphics

include(joinpath(@__DIR__, "func.jl"))

cd(@__DIR__)
mkpath("out"); mkpath("imgs")

# ZB23.R61 is not published on the server, so it is skipped throughout.
const SOLUTIONS = [1:60; 62:64]

# The two time windows reported in the paper (kyr, negative = before present).
const WINDOW_SMALL = (-1205e3, -1200e3)
const WINDOW_BIG   = (-1250e3, -1200e3)

"Run `ecc_wrapper_julia` over every solution and stack the results."
function all_solutions(tmin, tmax; all_peaks = false, solutions = SOLUTIONS)
    reduce(vcat,
           ecc_wrapper_julia(s, mean, tmin, tmax; all_peaks = all_peaks)
           for s in solutions)
end

"""
    summarize_durations(df; value = :avg_duration)

Mean cycle duration per target with a 95% interval. The interval is the
empirical 2.5th and 97.5th percentile of the values in `value` (one per
solution per 1 Myr block); no bootstrap is involved, so the pipeline is fully
deterministic.

Julia port of the `dplyr::summarize()` call in astronomical_boot.jl; R's
`quantile()` default (type 7) matches `Statistics.quantile`, so the mean and
the interval are identical to the R version.

`n_rows` is the number of values the mean is taken over (one per solution per
1 Myr block); `total_durations` is the number of peak-to-peak durations
underlying them.
"""
function summarize_durations(df; value = :avg_duration)
    counted = :n_durations in propertynames(df) ?
                  (:n_durations => sum => :total_durations) :
                  (nrow => :total_durations)
    out = combine(groupby(df, :variable),
                  value => mean => :mean,
                  value => (x -> quantile(x, 0.05 / 2)) => :lwr,
                  value => (x -> quantile(x, 1 - 0.05 / 2)) => :upr,
                  nrow => :n_rows,
                  counted)
    out.lwr_ci = out.mean .- out.lwr
    out.upr_ci = out.upr .- out.mean
    sort!(out, :variable)
    return out
end

# ---------------------------------------------------------------- durations --

@info "1 Myr block means, $(WINDOW_SMALL[1]/1e3) to $(WINDOW_SMALL[2]/1e3) Myr"
ms = all_solutions(WINDOW_SMALL...)
CSV.write("out/ZB23.N64_filtered_durations_-1205_-1200_.csv", ms)

@info "1 Myr block means, $(WINDOW_BIG[1]/1e3) to $(WINDOW_BIG[2]/1e3) Myr"
mes = all_solutions(WINDOW_BIG...)
CSV.write("out/ZB23.N64_filtered_durations_-1250_-1200_.csv", mes)

avgs     = summarize_durations(ms)
avgs_big = summarize_durations(mes)

CSV.write("out/summary_durations_-1205_-1200.csv", avgs)
CSV.write("out/summary_durations_-1250_-1200.csv", avgs_big)

@info "cycle durations, -1205 to -1200 Myr"
show(stdout, MIME("text/plain"), avgs)
println()
@info "cycle durations, -1250 to -1200 Myr"
show(stdout, MIME("text/plain"), avgs_big)
println()

# Duration ratios between targets, per solution and per 1 Myr block.
# (out/ratio.csv from the exploratory R path in astronomical_boot.jl is the
# older equivalent; this one comes from the Julia durations. Neither is shipped.)
ratios = unstack(ms, [:solution, :grp], :variable, :avg_duration)
CSV.write("out/ratios_julia_-1205_-1200.csv", ratios)

# ------------------------------------------------------------- all peaks ----
# Every individual peak-to-peak duration rather than the 1 Myr block means.
# Large (the -1250 file is ~63 MB), so this is opt-in:
#     julia --project=. run_analysis.jl --all-peaks

if "--all-peaks" in ARGS
    @info "all individual peak durations (this is slow and writes ~70 MB)"
    meall = all_solutions(WINDOW_SMALL...; all_peaks = true)
    CSV.write("out/ZB23.N64_filtered_all_duration_-1205_-1200.csv", meall)
    mall = all_solutions(WINDOW_BIG...; all_peaks = true)
    CSV.write("out/ZB23.N64_filtered_all_duration_-1250_-1200.csv", mall)

    avgs_all     = summarize_durations(meall; value = :duration)
    avgs_all_big = summarize_durations(mall;  value = :duration)
    CSV.write("out/summary_all_durations_-1205_-1200.csv", avgs_all)
    CSV.write("out/summary_all_durations_-1250_-1200.csv", avgs_all_big)
end

# ------------------------------------------------------------------ figures --

# Both figures stack the 9 targets as facet rows and the targets have wildly
# different y-scales (E ~405 kyr down to p2 ~13 kyr), so linkyaxes = :none is
# required *and* the figures have to be tall: at Makie's default 1200x900 nine
# rows collapse into unreadable slivers.

# Duration through time, per solution, with per-1-Myr violins.
# `solution` is mapped as a number so the 63 ensemble members render as one
# compact colorbar instead of a 63-entry legend that crowds out the panels.
mes_fig = copy(mes)
mes_fig.solnum = [parse(Int, match(r"R0*(\d+)", s)[1]) for s in mes_fig.solution]

dt = AlgebraOfGraphics.data(mes_fig)
mp = mapping(:avg_time => (x -> x / 1e3) => "Time (Myr)",
             :avg_duration => "Peak duration (kyr)",
             color = :solnum => "ZB23.R## solution",
             group = :solnum => nonnumeric,
             row = :variable)

plt_dur = dt * mp * visual(Lines, alpha = 0.35, linewidth = 0.8) +
    dt * mapping(:grp => (x -> (x .+ 0.5)) => "Time (Myr)",
                 :avg_duration => "Peak duration (kyr)",
                 group = :grp => nonnumeric,
                 row = :variable) *
        visual(Violin, color = (:grey25, 0.45), width = 0.9) +
    mapping([-1200, -1205, -1250]) *
        visual(VLines, color = :black, linewidth = 1.2, linestyle = :dash)

f = draw(plt_dur,
         # continuous colour scales take `colormap`, not `palette`
         scales(Color = (; colormap = :viridis)),
         facet = (; linkyaxes = :none),
         axis = (; yticks = WilkinsonTicks(5),
                 xticks = (WINDOW_BIG[1] / 1e3):5:(WINDOW_BIG[2] / 1e3)),
         figure = (; size = (1500, 1900)))
save("imgs/duration_vs_time.png", f)
@info "wrote imgs/duration_vs_time.png"

# Duration densities for both windows, with mean and 95% interval.
sml_dens = AlgebraOfGraphics.data(ms) *
    mapping(:avg_duration, color = direct("−1200 to −1205 Myr"), row = :variable) *
    visual(Density) +
    AlgebraOfGraphics.data(avgs) *
    mapping([:mean, :lwr, :upr], color = direct("−1200 to −1205 Myr"), row = :variable) *
    visual(VLines, linewidth = 3)

big_dens = AlgebraOfGraphics.data(mes) *
    mapping(:avg_duration, color = direct("−1200 to −1250 Myr"), row = :variable) *
    visual(Density) +
    AlgebraOfGraphics.data(avgs_big) *
    mapping([:mean, :lwr, :upr], color = direct("−1200 to −1250 Myr"), row = :variable) *
    visual(VLines, linewidth = 3)

f, ax = draw(big_dens + sml_dens,
             scales(Color = (; palette = [(:purple, 0.4), (:orange, 0.4)])),
             facet = (; linkxaxes = :none, linkyaxes = :none),
             axis = (; yticks = WilkinsonTicks(3)),
             figure = (; size = (1200, 1900)),
             legend = (; position = :top))
# AlgebraOfGraphics sets xlabelvisible/ylabelvisible = false on every facet, so
# assigning the label text alone (the old `f.content[9].xlabel = ...`) drew
# nothing; visibility has to be switched back on as well.
axs = [c for c in f.content if c isa Makie.Axis]
axs[end].xlabel = "Cycle duration (kyr)"
axs[end].xlabelvisible = true
axs[cld(length(axs), 2)].ylabel = "Density"
axs[cld(length(axs), 2)].ylabelvisible = true
save("imgs/cycle_duration_means_both.png", f)
@info "wrote imgs/cycle_duration_means_both.png"

@info "done"
