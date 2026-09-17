#!/usr/bin/env julia
#
# OPTIONAL: verify that the Julia `bandpass` in func.jl reproduces
# `astrochron::bandpass` (via CretaceousConstraints::bandpass_filter).
#
#     julia --project=crosscheck crosscheck_astrochron.jl
#
# The Julia filter in func.jl was written to mimic astrochron's, and the
# published durations use the Julia one. This script is the evidence that the
# two agree; it is NOT needed to reproduce any published number.
#
# Requires, beyond the main environment:
#   - R on PATH
#   - R packages: tidyverse, astrochron, and
#       remotes::install_github("japhir/CretaceousConstraints")

using Pkg
Pkg.activate(joinpath(@__DIR__, "crosscheck"))

using Downloads, ZipFile, CSV, Arrow, DataFrames, SHA
using Statistics, StatsBase, Distributions, Bootstrap
using DSP, Peaks
using CairoMakie, AlgebraOfGraphics
using RCall

cd(@__DIR__)
mkpath("imgs")

R"""
suppressPackageStartupMessages({
  library(tidyverse)
  if (!requireNamespace("CretaceousConstraints", quietly = TRUE)) {
    stop("Install it with: remotes::install_github('japhir/CretaceousConstraints')")
  }
  library(CretaceousConstraints)  # provides bandpass_filter() / taner_filter()
})
"""

include(joinpath(@__DIR__, "func.jl"))    # pure-Julia core
include(joinpath(@__DIR__, "func_R.jl"))  # astrochron wrappers

# --- one solution, one window -------------------------------------------------
solution, tmin, tmax = 1, -1205e3, -1200e3
x = get_ZB23(solution, tmin, tmax)

# ETP, exactly as in ecc_wrapper_julia (weights from Figure S12 caption)
etp_weights = [1.5, 1.2, 1.2]
x.etp .= etp_weights[1] .* zscore(x.ecc) .+
         etp_weights[2] .* zscore(x.epl) .+
         etp_weights[3] .* zscore(x.cp)

# the passbands the published numbers use: one source of truth (func.jl)
filter_freqs = sort(DataFrame(target = string.(keys(Milankovitch_targets)),
                              flow   = first.(values(Milankovitch_targets)),
                              fhigh  = last.(values(Milankovitch_targets))),
                    :target)

# --- R (astrochron) -----------------------------------------------------------
Rflt = rcopy(R"""
$(x) |> bandpass_filter(frequencies = $(filter_freqs), x = time, y = etp)
""")

# --- Julia (func.jl) ----------------------------------------------------------
jl = DataFrame()
for r in eachrow(filter_freqs)
    append!(jl, DataFrame(target = r.target,
                          time = x.time,
                          filter = bandpass(x.etp, 0.4, r.flow, r.fhigh)))
end

# --- compare ------------------------------------------------------------------
cmp = innerjoin(rename(select(Rflt, :target, :time, :filter), :filter => :R),
                rename(jl, :filter => :julia),
                on = [:target, :time])
# the join is on Float64 times; if R re-gridded them nothing would match
nrow(cmp) == nrow(jl) || error("join dropped rows: $(nrow(cmp)) of $(nrow(jl)) matched")

stats = combine(groupby(cmp, :target),
                [:R, :julia] => ((a, b) -> cor(a, b)) => :correlation,
                [:R, :julia] => ((a, b) -> sqrt(mean((a .- b) .^ 2))) => :rmse,
                :R => (a -> maximum(abs, a)) => :amplitude)
stats.rmse_pct = 100 .* stats.rmse ./ stats.amplitude

@info "astrochron vs. Julia bandpass, ZB23.R$(lpad(solution, 2, '0')), $(tmin/1e3) to $(tmax/1e3) Myr"
show(stdout, MIME("text/plain"), stats); println()
CSV.write("out/crosscheck_astrochron_vs_julia.csv", stats)

# record the R side of the comparison
versions = rcopy(R"""
c(R = R.version.string,
  astrochron = as.character(packageVersion("astrochron")),
  CretaceousConstraints = as.character(packageVersion("CretaceousConstraints")))
""")
open("out/crosscheck_versions.txt", "w") do io
    println(io, "Julia ", VERSION)
    foreach(p -> println(io, p.first, " ", p.second), pairs(versions))
end
@info "wrote out/crosscheck_versions.txt"

plt = AlgebraOfGraphics.data(Rflt) *
        mapping(:time, :filter, row = :target) *
        visual(Lines, color = :cyan, label = "astrochron") +
      AlgebraOfGraphics.data(jl) *
        mapping(:time, :filter, row = :target) *
        visual(Lines, color = :red, linestyle = :dash, label = "Julia")

f = draw(plt, facet = (; linkyaxes = :none))
save("imgs/crosscheck_astrochron_vs_julia.png", f)
@info "wrote imgs/crosscheck_astrochron_vs_julia.png"
