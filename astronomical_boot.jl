# EXPLORATORY RECORD. This is the original working script. The portability
# edits are marked "reproducibility fix" below; the density comparison of the
# all-peak durations (around `mall` / `meall`) was added after the original
# analysis. It is meant to be stepped through interactively, not run
# top-to-bottom: several blocks are mutually exclusive experiments and some
# reference variables whose defining block was later commented out or lost
# (`ml`, `cpks`, `bsi`, `md`, `ci`, `hpdi`, `sm`, `sstd`, `sci`, `megaboot`).
# The reproducible pipeline is run_analysis.jl.
#
# It uses RCall, so it needs the crosscheck/ environment (see README.md).

using Pkg
# reproducibility fix: was Pkg.activate("../../prj/2023-06-03_bayesian_snvec/"),
# which borrowed the environment of the unpublished SNVec.jl repo. Nothing here
# ever used SNVec itself -- only the packages that environment happened to list.
Pkg.activate(joinpath(@__DIR__, "crosscheck"))

# reading/writing data files
using Downloads
using ZipFile
using DataFrames
using CSV
using Arrow
using SHA # checksums of the ZB23 downloads

import MCMCChains: hpd, Chains
using StatsBase
using Distributions
import Bootstrap: bootstrap
using DSP # for bandpass filtering
using Peaks # to identify peaks

# reproducibility fix: was GLMakie, which needs a GPU and a display. CairoMakie
# renders the same figures headlessly. For interactive use (DataInspector etc.)
# run Pkg.add("GLMakie") and swap this back.
using CairoMakie
const GLMakie = CairoMakie  # the qualified GLMakie.density(...) calls below
using AlgebraOfGraphics

using RCall # for astrochron code
#R"library(astrochron)"
R"library(tidyverse)"
# reproducibility fix: was devtools::load_all() on a local path that no longer
# exists. The package is public; install it once with
#   remotes::install_github("japhir/CretaceousConstraints")
R"library(CretaceousConstraints)"

# reproducibility fix: func.jl was split so the pure-Julia core loads without
# RCall. func_R.jl holds the astrochron wrappers used further down. (Was
# Revise.includet; Revise is not part of the crosscheck environment, add it
# to your global environment if you want live reloading.)
include("func.jl")
include("func_R.jl")

# read full 3.5 Gyr files (~8.75 million rows each, hundreds of MB of download).
# Not needed for the analysis; only used for the obliquity overview plot below.
# dat1 = get_ZB23_full(1)
# dat2 = get_ZB23_full(2)
# dat3 = get_ZB23_full(3)

# d = vcat(dat1, dat2)

# using AlgebraOfGraphics
# pl = data(d) *
#     mapping(:time => (x -> x * 1e-3) => "Time (Myr)",
#             :epl => (x -> x * 180/pi) => "Obliquity (°)",
#             color = :sol => "") *
#                 visual(Lines, alpha = 0.6)
# draw(pl)

# limit to certain age range (in kyr)
#dat2 = get_ZB23(1, -337e3, -300e3)


# plot single solution
# fig = Figure()
# ax = Axis(fig[1,1], xlabel = "Time (kyr)", ylabel = "Eccentricity (-)")
# #lines!(ax, dat.time, dat.ecc, label = unique(dat.sol))
# #xlims!(ax, (-337e3, -300e3))
# lines!(ax, dat2.time, dat2.ecc, color = dat2.sol, label = dat2.sol)


# all of them
# 61 is missing!!!
# allsols = vcat(get_ZB23.(1:60, -337e3, -300e3))
# allsols = vcat(get_ZB23.(62:64, -337e3, -300e3))
# allsols = vcat(get_ZB23.(1:60, -337e3, -300e3))

# for Dan Lunt: every full 3.5 Gyr solution stacked (63 downloads, several GB)
# fullsols = reduce(vcat, get_ZB23_full.([1:60; 62:64]))



# read from Arrow cache
# nums = vcat(collect(1:60), collect(62:64))
# sols = repeat(["ZB23.R01"], length(nums))
# for (i,sl) in enumerate(nums)
#     sols[i] = "ZB23.R$(lpad(sl, 2, '0'))"
# end
# files = "/home/japhir/SurfDrive/Postdoc1/prj/2023-06-03_bayesian_snvec/out/" .* sols .* "_-337000.0--300000.0.arrow"

# out = DataFrame()
# for (i, fl) in enumerate(files)
#     out = vcat(out,  DataFrame(Arrow.Table(fl)))
# end
# Arrow.write("out/ZB23.R1-64_-337000.0--300000.0.arrow", out)

out = Arrow.Table("out/ZB23.R1-64_-337000.0--300000.0.arrow") |> DataFrame


# data(out) * mapping(:time, :cp, color = :sol, group = :sol) * visual(Lines, alpha = 0.1)
pl = data(out) * mapping(:time, :ecc, color = :sol, group = :sol) * visual(Lines, alpha = 0.1)

# draw(pl)

# spec = rcopy(R"""
# nested_spectral_analysis($(out), nest = "sol", method = "FFT", x = time, y = ecc)
# """)

#Arrow.write("out/ZB23.R1-64_-337000.0--300000.0_spectral.arrow", spec)

spec = Arrow.Table("out/ZB23.R1-64_-337000.0--300000.0_spectral.arrow") |> DataFrame




plspec = data(spec) * mapping(:frequency, :power, color = :sol, group = :sol) * visual(Lines)
draw(plspec)
DataInspector()



# flt = rcopy(R"""
# freqs = data.frame(target = "LEC", flow = 1/405 - 0.3 * 1/405, fhigh = 1/405 + 0.3 * 1/405)
# nested_taner_filter($(out), frequencies = freqs, nest = "sol", x = time, y = ecc)
# """)

# Arrow.write("out/ZB23.R1-64_-337000.0--300000.0_flt.arrow", flt)
flt = Arrow.Table("out/ZB23.R1-64_-337000.0--300000.0_flt.arrow") |> DataFrame

plflt = data(flt) * mapping(:time, :filter, group = :sol, color = :sol) * visual(Lines)

(pl + plflt) |> draw

plflt |> draw

pks = rcopy(R"""
peaks = $(flt) |>
  filter(sol == "ZB23.R01") |>
  select(time, filter) |>
  astrochron::peak() |>
  mutate(diff = Location - lag(Location))
""")

hist(disallowmissing(pks.diff[2:end]),
     bins = 40)


data(pks) * mapping(:diff) * visual(Histogram)







# reproduce what Margriet does
x = get_ZB23(1, -1205e3, -1200e3)
# ml = get_ZB23(2, -1205e3, -1200e3)

# calculate ETP
etp_weights = [1.5, 1.2, 1.2]
# uses StatsBase.zscore
x.necc .= zscore(x.ecc)
x.nobl .= zscore(x.epl)
x.nprec .= zscore(x.cp)
@. x.etp .= etp_weights[1] * x.necc + etp_weights[2] * x.nobl + etp_weights[3] * x.nprec

# first do it in R, because whe know it
filter_freqs = DataFrame(
    target = ["p1", "p2",  "p", "o1", "o2", "o12", "o3", "e", "E"
              ],
    flow =  [0.0653, 0.075, 0.065, 0.0415, 0.0472, 0.042, 0.057, 1/155, 1/440
                 ],
    fhigh = [0.07, 0.0785, 0.0785, 0.044, 0.0505, 0.0505, 0.06,# 0.0594, # OR 0.06 for some!
                 1/75, 1/370
                 ]
)

window_size = 1e3

Rflt = rcopy(R"""
cp_peaks = $(ml) |>
  bandpass_filter(frequencies = $(filter_freqs), x = time, y = etp) #|>
  # select(-flow, -fhigh) #|>
  # pivot_wider(names_from = target, values_from = filter)
  # astrochron::peak(plateau = FALSE) |>
  # mutate(diff = Location - lead(Location)) |>
  # drop_na()
""")


# filter
# rather than use astrochron::bandpass, let's see if we can use DSP
filtered_components = Dict()
for (name, (flow, fhigh)) in Milankovitch_targets
    filtered_components[name] = bandpass(ml.etp, 0.4, flow, fhigh)
end
# append to df
for (name, component) in filtered_components
    ml[!, Symbol(name)] = component
end
filt_long = stack(ml[!, Not([:ecc, :inc, :epl, :cp, :sol, :necc, :nobl, :nprec])],
                  Not([:time, :etp]),
                  variable_name = :target, value_name = :amplitude)

# find peaks
pks = rcopy(R"""
cp_peaks = $(Rflt) |>
  filter(target == "p1") |>
  select(time, filter) |>
  astrochron::peak(plateau = FALSE) |>
  mutate(diff = Location - lead(Location)) |>
  drop_na()
""")

pk = findmaxima(ml.p1)
jl_peaks = DataFrame(
    time = ml.time[pk.indices],
    peak = pk.heights
)

# make a figure
# fig = Figure()
# ax = Axis(fig[1, 1], xlabel="Time (kyr)", ylabel="Amplitude", title="Milankovitch Components")

rw = AlgebraOfGraphics.data(ml) * mapping(:time, :etp) * visual(Lines)
# rw |> draw

flt = AlgebraOfGraphics.data(Rflt) *
    mapping(:time, :filter, row = :target) * visual(Lines, color = :cyan, label = "R")
# flt |> draw

julflt = AlgebraOfGraphics.data(filt_long) *
    mapping(:time, :amplitude, row = :target) * visual(Lines, color = :red, label = "Julia")
# julflt |> draw


f, ax = rw + flt + julflt  |> draw

# add peaks
scatter!(f.content[8], pks.Location, pks.Peak_Value, label = "R", color = :cyan)
scatter!(f.content[8], jl_peaks.time, jl_peaks.peak, label = "Julia")


fig, ax, ln = lines(ml.time, ml.cp, label = ml.sol)
scatter!(ax, cpks.Location, cpks.Peak_Value)
scatter!(ax, cpks.Location[1:70], cpks.Peak_Value[1:70])
scatter!(ax, cpks.Location[71:140], cpks.Peak_Value[71:140])
scatter!(ax, cpks.Location[141:211], cpks.Peak_Value[141:211])
scatter!(ax, cpks.Location[212:282], cpks.Peak_Value[212:282])
scatter!(ax, cpks.Location[283:end], cpks.Peak_Value[283:end])

# calculate the mean and std
mn = mean(cpks.diff)
st = std(cpks.diff)
# std(cpks.diff; corrected = true) # equivalent to default
# std(cpks.diff; corrected = true, mean = mn) # doesn't change anything

# std(cpks.diff; corrected = false) # slightly lower value
# std(cpks.diff; corrected = false, mean = mn) # doesn't change anything
movingaverage(g, n) = [i < n ? mean(g[begin:i]) : mean(g[i-n+1:i]) for i in 1:length(g)]
cpks.rnmn = movingaverage(cpks.diff, 70)
rmn = mean(cpks.rnmn)
rst = std(cpks.rnmn)

# what margriet does
movingsum(g, n) = [i < n ? sum(g[begin:i]) : sum(g[i-n+1:i]) for i in 1:length(g)]
cpks.rsm = movingsum(cpks.diff, 70)
mmn = mean(cpks.rsm)/70
mst = std(cpks.rsm)/70



# bootstrapped mean
bs = bootstrap(mean, cpks.diff, BasicSampling(1e6))
bci = confint(bs, BasicConfInt(0.95)) |> Iterators.flatten |> collect

# bootstrapped standard deviation
#bsi = bootstrap(std, cpks.diff, BasicSampling(1e6))




# this just gets the actual peak diffs for each solution
mvs1 = my_very_simple_wrapper(1)
mvs2 = my_very_simple_wrapper(6)
mvs3 = my_very_simple_wrapper(11)
mvs4 = my_very_simple_wrapper(16)
mvs5 = my_very_simple_wrapper(21)
mvs6 = my_very_simple_wrapper(26)
mvs7 = my_very_simple_wrapper(31)
mvs8 = my_very_simple_wrapper(36)
mvs9 = my_very_simple_wrapper(41)
mvs10 = my_very_simple_wrapper(46)

# combine
mvs = vcat(mvs1, mvs2, mvs3, mvs4, mvs5, mvs6, mvs7, mvs8, mvs9, mvs10)

# calculate summary stats
vsm = mean(mvs)
vsstd = std(mvs)
vsstd = std(mvs; corrected = true, mean = vsm)
vsstd = std(mvs; corrected = false, mean = vsm)
vsmi = 14.132 # copied from excel vs. vsm = 14.1408
vsci = 0.3270688235 # copied from excel =confidence(0.05, vsstd, 75)

@info "$(vsm) ± $(vsstd) (mean ± σ)"
# 14.14 ± 1.544
@info "$(vsm) ± $(2vsstd) (mean ± 2σ)"

# NOTE!!! this does not take into account
# the limited total number of cycles correctly

# the normal way of calculating the 95% confidence interval
vsci = vsstd / sqrt(length(mvs)) * quantile(TDist(length(mvs)-1), 1-0.05/2)
# note T-distribution with pretty high N, so should be the same as from Normal
vsci = vsstd / sqrt(length(mvs)) * quantile(Normal(), 1-0.05/2)
# vsstd / sqrt(75) * quantile(Normal(), 1-0.05/2)
@info "$(vsm) ± $(vsci) (mean ± 95% CI)"
# 14.14 ± 0.051

# bootstrapped mean of all diffs
bmas = bootstrap(mean, mvs, BasicSampling(10_000))
bcis = confint(bmas, BasicConfInt(0.95))[1]
confint(bmas, PercentileConfInt(0.95))[1]
vci = [quantile(bmas.t1[1], 0.05/2), quantile(bmas.t1[1], 1-0.05/2)]
vhpdi = DataFrame(hpd(Chains(bmas.t1[1]), alpha = 0.05))[1, 2:3] |> collect
@info "$(bcis[1]) - $(bcis[3]-bcis[1]) +$(bcis[1]-bcis[2])"
# 14.14 - 0.0525 + 0.0517

# another way to calculate the bootstrapped mean
# that takes into account the fact that they are time series
# but is still caculating it for 5 Myr instead of 1 Myr
# so the uncertainty is too small
# it's a bit slow
bma = bootstrap(mean, mvs, MaximumEntropySampling(999))
bci = confint(bma, BCaConfInt(0.95))[1]
@info "$(bci[1]) - $(bci[3]-bci[1]) +$(bci[1]-bci[2])"
# 14.14 - 0.050 + 0.054

# all are very close to 14.14 ± 0.05 kyr!



# replicate moving block bootstrap
# https://en.wikipedia.org/wiki/Bootstrapping_(statistics)#Time_series:_Simple_block_bootstrap
# this allows us to create rolling subsets of length 75
# draw n/b = 5 of the rolling mean blocks, with replacement
# ss = rand(1:length(mvs1), 5) # euhh this is still not correct, it should do that for each little one
# # convert those 5 numbers to consecutive 75 runs (maybe I'm unlucky and it hits the end?)
# idx = [ s .+ collect(1:70) for s in ss ] |> Iterators.flatten |> collect
# # Then aligning these n/b blocks in the order they were picked, will give the bootstrap observations.
# boots1 = mvs1[idx]
# block_boot = bootstrap(mean, boots, BasicSampling(10_000))
# block_ci = confint(block_boot, BasicConfInt(0.95))[1]
# @info "$(block_ci[1]) - $(block_ci[3]-block_ci[1]) + $(block_ci[1]-block_ci[2])"
# => now integrated in my_simple_wrapper!

# moving block bootstrap
b = 69
N = 10_000
mb1 = my_simple_wrapper(1; b = b, N = N)
mb2 = my_simple_wrapper(6; b = b, N = N)
mb3 = my_simple_wrapper(11; b = b, N = N)
mb4 = my_simple_wrapper(16; b = b, N = N)
mb5 = my_simple_wrapper(21; b = b, N = N)
mb6 = my_simple_wrapper(26; b = b, N = N)
mb7 = my_simple_wrapper(31; b = b, N = N)
mb8 = my_simple_wrapper(36; b = b, N = N)
mb9 = my_simple_wrapper(41; b = b, N = N)
mb10 = my_simple_wrapper(46; b = b, N = N)

mb = vcat(mb1, mb2, mb3, mb4, mb5, mb6, mb7, mb8, mb9, mb10)

# summary stats
mb_m = mean(mb)
mb_s = std(mb)
# confidence intervals
#mbbci = sstd / sqrt(75) * quantile(TDist(75-1), 1-0.05/2)
# mbbci = sstd / sqrt(75) * quantile(TDist(75-1), 1-0.05/2)
mb_ci = [quantile(mb, 0.05/2), quantile(mb, 1-0.05/2)]
@info "$(mb_m) - $(mb_m - mb_ci[1]) + $(mb_ci[2]-mb_m) (mean ± 95% CI)"
# that's 14.14 +0.56 -0.37

# this is a "simple block bootstrap"
# get solution
# find peaks
# calc distances
# chop up into 5 non-overlapping blocks
# calculate bootstrapped mean for each
# concatenate them together
sb1 = my_wrapper(1, b = b, N = N) # 353 peaks
# sb1_2 = my_wrapper(2, b = b, N = N)
# sb1_3 = my_wrapper(3, b = b, N = N)
# sb1_4 = my_wrapper(4, b = b, N = N)
# sb1_5 = my_wrapper(5, b = b, N = N)
sb2 = my_wrapper(6, b = b, N = N)
# sb2_7 = my_wrapper(7, b = b, N = N)
# sb2_8 = my_wrapper(8, b = b, N = N)
# sb2_9 = my_wrapper(9, b = b, N = N)
# sb2_10 = my_wrapper(10, b = b, N = N)
sb3 = my_wrapper(11, b = b, N = N)
sb4 = my_wrapper(16, b = b, N = N)
sb5 = my_wrapper(21, b = b, N = N)
sb6 = my_wrapper(26, b = b, N = N)
sb7 = my_wrapper(31, b = b, N = N)
sb8 = my_wrapper(36, b = b, N = N)
sb9 = my_wrapper(41, b = b, N = N)
sb10 = my_wrapper(46, b = b, N = N)

sb = vcat(sb1, sb2, sb3, sb4, sb5, sb6, sb7, sb8, sb9, sb10)

# calculate some summary stats
sb_m = mean(sb)
sb_s = std(sb)
# sb_med = median(sb)
# mode(mb) == md
# sb_hpdi = DataFrame(hpd(Chains(mb)))[1, 2:3] |> collect
sb_ci = [quantile(sb, 0.05/2), quantile(sb, 1-0.05/2)]
@info "$(sb_m) -$(sb_m - sb_ci[1]) + $(sb_ci[2]-sb_m)"
# so that's 14.14 - 0.77 + 0.89


# what if we don't do any bootstrapping, just use the solutions as "samples"?
s1 = wrapper(1)
s_2 = wrapper(2)
s_3 = wrapper(3)
s_4 = wrapper(4)
s_5 = wrapper(5)
s2 = wrapper(6)
s_7 = wrapper(7)
s_8 = wrapper(8)
s_9 = wrapper(9)
s_10 = wrapper(10)
s3 = wrapper(11)
s_12 = wrapper(12)
s_13 = wrapper(13)
s_14 = wrapper(14)
s_15 = wrapper(15)
s4 = wrapper(16)
s_17 = wrapper(17)
s_18 = wrapper(18)
s_19 = wrapper(19)
s_20 = wrapper(20)
s5 = wrapper(21)
s_22 = wrapper(22)
s_23 = wrapper(23)
s_24 = wrapper(24)
s_25 = wrapper(25)
s6 = wrapper(26)
s_27 = wrapper(27)
s_28 = wrapper(28)
s_29 = wrapper(29)
s_30 = wrapper(30)
s7 = wrapper(31)
s_32 = wrapper(32)
s_33 = wrapper(33)
s_34 = wrapper(34)
s_35 = wrapper(35)
s8 = wrapper(36)
s_37 = wrapper(37)
s_38 = wrapper(38)
s_39 = wrapper(39)
s_40 = wrapper(40)
s9 = wrapper(41)
s_42 = wrapper(42)
s_43 = wrapper(43)
s_44 = wrapper(44)
s_45 = wrapper(45)
s10 = wrapper(46)
s_47 = wrapper(47)
s_48 = wrapper(48)
s_49 = wrapper(49)
s_50 = wrapper(50)
s_51 = wrapper(51)
s_52 = wrapper(52)
s_53 = wrapper(53)
s_54 = wrapper(54)
s_55 = wrapper(55)
s_56 = wrapper(56)
s_57 = wrapper(57)
s_58 = wrapper(58)
s_59 = wrapper(59)
s_60 = wrapper(60)
# s_61 = wrapper(61)
s_62 = wrapper(62)
s_63 = wrapper(63)
s_64 = wrapper(64)

# same as Margriet
# s = vcat(s1, s2, s3, s4, s5, s6, s7, s8, s9, s10)

# all 64 solutions for p and t, TODO e and E!
ms = vcat(
    s1, s_2, s_3, s_4, s_5, s2, s_7, s_8, s_9, s_10,
    s3, s_12, s_13, s_14, s_15, s4, s_17, s_18, s_19, s_20,
    s5, s_22, s_23, s_24, s_25, s6, s_27, s_28, s_29, s_30,
    s7, s_32, s_33, s_34, s_35, s8, s_37, s_38, s_39, s_40,
    s9, s_42, s_43, s_44, s_45, s10, s_47, s_48, s_49, s_50,
    s_51, s_52, s_53, s_54, s_55, s_56, s_57, s_58, s_59, s_60,
    # s_61, # this one doesn't exist on the website
    s_62, s_63, s_64
)

# Arrow.write("out/ms.arrow", ms)
# Arrow.write("out/ms_taner.arrow", ms)

ms = DataFrame(Arrow.Table("out/ms.arrow"))
# ms = DataFrame(Arrow.Table("out/ms_taner.arrow"))


bt = combine(groupby(ms, :target),
             :mean_diff_1Myr =>
                 x -> bootstrap(mean,
                                # identity, # this just turns it into a distribution with 10000 samples
                                x,
                                BasicSampling(10_000)))

rename!(bt, :mean_diff_1Myr_function => :boot)

bt.mean = zeros(nrow(bt))
bt.lwr = zeros(nrow(bt))
bt.upr = zeros(nrow(bt))
# bt[!, :smp] .= repeat([zeros(10_000)], nrow(bt))
for i in 1:nrow(bt)
    ci = confint(bt.boot[i], BasicConfInt(0.95))[1]
    bt.mean[i] = ci[1]
    bt.lwr[i] = ci[2]
    bt.upr[i] = ci[3]
    # bt.smp[i] .= bt.boot[i].t1[1]
end
bt.upr_ci = bt.upr .- bt.mean
bt.lwr_ci = bt.mean .- bt.lwr

# for copy paste
select(bt, :target, :mean, :lwr, :upr, :lwr_ci, :upr_ci)

smp = combine(groupby(bt, :target),
              :boot => ByRow(x -> x.t1[1]))
smp = combine(groupby(smp, :target), :boot_function => (x -> x[1]) => AsTable)


# ratio's?
ratio = rcopy(R"""
      $(ms) |>
   # select(-grp, -solution) |>
   pivot_wider(id_cols = c(solution, grp), names_from = target, values_from = c(# n,
mean_diff_1Myr))
""")

CSV.write("out/ratio.csv", ratio)


ecc_wrapper = ecc_wrapper_julia

# do the same for eccentricity
es1 = ecc_wrapper(1)
# CSV.write("out/es1.csv", es1)
# es1 = CSV.read("out/es1.csv", DataFrame)
# es1 = ecc_wrapper_julia(1) # re-implemented in Julia
es_2 = ecc_wrapper(2)
# CSV.write("out/es_2.csv", es_2)
es_3 = ecc_wrapper(3)
es_4 = ecc_wrapper(4)
es_5 = ecc_wrapper(5)
es2 = ecc_wrapper(6)
es_7 = ecc_wrapper(7)
es_8 = ecc_wrapper(8)
es_9 = ecc_wrapper(9)
es_10 = ecc_wrapper(10)
es3 = ecc_wrapper(11)
es_12 = ecc_wrapper(12)
es_13 = ecc_wrapper(13)
es_14 = ecc_wrapper(14)
es_15 = ecc_wrapper(15)
es4 = ecc_wrapper(16)
es_17 = ecc_wrapper(17)
es_18 = ecc_wrapper(18)
es_19 = ecc_wrapper(19)
es_20 = ecc_wrapper(20)
es5 = ecc_wrapper(21)
es_22 = ecc_wrapper(22)
es_23 = ecc_wrapper(23)
es_24 = ecc_wrapper(24)
es_25 = ecc_wrapper(25)
es6 = ecc_wrapper(26)
es_27 = ecc_wrapper(27)
es_28 = ecc_wrapper(28)
es_29 = ecc_wrapper(29)
es_30 = ecc_wrapper(30)
es7 = ecc_wrapper(31)
es_32 = ecc_wrapper(32)
es_33 = ecc_wrapper(33)
es_34 = ecc_wrapper(34)
es_35 = ecc_wrapper(35)
es8 = ecc_wrapper(36)
es_37 = ecc_wrapper(37)
es_38 = ecc_wrapper(38)
es_39 = ecc_wrapper(39)
es_40 = ecc_wrapper(40)
es9 = ecc_wrapper(41)
es_42 = ecc_wrapper(42)
es_43 = ecc_wrapper(43)
es_44 = ecc_wrapper(44)
es_45 = ecc_wrapper(45)
es10 = ecc_wrapper(46)
es_47 = ecc_wrapper(47)
es_48 = ecc_wrapper(48)
es_49 = ecc_wrapper(49)
es_50 = ecc_wrapper(50)
es_51 = ecc_wrapper(51)
es_52 = ecc_wrapper(52)
es_53 = ecc_wrapper(53)
es_54 = ecc_wrapper(54)
es_55 = ecc_wrapper(55)
es_56 = ecc_wrapper(56)
es_57 = ecc_wrapper(57)
es_58 = ecc_wrapper(58)
es_59 = ecc_wrapper(59)
es_60 = ecc_wrapper(60)
# es_61 = ecc_wrapper(61)
es_62 = ecc_wrapper(62)
es_63 = ecc_wrapper(63)
es_64 = ecc_wrapper(64)

# same as Margriet
# es = vcat(es1, es2, es3, es4, es5, es6, es7, es8, es9, es10)

# all 64 esolutiones for p and t, TODO e and E!
mes = vcat(
    es1, es_2, es_3, es_4,
    es_5, es2, es_7, es_8, es_9, es_10,
    es3, es_12, es_13, es_14, es_15,
    es4, es_17, es_18, es_19, es_20,
    es5, es_22, es_23, es_24, es_25,
    es6, es_27, es_28, es_29, es_30,
    es7, es_32, es_33, es_34, es_35,
    es8, es_37, es_38, es_39, es_40,
    es9, es_42, es_43, es_44, es_45,
    es10, es_47, es_48, es_49, es_50,
    es_51, es_52, es_53, es_54, es_55,
    es_56, es_57, es_58, es_59, es_60,
    # es_61, # this one doesn't exist on the website
    es_62, es_63, es_64
)

# save results
# ms = mes
# CSV.write("out/ZB23.N64_filtered_durations_-1205_-1200_.csv", ms)
# the one with _ has been revised to also calculate avg time and avg amplitude
# CSV.write("out/ZB23.N64_filtered_durations_-1250_-1200_.csv", mes)
# mall = mes
# CSV.write("out/ZB23.N64_filtered_all_duration_-1250_-1200.csv", mall)
# CSV.write("out/ZB23.N64_filtered_all_duration_-1205_-1200.csv", meall)

# ms = CSV.read("out/ZB23.N64_filtered_durations_-1205_-1200.csv", DataFrame)
# mes = CSV.read("out/ZB23.N64_filtered_durations_-1250_-1200.csv", DataFrame)
# fixed number of groups by calculating avg time of duration correctly
ms = CSV.read("out/ZB23.N64_filtered_durations_-1205_-1200_.csv", DataFrame)
mes = CSV.read("out/ZB23.N64_filtered_durations_-1250_-1200_.csv", DataFrame)

# combine(groupby(combine(groupby(mes, :variable), :avg_duration => unique), :variable), :avg_duration_unique => length => :n_unique)

# every individual peak-to-peak duration; written by
#     julia --project=. run_analysis.jl --all-peaks
mall = CSV.read("out/ZB23.N64_filtered_all_duration_-1250_-1200.csv", DataFrame)
meall = CSV.read("out/ZB23.N64_filtered_all_duration_-1205_-1200.csv", DataFrame)

avgs = rcopy(R"""
  $(ms) |>
   summarize(.by = c(variable), mean = mean(avg_duration),
             lwr = quantile(avg_duration, 0.05/2),
             upr = quantile(avg_duration, 1-0.05/2),
             lwr_ci = mean - lwr,
             upr_ci = upr - mean,
             total_durations = sum(n_durations),
             total_blocks = n())
""")


avgs_big = rcopy(R"""
  $(mes) |>
   summarize(.by = c(variable), mean = mean(avg_duration),
             lwr = quantile(avg_duration, 0.05/2),
             upr = quantile(avg_duration, 1-0.05/2),
             lwr_ci = mean - lwr,
             upr_ci = upr - mean,
             total_durations = sum(n_durations),
             total_blocks = n())
""")

avgs_all_big = rcopy(R"""
  $(mall) |>
   summarize(.by = c(variable), mean = mean(duration),
             lwr = quantile(duration, 0.05/2),
             upr = quantile(duration, 1-0.05/2),
             lwr_ci = mean - lwr,
             upr_ci = upr - mean,
             # total_durations = n(),
             total_peaks = n()) |>
   arrange(variable)
""")

avgs_all = rcopy(R"""
  $(meall) |>
   summarize(.by = c(variable), mean = mean(duration),
             lwr = quantile(duration, 0.05/2),
             upr = quantile(duration, 1-0.05/2),
             lwr_ci = mean - lwr,
             upr_ci = upr - mean,
             # total_durations = n(),
             total_peaks = n()) |>
   arrange(variable)
""")

dt = AlgebraOfGraphics.data(mes)
mp = mapping(:avg_time => (x -> x / 1e3) => "Time (Myr)",
                        :avg_duration => "Peak duration (kyr)",
                        color = :solution,
                        row = :variable)

plt_dur = dt * mp * # visual(Scatter) +
    visual(Lines, alpha = 0.2, legend = (;alpha = 1)) +
    dt *
    mapping(:grp => (x -> (x .+ 0.5)) => "Time (Myr)",
            :avg_duration => "Peak duration (kyr)",
            group = :grp => nonnumeric,
            row = :variable) *
                visual(Violin)  +
    mapping([-1200, -1205, -1250]) * visual(VLines)
f = draw(plt_dur, scales(Color = (; palette = from_continuous(:viridis))),
         facet = (;linkyaxes = :none))
save("imgs/duration_vs_time.png", f)
# the weird one with some peak durations up to 499.2 kyr is ZB23.R28

sml_dens = AlgebraOfGraphics.data(ms) *
    mapping(:avg_duration,
            color = direct("−1200 to −1205 Myr"),
            # color = :solution,
            # color = :grp,
            # stack = :solution,
            row = :variable) *
    visual(Density) +
    # and vlines for mean and 95% CI
    AlgebraOfGraphics.data(avgs) *
    mapping([:mean, :lwr, :upr],
            color = direct("−1200 to −1205 Myr"),
            row = :variable) *
    visual(VLines, linewidth = 3)

big_dens =  AlgebraOfGraphics.data(mes) *
   mapping(:avg_duration,
            color = direct("−1200 to −1250 Myr"),
            # color = :solution,
            # color = :grp,
            # stack = :solution,
            row = :variable) *
    visual(Density) +
    # and vlines for mean and 95% CI
    AlgebraOfGraphics.data(avgs_big) *
    mapping([:mean, :lwr, :upr],
            color = direct("−1200 to −1250 Myr"),
            row = :variable) *
    visual(VLines, linewidth = 3)

all_big_dens =  AlgebraOfGraphics.data(mall) *
   mapping(:duration,
            color = direct("−1200 to −1250 Myr all peaks"),
            # color = :solution,
            # color = :grp,
            # stack = :solution,
            row = :variable) *
    visual(Density, bandwidth = 0.2) +
    # and vlines for mean and 95% CI
    AlgebraOfGraphics.data(avgs_all_big) *
    mapping([:mean, :lwr, :upr],
            color = direct("−1200 to −1250 Myr all peaks"),
            row = :variable) *
    visual(VLines, linewidth = 3)

all_dens =  AlgebraOfGraphics.data(meall) *
   mapping(:duration,
            color = direct("−1200 to −1205 Myr all peaks"),
            # color = :solution,
            # color = :grp,
            # stack = :solution,
            row = :variable) *
    visual(Density, bandwidth = 0.4) +
    # and vlines for mean and 95% CI
    AlgebraOfGraphics.data(avgs_all) *
    mapping([:mean, :lwr, :upr],
            color = direct("−1200 to −1205 Myr all peaks"),
            row = :variable) *
    visual(VLines, linewidth = 3)


# plt = big_dens
draw(all_big_dens + big_dens, facet = (;linkxaxes = :none, linkyaxes = :none))
draw(all_dens + sml_dens, facet = (;linkxaxes = :none, linkyaxes = :none))


plt = # all_big_dens + all_dens +
    big_dens + sml_dens
f, ax = draw(plt, scales(Color = (; palette = [(:purple, 0.4), (:orange, 0.4)#,
                                               # (:cyan, 0.4), (:green, 0.4)
                                               ])),
             facet=(; linkxaxes=:none,linkyaxes=:none),
             legend = (; position = :top))
f.content[9].xlabel = "Cycle duration (kyr)"

# the same figure was saved under different names as the plot above was
# varied (taner filter, eccentricity-only, all peaks, ...). Only
# cycle_duration_means_both.png is regenerated by run_analysis.jl.
save("imgs/cycle_duration_means_both.png", f)
# save("imgs/cycle_duration.png", f)
# save("imgs/cycle_duration_taner.png", f)
# save("imgs/cycle_duration_means.png", f)
# save("imgs/cycle_duration_means_ecc.png", f)
# save("imgs/cycle_duration_means_all.png", f)

# for (i, trg) in enumerate(targets)
#     hist!(ax,
#           s.mean_diff_1Myr[s.target.==trg],
#           offset = i, direction = :x)
# end

# the ~1 Myr uncertainty is quite a bit larger
# than the if you calculate it for full 5 Myr (see above)


# plot
f = Figure()
ax = Axis(f[1,1], xlabel =  "ΔAge between peaks of cp (kyr)")
hist!(mvs, color = (:brown, .2), bins = 50,
     normalization = :pdf,
     label = "ZB23.N64 subset")
#hist!(mvs1, bins = 50, label = "ZB23.R01", normalization = :pdf)
#hist!(mvs2, label = "ZB23.R06", bins = 50, normalization = :pdf)
# hist!(mvs3, label = "ZB23.R11", normalization = :pdf)
# hist!(mvs4, label = "ZB23.R16", normalization = :pdf)
# hist!(mvs5, label = "ZB23.R21", normalization = :pdf)
# hist!(mvs6, label = "ZB23.R26", normalization = :pdf)
# hist!(mvs7, label = "ZB23.R31", normalization = :pdf)
# hist!(mvs8, label = "ZB23.R36", normalization = :pdf)
# hist!(mvs9, label = "ZB23.R41", normalization = :pdf)
# hist!(mvs10, label = "ZB23.R46", normalization = :pdf)

# hist!(mb1, label = "ZB23.R01", normalization = :pdf)
# hist!(mb2, label = "ZB23.R06", normalization = :pdf)
# hist!(mb3, label = "ZB23.R11", normalization = :pdf)
# hist!(mb4, label = "ZB23.R16", normalization = :pdf)
# hist!(mb5, label = "ZB23.R21", normalization = :pdf)
# hist!(mb6, label = "ZB23.R26", normalization = :pdf)
# hist!(mb7, label = "ZB23.R31", normalization = :pdf)
# hist!(mb8, label = "ZB23.R36", normalization = :pdf)
# hist!(mb9, label = "ZB23.R41", normalization = :pdf)
# hist!(mb10, label = "ZB23.R46", normalization = :pdf)

hist!(sb, normalization = :pdf, label = "1 Myr simple block")
hist!(mb, normalization = :pdf, label = "1 Myr moving block")
vlines!([sb_m, sb_ci[1], sb_ci[2]], label = "1 Myr simple block")
vlines!([mb_m, mb_ci[1], mb_ci[2]], label = "1 Myr moving block")

# bootstrapped full thing
hist!(bmas.t1[1], normalization = :pdf, label = "5 Myr bootstrapped mean")
vlines!(bcis |> collect, label = "5 Myr bootstrapped mean")


# very simple means
vlines!([vsm,
         vsm + vsstd, vsm - vsstd,
         vsm + 2vsstd, vsm - 2vsstd,
         vsm + vsci, vsm - vsci],
        color = [(:red, 1),
                 (:red, 0.6), (:red, 0.6),
                 (:red, 0.4), (:red, 0.4),
                 (:brown, .1), (:brown, .1)], label = "very simple means")

# our summary stats
# the ones we like?
vlines!([mn, md, ci[1], ci[2], hpdi[1], hpdi[2]],
        # color = [:red, :orange, :blue, :blue, :cyan, :cyan],
        label = ["mean", "median", "95% CI", "95% CI", "95% HPDI", "95% HPDI"])


# simple means
vlines!([sm,
         sm - sstd, sm + sstd,
         sm - 2sstd, sm + 2sstd,
         sm - sci, sm + sci],
        # color = [:purple, :purple, :purple,
        #          (:purple, 0.1), (:purple, 0.1),
        #          (:purple , .5), (:purple , .5)],
        label = ["mean", "-σ", "σ", "-2σ", "2σ", "95% lwr", "95% upr"])

# bootstrapped mean for all
vlines!(bci[1] |> collect, color = :yellow, linewidth = 5)
vlines!(bcis[1] |> collect, color = :orange, linewidth = 5)


# margriet's means etc.
mm = 14.130
mlsd = 0.182
vlines!([mm,
         mm - mlsd, mm + mlsd,
         mm - 2*mlsd, mm + 2*mlsd],
        # color = [:green, :darkgreen, :darkgreen, :darkgreen, :darkgreen],
        label = "Margriet")

axislegend()


bsm = bootstrap(mean, cpks.diff, MaximumEntropySampling(N))

bsm1 = bootstrap(mean, cpks.diff[1:70], MaximumEntropySampling(N)) # 13.88
bsm2 = bootstrap(mean, cpks.diff[71:140], MaximumEntropySampling(N)) # 14.1429
bsm3 = bootstrap(mean, cpks.diff[141:211], MaximumEntropySampling(N)) # 14.107
bsm4 = bootstrap(mean, cpks.diff[212:282], MaximumEntropySampling(N)) # 14.5746
bsm5 = bootstrap(mean, cpks.diff[283:end], MaximumEntropySampling(N)) # 14.5746

bcim = confint(bsm, PercentileConfInt(0.95)) |> Iterators.flatten |> collect

GLMakie.density(cpks.rnmn, color = :black)
GLMakie.density!(megaboot, color = (:pink, 0.6))
# incorrect! too large of an N
GLMakie.density!(bsm.t1 |> Iterators.flatten |> collect, color = :gray)
GLMakie.density!(bsm1.t1 |> Iterators.flatten |> collect)
GLMakie.density!(bsm2.t1 |> Iterators.flatten |> collect)
GLMakie.density!(bsm3.t1 |> Iterators.flatten |> collect)
GLMakie.density!(bsm4.t1 |> Iterators.flatten |> collect)
GLMakie.density!(bsm5.t1 |> Iterators.flatten |> collect)


x = vcat(
bsm1.t1 |> Iterators.flatten |> collect,
bsm2.t1 |> Iterators.flatten |> collect,
bsm3.t1 |> Iterators.flatten |> collect,
bsm4.t1 |> Iterators.flatten |> collect,
bsm5.t1 |> Iterators.flatten |> collect)

GLMakie.density!(x, color = :gray)

# # bootstrapped 70 cycle mean
# function bmi(x, N::Int = 70)
#     n = length(x)
#     rand(1:length(x))

# end

# plot distribution and mean + std
#hist(cpks.diff, bins = 50)
fig, ax = GLMakie.density(cpks.diff, label = "density of raw values")
density!(cpks.rnmn) # this is what Margriet gets
scatter!(ax, cpks.diff, repeat([0], length(cpks.diff)), color = :black, label = "raw values")
hist!(ax, cpks.diff, bins = 70; normalization = :pdf, label = "hist of raw values")
vlines!(ax, [mn - 2 * st, mn - st, mn, mn + st, mn + 2 * st],
        color = -2:2,
        label = ["-2σ", "-σ", "µ", "σ", "2σ"], colormap = [:red, :white, :red])
vlines!(ax, bci, color = [:red, :orange, :blue], label = ["bootstrapped median", "bootstrapped lower", "bootstrapped upper"])
#hist!(ax, bs.t1 |> Iterators.flatten |> collect; normalization = :pdf, label = "density of bootstrapped median")
density!(ax, bs.t1 |> Iterators.flatten |> collect, label = "density of bootstrapped mean")
density!(ax, bsm.t1 |> Iterators.flatten |> collect, label = "density of smart bootstrapped mean")
density!(ax, mn .+ 2*(bsi.t1 |> Iterators.flatten |> collect), label = "density of bootstrapped standard deviation")
density!(ax, mn .- 2*(bsi.t1 |> Iterators.flatten |> collect), label = "density of bootstrapped standard deviation")
vlines!(ax, [rmn - 2*rst, rmn - rst, rmn, rmn + rst, rmn + 2*rst], color = -2:2, colormap = [:red, :white, :red])
#axislegend()


# same for obliquity

opks = rcopy(R"""
obl_peaks = $(ml) |>
  select(time, epl) |>
  astrochron::peak(plateau = FALSE) |>
  mutate(diff = Location - lead(Location)) |>
  drop_na()
""")

mn = mean(opks.diff)
st = std(opks.diff)
bs = bootstrap(mean, opks.diff, BasicSampling(1e6))
bci = confint(bs, BasicConfInt(0.95)) |> Iterators.flatten |> collect


fig, ax, ln = lines(ml.time, ml.epl, label = ml.sol)
scatter!(ax, opks.Location, opks.Peak_Value)


# plot distribution and mean + std
#hist(cpks.diff, bins = 50)
fig, ax = GLMakie.density(opks.diff, label = "density of raw values")
scatter!(ax, opks.diff, repeat([0], length(opks.diff)), color = :black, label = "raw values")
hist!(ax, opks.diff, bins = 70; normalization = :pdf, label = "hist of raw values")
vlines!(ax, [mn - 2 * st, mn - st, mn, mn + st, mn + 2 * st],
        color = -2:2,
        label = ["-2σ", "-σ", "µ", "σ", "2σ"], colormap = [:red, :white, :red])

vlines!(ax, bci, color = [:red, :orange, :gray], label = ["bootstrapped median", "bootstrapped lower", "bootstrapped upper"])
#hist!(ax, bs.t1 |> Iterators.flatten |> collect; normalization = :pdf, label = "density of bootstrapped median")
density!(ax, bs.t1 |> Iterators.flatten |> collect, label = "density of bootstrapped mean")
density!(ax, mn .+ 2*(bsi.t1 |> Iterators.flatten |> collect), label = "density of bootstrapped standard deviation")
density!(ax, mn .- 2*(bsi.t1 |> Iterators.flatten |> collect), label = "density of bootstrapped standard deviation")
#axislegend()
