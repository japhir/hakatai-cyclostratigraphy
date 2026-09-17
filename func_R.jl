# OPTIONAL: astrochron-based cross-check wrappers.
#
# The published durations do NOT come from this file -- they come from
# `ecc_wrapper_julia` in func.jl. These wrappers are kept because the Julia
# `bandpass` in func.jl was written to mimic `astrochron::bandpass`, and they
# let you verify that the two give equivalent results (see
# crosscheck_astrochron.jl).
#
# Requires, on top of func.jl:
#   - RCall.jl (use the `crosscheck/` environment, which adds it)
#   - R, with: tidyverse, astrochron, and
#     remotes::install_github("japhir/CretaceousConstraints")
#       -> provides `bandpass_filter()` / `taner_filter()`
#
# `include("func.jl")` must have run first.


"gets the ZB23.RXX solution between tmin and tmax, extracts parameter, finds peaks, calculates differences"
my_very_simple_wrapper = function(solution = 1, tmin = -1205e3, tmax = -1200e3;
                                  b = 75, N = 999, parameter = "cp")
    x = get_ZB23(solution, tmin, tmax)
    pks = rcopy(R"""
cp_peaks = $(x) |>
  dplyr::select(all_of(c("time", $(parameter)))) |>
  astrochron::peak(plateau = FALSE, genplot = FALSE) |>
  dplyr::mutate(diff = Location - lead(Location)) |>
  tidyr::drop_na()
""")
    @info "identied peaks"

    return pks.diff
end

"rolling block bootstrap"
my_simple_wrapper = function(solution = 1, tmin = -1205e3, tmax = -1200e3;
                             b = 75, N = 999, parameter = "cp")
    x = get_ZB23(solution, tmin, tmax)
    pks = rcopy(R"""
cp_peaks = $(x) |>
  dplyr::select(all_of(c("time", $(parameter)))) |>
  astrochron::peak(plateau = FALSE, genplot = FALSE) |>
  dplyr::mutate(diff = Location - lead(Location)) |>
  tidyr::drop_na()
""")
    @info "identied peaks"
    megaboot = moving_block_bootstrap( pks.diff, mean; b = b, N = N )
    @info "booted"

    return megaboot
end

my_wrapper = function(solution = 1, fun = mean, tmin = -1205e3, tmax = -1200e3;
                      b = 75, N = 999, parameter = "cp")
    x = get_ZB23(solution, tmin, tmax)
    pks = rcopy(R"""
cp_peaks = $(x) |>
  dplyr::select(all_of(c("time", $(parameter)))) |>
  astrochron::peak(plateau = FALSE, genplot = FALSE) |>
  dplyr::mutate(diff = Location - lead(Location)) |>
  tidyr::drop_na()
""")
    @info "identied peaks"
    bt = my_boot(pks.diff, fun; b = b, N = N)
    @info "booted"
    return bt
end

wrapper = function(solution = 1, fun = mean, tmin = -1205e3, tmax = -1200e3)
    x = get_ZB23(solution, tmin, tmax)

    # create ETP with weights
    # 1.5, 1.2, 1.2 from Figure S12 caption
    etp_weights = [1.5, 1.2, 1.2]
    # uses StatsBase.zscore
    necc = zscore(x.ecc)
    nobl = zscore(x.epl)
    nprec = zscore(x.cp)
    @. x.etp .= etp_weights[1] * necc + etp_weights[2] * nobl + etp_weights[3] * nprec
    @info "calculated ETP with weights $etp_weights"

    filter_freqs = DataFrame(
    target = ["p1", "p2",  "p", "o1", "o2", "o12", "o3", "e", "E"],
        flow =  [0.0653, 0.075, 0.065, 0.0415, 0.0472, 0.042, 0.057, 1/155, 1/440],
        # NB: o3 upper bound 0.06 here (exploratory R runs); the published
        # pipeline in func.jl uses 0.0594, see Milankovitch_targets.
        fhigh = [0.07, 0.0785, 0.0785, 0.044, 0.0505, 0.0505, 0.06,
                 1/75, 1/370
                 ]
    )

    window_size = 1e3
    cut_groups = tmin:window_size:tmax |> collect
    cut_groups[end] = cut_groups[end] + 1

R"""
my_peak <- function(data) {
    data |>
        dplyr::select(all_of(c("time", "filter"))) |>
        astrochron::peak(plateau = FALSE, genplot = FALSE, verbose = FALSE) |>
        dplyr::rename(peak_age = Location, peak_value = Peak_Value) |>
        dplyr::mutate(diff = peak_age - dplyr::lag(peak_age)) |>
        tidyr::drop_na()
}
"""

    pks = rcopy(R"""
cp_flt = $(x) |>
  # this is from my R package CretaceousConstraints
  bandpass_filter(frequencies = $(filter_freqs), x = time, y = etp) |>
  # taner_filter(frequencies = $(filter_freqs), x = time, y = etp, roll = 1e3, verbose = FALSE, genplot = FALSE) |>
  tidyr::nest(.by = "target") |>
  dplyr::mutate(pk = purrr::map(data, my_peak)) |>
  dplyr::select(-"data") |>
  tidyr::unnest(cols = "pk") |>
  mutate(grp = cut(peak_age, $(cut_groups)))
""")
    @info "filtered frequencies"
    @info "identied peaks"

    means = rcopy(R"""
cp_flt |>
   summarize(.by = c(grp, target),
             mean_age = mean(peak_age),
             mean_diff_1Myr = mean(diff),
             n = n()) |>
   mutate(solution = $(x.sol[1]))
""")
    # TODO: do something else so it become
    # for (i, gr), (j, target) in enumerate(pks.grp), enumerate(pks.target)
    #     means = bootstrap(mean, pks.diff, MaximumEntropySampling(999))
    # end

    @info "calculated means"
    return means
end

ecc_wrapper = function(solution = 1, fun = mean, tmin = -1250e3, tmax = -1200e3)
    x = get_ZB23(solution, tmin, tmax)

    # create ETP with weights
    # 1.5, 1.2, 1.2 from Figure S12 caption
    etp_weights = [1.5, 1.2, 1.2]
    # uses StatsBase.zscore
    necc = zscore(x.ecc)
    nobl = zscore(x.epl)
    nprec = zscore(x.cp)
    @. x.etp .= etp_weights[1] * necc + etp_weights[2] * nobl + etp_weights[3] * nprec
    @info "calculated ETP with weights $etp_weights"

    filter_freqs = DataFrame(
    target = ["p1", "p2",  "p", "o1", "o2", "o12", "o3", "e", "E"],
        flow =  [0.0653, 0.075, 0.065, 0.0415, 0.0472, 0.042, 0.057, 1/155, 1/440],
        # NB: o3 upper bound 0.06 here (exploratory R runs); the published
        # pipeline in func.jl uses 0.0594, see Milankovitch_targets.
        fhigh = [0.07, 0.0785, 0.0785, 0.044, 0.0505, 0.0505, 0.06,
                 1/75, 1/370]
    )

    window_size = 1e3
    cut_groups = tmin:window_size:tmax |> collect
    cut_groups[end] = cut_groups[end] + 1

R"""
my_peak <- function(data) {
    data |>
        dplyr::select(all_of(c("time", "filter"))) |>
        astrochron::peak(plateau = FALSE, genplot = FALSE, verbose = FALSE) |>
        dplyr::rename(peak_age = Location, peak_value = Peak_Value) |>
        dplyr::mutate(diff = peak_age - dplyr::lag(peak_age)) |>
        tidyr::drop_na()
}
"""

    pks = rcopy(R"""
cp_flt = $(x) |>
  # this is from my R package CretaceousConstraints
  bandpass_filter(frequencies = $(filter_freqs), x = time, y = etp) |>
  # taner_filter(frequencies = $(filter_freqs), x = time, y = etp, roll = 1e3, verbose = FALSE, genplot = FALSE) |>
  tidyr::nest(.by = "target") |>
  dplyr::mutate(pk = purrr::map(data, my_peak)) |>
  dplyr::select(-"data") |>
  tidyr::unnest(cols = "pk") |>
  mutate(grp = cut(peak_age, $(cut_groups)))
""")
    @info "filtered frequencies"
    @info "identied peaks"

    means = rcopy(R"""
cp_flt |>
   summarize(.by = c(grp, target),
             mean_age = mean(peak_age),
             mean_diff_1Myr = mean(diff),
             n = n()) |>
   mutate(solution = $(x.sol[1]))
""")
    # TODO: do something else so it become
    # for (i, gr), (j, target) in enumerate(pks.grp), enumerate(pks.target)
    #     means = bootstrap(mean, pks.diff, MaximumEntropySampling(999))
    # end

    @info "calculated means"
    return means
end
