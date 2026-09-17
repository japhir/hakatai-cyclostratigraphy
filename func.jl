# Pure-Julia core routines for the Hakatai cyclostratigraphy analysis.
#
# This file has NO dependency on R / RCall: it can be `include`d with only the
# packages in the top-level Project.toml. The astrochron-based cross-check
# wrappers live in `func_R.jl`, which is optional (see README.md).
#
# Expects the following to be in scope (see run_analysis.jl):
#   Downloads, ZipFile, CSV, Arrow, DataFrames, DSP, Peaks, StatsBase, Bootstrap, SHA

"""
    data_checksum(dat)

sha256 over the raw bytes of the five data columns of a ZB23 window
(time, ecc, inc, epl, cp), so a re-download can be compared with the data the
published results used independently of the container format.
"""
function data_checksum(dat::DataFrame)
    ctx = SHA256_CTX()
    for c in (:time, :ecc, :inc, :epl, :cp)
        update!(ctx, reinterpret(UInt8, Vector{Float64}(dat[!, c])))
    end
    bytes2hex(digest!(ctx))
end

"""
    verify_checksum(dat, key; list = "CHECKSUMS.txt")

Compare `data_checksum(dat)` with the entry for `key` (the cache file name) in
CHECKSUMS.txt next to func.jl. Warns, never throws: the ZB23 server may
legitimately update its files, but then the numbers here may not reproduce.
"""
function verify_checksum(dat::DataFrame, key::AbstractString;
                         list = joinpath(@__DIR__, "CHECKSUMS.txt"))
    isfile(list) || return nothing
    expected = Dict(String(strip(p[2])) => String(p[1])
                    for p in (split(l, r"\s+", limit = 2) for l in eachline(list) if !startswith(l, "#"))
                    if length(p) == 2)
    if !haskey(expected, key)
        @warn "no checksum listed for $key in CHECKSUMS.txt"
    elseif (actual = data_checksum(dat)) != expected[key]
        @warn "checksum mismatch for $key: the downloaded data differ from the data the published results used" expected = expected[key] actual
    else
        @info "checksum ok"
    end
    return nothing
end

"""
    get_ZB23_full(N)

Download (or load from `out/` cache) the complete 3.5 Gyr ZB23.R`N` solution.
Large: ~8.75 million rows per solution. `get_ZB23(N, tmin, tmax)` below reads
only a time window and is what the analysis uses.
"""
function get_ZB23_full(N::Integer = 1)
    if N > 64 || N < 1
        throw(ArgumentError("N should be between 1 and 64, got $N"))
    end

    base_url = "https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/"

    number = lpad(N, 2,'0')
    astronomical_solution = "ZB23.R$(number)"
    @info astronomical_solution

    mkpath("out")
    outfile = "out/$(astronomical_solution).arrow"
    if isfile(outfile)
        @info "loading file from cache"
        return Arrow.Table(outfile) |> DataFrame
    end

    url = "$(base_url)3.5Gyr/ZB23-N64-eiop/$(astronomical_solution).eiop.dat.zip"

    dwn = Downloads.download(url)
    @info "downloaded"
    unz = ZipFile.Reader(dwn)
    @info "unzipped"

    cols = [
        "time", # time in kyr, negative
        "ecc", # eccentricity
        "inc", # inclination
        "epl", # obliquity
        "cp" # climatic precession e*sin(omegabar)
    ]
    dat = CSV.File(unz.files[1];
                   header = cols,
                   comment = "%",
                   ignorerepeated = true,
                   delim = " ",
                   stripwhitespace = true) |>
                       DataFrame
    @info "read"

    dat.sol .= astronomical_solution
    @info "added name"

    Arrow.write(outfile, dat)
    @info "saved Arrow"

    rm(dwn)
    @info "removed dowloaded file"

    return dat
end

function get_ZB23(N::Integer = 1,
                  tmin = -3.5e6,
                  tmax = 0)
    # these are properties of the solutions
    tend = -3.5e6
    tres = -0.4

    if N > 64 || N < 1
        throw(ArgumentError("N should be between 1 and 64, got $N"))
    end
    skp = convert(Int, floor(tmax / tres))
    ftskp = convert(Int, floor(tend / tres - tmin / tres))

    base_url = "https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro/"

    number = lpad(N, 2,'0')
    astronomical_solution = "ZB23.R$(number)"
    @info astronomical_solution

    mkpath("out")
    outfile = "out/$(astronomical_solution)_$(tmin)-$(tmax).arrow"
    if isfile(outfile)
        @info "loading file from cache"
        return Arrow.Table(outfile) |> DataFrame
    end

    # # see if we can reuse a bigger one?
    # fls = readdir("out")
    # my_parse = function(x)
    #     if isnothing(x.captures)
    #         return nothing
    #     else
    #         parse(Float64, x.captures)
    #     end
    # end

    # mins = my_parse.(match.(r"_(-[0-9.e]+)", fls))
    # maxs = parse.(Float64, map(x -> x.captures, match.(r"--([0-9e.]+).arrow", fls)))
    # rx = astronomical_solution * r"_-[0-9.e]+--[-0-9.e]+.arrow$"
    # mtch = match.(rx, fls)
    # for (i, fl) in enumerate(fls)
    #     if isnothing(mtch[i])
    #         break
    #     else
    #         rx2 = "_-[0-9.e]+"
    #         fl

    url = "$(base_url)3.5Gyr/ZB23-N64-eiop/$(astronomical_solution).eiop.dat.zip"

    dwn = Downloads.download(url)
    @info "downloaded"
    unz = ZipFile.Reader(dwn)
    @info "unzipped"

    cols = [
        "time", # time in kyr, negative
        "ecc", # eccentricity
        "inc", # inclination
        "epl", # obliquity
        "cp" # climatic precession e*sin(omegabar)
    ]
    dat = CSV.File(unz.files[1];
                   header = cols,
                   comment = "%",
                   ignorerepeated = true,
                   delim = " ",
                   skipto = skp, footerskip = ftskp,
                   stripwhitespace = true) |>
                       DataFrame
    @info "read"
    if tmin != tend && tmax != 0
        @info "subset"
    end

    rm(dwn)
    @info "removed file"

    dat.sol .= astronomical_solution
    @info "added name"

    verify_checksum(dat, basename(outfile))

    Arrow.write(outfile, dat)
    @info "saved Arrow"

    return dat
end

# Filter passbands (flow, fhigh) in 1/kyr for each Milankovitch target. This is
# the single definition used by ecc_wrapper_julia below.
const Milankovitch_targets = Dict(
    :p1 => (0.0653, 0.07),
    :p2 => (0.075, 0.0785),
    :p => (0.065, 0.0785),
    :o1 => (0.0415, 0.044),
    :o2 => (0.0472, 0.0505),
    :o12 => (0.042, 0.0505),
    :o3 => (0.057, 0.0594),
    :e => (1/155, 1/75),
    :E => (1/440, 1/370)
)


"""
    bandpass(data, dt, flow, fhigh; filter_order = 4, pad_fac = 2)

Zero-phase Butterworth bandpass of a regularly sampled series, written to
mimic `astrochron::bandpass`.

- `data`: the input series (Vector{Float64}).
- `dt`: sampling interval, in kyr.
- `flow`, `fhigh`: passband edges, in 1/kyr.
- `filter_order`: Butterworth order (default 4).
- `pad_fac`: the series is zero-padded on both sides by `pad_fac * length(data)`
  samples before `filtfilt`, and the padding is removed afterwards (default 2,
  i.e. the padded series is five times the input length).

Because of the padding and because the ETP is z-scored over the analysis
window, the filtered values near the window ends depend on the window; see
README "Notes for re-running".

(Initial version generated by DeepSeek, then revised.)
"""
function bandpass(data::Vector{Float64}, dt::Float64, flow::Float64, fhigh::Float64; filter_order::Int=4, pad_fac::Int=2)
    # Calculate sampling frequency
    fs = 1.0 / dt  # Sampling frequency in 1/kyr

    # Normalize frequencies to the Nyquist frequency (fs/2)
    nyquist = fs / 2
    f_norm_low = flow / nyquist
    f_norm_high = fhigh / nyquist

    # Design a Butterworth bandpass filter
    design = Butterworth(filter_order)
    response = Bandpass(f_norm_low, f_norm_high)

    # Create the digital filter
    filt = digitalfilter(response, design)

    # Pad the signal
    pad_length = length(data) * pad_fac
    # I tested that:
    # with zeros the output is visually very similar
    # to astrochron::bandpass(window = 0) (rectangular)
    data_padded = vcat(#reverse(data[1:pad_length]),
                       zeros(pad_length),
                       data,
                       # reverse(data[end-pad_length+1:end])
                       zeros(pad_length)
                       )

    # Apply the filter using zero-phase filtering
    filtered_padded = filtfilt(filt, data_padded)

    # Remove the padding
    filtered_data = filtered_padded[pad_length+1:end-pad_length]

    return filtered_data
end

# bootstrapping in 5 non-overlapping chunks
"simple block bootstrap"
my_boot = function(x, fun = mean; b = 69, N = 999)
    n = length(x)
    nblocks = convert(Int, floor(n/b))
    @show nblocks

    megaboot = zeros(N * nblocks)
    # srt = rand(1:nblocks,nblocks)

    cache = zeros(N)
    for i in 1:nblocks # all blocks
        idxmin = (i-1)*b+1
        idxmax = (i-1)*b+b
        cache .= bootstrap(fun,
                           # subset data to rolling blocks of size b
                           x[idxmin:idxmax],
                           BasicSampling(N)).t1 |>
                               Iterators.flatten |>
                               collect
        megaboot[(i-1)*N+1:i*N] .= cache
    end
    return megaboot
end


# bootstrapping for time series
moving_block_bootstrap = function(x, fun = mean; b = 69, N = 999)
    # I did N = 10_000 first but the paper mentions 999, say

    n = length(x)
    rolling_blocks = n - b + 1 # number of possible block start positions

    nblocks = convert(Int, floor(n/b))
    @show nblocks # should be about 5
    # pick nblocks at random
    ss = rand(1:rolling_blocks, nblocks)
    # draw n/b of the blocks, with replacement; block starting at s covers s:s+b-1
    idx = [ s .+ (0:b-1) for s in ss ] |> Iterators.flatten |> collect
    # Then aligning these n/b blocks in the order they were picked, will give the bootstrap observations.
    boot_obs = x[idx]
    boot = bootstrap(fun, boot_obs, BasicSampling(N))
    # megaboot = zeros(N * nblocks)
    # cache = zeros(N)
    # for i in 1:rolling_blocks
    #     cache .= bootstrap(fun,
    #                        # subset data to rolling blocks of size b
    #                        boots[i:(i+b-1)],
    #                        MaximumEntropySampling(N)).t1 |>
    #                            Iterators.flatten |>
    #                            collect
    #     megaboot[(i-1)*N+1:i*N] .= cache
    # end
    return boot.t1 |> Iterators.flatten |> collect
end


"""
    ecc_wrapper_julia(solution, fun = mean, tmin, tmax; all_peaks = false)

The published pipeline for one ZB23 solution: load the window, build the
weighted ETP, bandpass at every `Milankovitch_targets` band, detect peaks,
and return per-1 Myr-block `fun` of the peak-to-peak durations (or, with
`all_peaks = true`, every individual duration).
"""
ecc_wrapper_julia = function(solution = 1,
                             fun = mean,
                             tmin = -1250e3, tmax = -1200e3;
                             all_peaks = false
                             )
    x = get_ZB23(solution, tmin, tmax)

    # create ETP with weights
    # 1.5, 1.2, 1.2 from Figure S12 caption
    etp_weights = [1.5, 1.2, 1.2]
    # uses StatsBase.zscore
    necc = zscore(x.ecc)
    nobl = zscore(x.epl)
    nprec = zscore(x.cp)
    @. x.etp .= etp_weights[1] * necc +
        etp_weights[2] * nobl +
        etp_weights[3] * nprec
    @info "calculated ETP with weights $etp_weights"

    # passbands: the top-level Milankovitch_targets Dict

    # # using_a_dict = function()
    # filtered_targets = Dict{Symbol, Vector{Float64}}()
    # for (k, (flow, fhigh)) in Milankovitch_targets
    #     filtered_targets[k] = bandpass(x.etp, 0.4, flow, fhigh)
    # end
    #     # end
    # using_a_df = function()
    filtered_targets = DataFrame(etp = x.etp)
    for (name, (flow, fhigh)) in Milankovitch_targets
        filtered_targets[!, name] = bandpass(x.etp, 0.4, flow, fhigh)
    end
    df = stack(filtered_targets[!, Not(:etp)], view = true)
        # end
    # @btime using_a_dict()
    # # after first run, 995 ms
    # @btime using_a_df()
    # # 991 ms (sliiightly faster)
    @info "filtered"

    # work on grouped df
    peaks = combine(groupby(df, :variable; sort = false),
                    :value => (z -> (pk = findmaxima(z);
                                     (amp = pk.heights, time = x.time[pk.indices]))) => AsTable)
    # append to df => just work with dicts for now
    # for (name, component) in filtered_targets
    #     x[!, Symbol(name)] = component
    # end
    # # TODO: reshaping is a bit slow, so could also not do it?
    # filt_long = stack(x[!, Not([:ecc, :inc, :epl, :cp, :sol])],
    #                   Not([:time, :etp]),
    #                   variable_name = :target, value_name = :amplitude)
    # peaks = Dict{Symbol, Vector{Float64}}()
    # peak_height = Dict{Symbol, Vector{Float64}}()
    # for (c) in enumerate(eachcol(filtered_targets[!,2:end]))
    #     @show names(c)
    #     fp = findmaxima(c)
    #     peaks[Milankovitch_targets.keys[i]]= x.time[fp.indices]
    # end
    @info "detected peaks"

    diffs = combine(groupby(peaks, :variable; sort = false),
                    :time => (z -> z[1:end-1] .+ (z[2:end] .- z[1:end-1]) ./ 2) => :time, # mean age of duration interval
                    :amp => (z -> z[1:end-1] .+ (z[2:end] .- z[1:end-1]) ./ 2) => :amp, # mean amplitude of the two bounding peaks
                    :time => (z -> z[1:end-1] .- z[2:end]) => :duration)

    @info "calculated durations"

    # # inspection plot
    # plt_raw = AlgebraOfGraphics.data(x) *
    #     mapping(:time, :etp) *
    #     visual(Lines, label = "data")

    # plt_flt = AlgebraOfGraphics.data(df) *
    #     mapping(direct(repeat(x.time, outer = Milankovitch_targets.count)), :value, row = :variable) *
    #     visual(Lines, color = :cyan, label = "filter")

    # plt_pks = AlgebraOfGraphics.data(peaks) *
    #     mapping(:time, :amp, row = :variable) *
    #     visual(Scatter, color = :red, label = "peaks")

    # plt_dff = AlgebraOfGraphics.data(diffs) *
    #     mapping(:time, :amp, :duration => (x -> x / 2), row = :variable) *
    #     visual(Errorbars, color = :orange, label = "duration", direction = :x, whiskerwidth = 15)

    # f = (plt_raw + plt_flt + plt_pks + plt_dff) |> draw
    # save("imgs/illustrate_algorithm.png", f)

    # plt_dur = AlgebraOfGraphics.data(diffs) *
    #     mapping(:time => (x -> x / 1000) => "Time (Myr)", :duration => "Peak duration (kyr)", row = :variable) * visual(Lines)

    # f = draw(plt_dur, facet = (;linkyaxes = :none))
    # save("imgs/discrete_diffs.png", f)
    # temporary return all peaks function
    # peaks.solution .= x.sol[1]
    # return peaks

    # something is up with this! seems like there are only very few
    # unique durations, maybe because of the Nyquist frequency?
    # durs = combine(groupby(diffs, :variable), :duration => unique => :udur)
    # sort!(combine(groupby(durs, :variable), nrow), :variable)


    # opt-in: return every individual peak-to-peak duration, rather than the
    # 1 Myr block means. Reproduces out/ZB23.N64_filtered_all_duration_*.csv
    if all_peaks
        diffs.solution .= x.sol[1]
        sort!(diffs, :variable)
        return diffs
    end

    block_size = 1000 # kyr = 1 Myr

    diffs.grp = floor.(Int, diffs.time / block_size)

    summ = combine(groupby(diffs, [:variable, :grp]; sort = false),
                   :time => fun => :avg_time,
                   :amp => fun => :avg_amp,
                   :duration => fun => :avg_duration,
                   # number of peak-to-peak durations this block mean is
                   # taken over
                   nrow => :n_durations)

    sort!(summ, :variable)

    # grp = Dict{Symbol, Vector{String}}()
    # ns = Dict{Symbol, Vector{Int}}()
    # means = Dict{Symbol, Vector{Float64}}()
    # nmeans = length(cut_groups)-1

    # # TODO
    #     grp[k] = repeat([""],nmeans)
    #     means[k] = zeros(nmeans)
    #     ns[k] = zeros(Int, nmeans)
    #     for i in 1:nmeans
    #         # which ages fall within the cut group?
    #         ids = time .>= cut_groups[i] .&&
    #             time .< cut_groups[i+1]
    #         grp[k][i] = "($(cut_groups[i]),$(cut_groups[i+1])]"
    #         means[k][i] = mean(v[ids])
    #         ns[k].n[i] = length(v[ids])
    #     end
    # end

    @info "calculated means"

    summ.solution .= x.sol[1]
    # return (x.sol[1], peaks, means, ns)
    return(summ)
end
