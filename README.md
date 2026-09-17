# Hakatai cyclostratigraphy — ZB23.N64 cycle-duration analysis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21038841.svg)](https://doi.org/10.5281/zenodo.21038841)

Julia code accompanying the Middle Hakatai Shale cyclostratigraphy manuscript.

The code downloads the [ZB23.N64 astronomical solutions](https://www.soest.hawaii.edu/oceanography/faculty/zeebe_files/Astro.html)
(Zeebe & Lantink, 2024, *The Astronomical Journal*: 64 solutions over the past
3.5 Gyr), builds a weighted ETP signal, bandpass-filters it at the Milankovitch
target frequencies, detects peaks, and reports mean cycle durations per 1 Myr
block with 95% intervals (empirical 2.5th and 97.5th percentiles across the 63
solutions and blocks; no bootstrap, the pipeline is deterministic).

## Quick start

Requires [Julia](https://julialang.org/downloads/) 1.12 (the pinned
`Manifest.toml` was generated on 1.12.5 and will not instantiate on older
releases). No R required.

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'   # once, installs pinned deps
cp -r out out.reference                                # keep the shipped tables to compare against
julia --project=. run_analysis.jl
```

The shipped `out/*.csv` and `imgs/*.png` are the reference output and are
overwritten by the run. The pipeline is deterministic, so the four regenerated
tables should be byte-identical to your copies (`diff out.reference out`).

**Cost of a cold run.** Each of the 63 solutions is a ~150 MB zip of the full
3.5 Gyr record, and it is downloaded once per time window: two windows, so 126
downloads and ~19 GB of transfer; allow hours. Only the requested window is
kept, cached as `out/ZB23.R*_<tmin>-<tmax>.arrow` (~440 MB in total), and a
checksum of each window's data is compared with `CHECKSUMS.txt` (a warning is
printed if the server's file has changed). Later runs use the cache and take
minutes. The caches are gitignored: they are reproducible downloads, not
results.

If your `startup.jl` calls `Pkg.activate(".")`, it overrides `--project` in
the `Pkg.instantiate()` step above; add `--startup-file=no` for that step. The
scripts themselves activate their environment explicitly.

Outputs:

| file | contents | shipped |
| --- | --- | --- |
| `out/ZB23.N64_filtered_durations_-1205_-1200_.csv` | 1 Myr block mean durations, −1205 to −1200 Myr | yes |
| `out/ZB23.N64_filtered_durations_-1250_-1200_.csv` | 1 Myr block mean durations, −1250 to −1200 Myr | yes |
| `out/summary_durations_*.csv` | mean duration per target with 95% interval | yes |
| `out/ratios_julia_-1205_-1200.csv` | durations pivoted wide, for target ratios | no |
| `imgs/duration_vs_time.png` | duration through time, per solution | yes |
| `imgs/cycle_duration_means_both.png` | duration densities for both windows | yes |

The trailing underscore in the two duration file names is deliberate: it
distinguishes this output (with `avg_time` and `avg_amp`) from an earlier
grouping used during exploration.

Add `--all-peaks` to also write every individual peak-to-peak duration and
their summaries (`out/ZB23.N64_filtered_all_duration_*.csv` and
`out/summary_all_durations_*.csv`; slow, ~70 MB):

```sh
julia --project=. run_analysis.jl --all-peaks
```

## What is in here

| file | role |
| --- | --- |
| `run_analysis.jl` | **the reproducible pipeline.** Pure Julia, runs top to bottom. |
| `func.jl` | core routines: `get_ZB23`, `bandpass`, `ecc_wrapper_julia` (plus bootstrap helpers the pipeline does not use). No R. |
| `func_R.jl` | *optional* astrochron wrappers. Needs RCall + R. |
| `crosscheck_astrochron.jl` | *optional* check that the Julia filter matches astrochron. |
| `astronomical_boot.jl` | the original exploratory working script, kept as a record. Uses RCall, so the `crosscheck` environment. Step through it; it does not run top to bottom. |
| `filter.jl` | scratch comparison of filter padding strategies, for the REPL. |
| `Project.toml` / `Manifest.toml` | pinned environment for the pipeline. `[compat]` holds exact pins matching the Manifest. |
| `crosscheck/Project.toml` / `Manifest.toml` | the same pins, plus RCall, for the optional R cross-check. |
| `CHECKSUMS.txt` | sha256 of the data in each downloaded ZB23 window. |
| `make_release.sh` | assemble, verify and zip a clean copy of just this code for sharing. |

The published cycle durations come from `ecc_wrapper_julia` in `func.jl` —
the **pure-Julia** path. R is not involved in any reported number.

The four duration/summary tables in `out/` and the two figures in `imgs/` are
tracked and shipped as the reference output. Everything else under `out/` (the
ZB23 caches, the all-peaks tables, the ratios) is regenerated and gitignored.

## The optional astrochron cross-check

The `bandpass` function in `func.jl` was written to reproduce
`astrochron::bandpass`. `crosscheck_astrochron.jl` is intended to demonstrate
that the two agree, and is the reason the R code is kept here rather than
deleted. It takes its passbands from `Milankovitch_targets` in `func.jl`, so it
tests the bands the published numbers use, and writes
`out/crosscheck_astrochron_vs_julia.csv`, `imgs/crosscheck_astrochron_vs_julia.png`
and `out/crosscheck_versions.txt` (the R, astrochron and CretaceousConstraints
versions used). These are shipped when present.

It needs a separate environment, because it adds RCall (which will not install
without a working R):

```sh
julia --project=crosscheck -e 'using Pkg; Pkg.instantiate()'
julia --project=crosscheck crosscheck_astrochron.jl
```

plus, in R:

```r
install.packages(c("tidyverse", "astrochron", "remotes"))
remotes::install_github("japhir/CretaceousConstraints")
```

[`CretaceousConstraints`](https://github.com/japhir/CretaceousConstraints)
provides `bandpass_filter()` and `taner_filter()`.

## Notes for re-running

- **ZB23.R61 is skipped throughout** — it is not published on the server.
- The ETP weights are `1.5, 1.2, 1.2` (eccentricity, obliquity, precession),
  from the Figure S12 caption.
- Filter passbands are in `Milankovitch_targets` in `func.jl`, in 1/kyr. The
  exploratory R-side runs in `func_R.jl` used an upper bound of 0.06 for o3
  where the pipeline uses 0.0594; the R path is not used for any published
  number.
- The solution sampling interval is 0.4 kyr; this is passed to `bandpass`
  explicitly and is not derived from the data.
- **The window read is two samples (0.8 kyr) wider than stated.** `get_ZB23`
  positions the read by line number and the solution files carry a comment
  line, so each window starts at −1 199 999.2 kyr rather than −1 200 000.0
  (12 503 rows instead of 12 501 for the 5 Myr window). The published tables
  were computed this way. Reading the exact window instead moves about 5% of
  the block means by one sample (0.4 kyr), which bounds the shift of any
  summary mean to 0.02 kyr, so it was left as is and documented; to change it,
  filter the loaded rows on `tmin <= time <= tmax` in `get_ZB23`.
- **Long-period targets depend on the window edges.** The ETP is z-scored over
  the window and the filter zero-pads it, so the filtered signal near the ends
  depends on the window. For E (405 kyr, only ~12 cycles in the 5 Myr window)
  the same solution and block can differ between the two windows by up to
  ~70 kyr, comparable to the reported 95% interval; for e up to ~13 kyr; for
  the obliquity and precession targets by less than 0.3 kyr. The 50 Myr window
  is the more robust one for E and e.
- Nothing in `run_analysis.jl` draws random numbers, so no seed is needed. The
  bootstrap helpers in `func.jl` are only exercised from the exploratory script.

## Relationship to the inversion spreadsheet

The archived deposit contains both these CSVs and the cyclostratigraphic
inversion model spreadsheet. The same numbers appear in both, deliberately:

| | |
| --- | --- |
| `out/ZB23.N64_filtered_durations_-1205_-1200_.csv` | = tab `raw ZB23 outputs 1200-1205 Ma` |
| `out/ZB23.N64_filtered_durations_-1250_-1200_.csv` | = tab `raw ZB23 outputs 1200-1250 Ma` |

**The CSVs are canonical.** They are the direct, unedited output of
`run_analysis.jl`. The spreadsheet tabs are a copy of them, pasted in, with the
405 kyr target renamed from `E` to `longecc`. Row order is identical, which
was verified row-for-row on (`grp`, `avg_duration`, `solution`) with zero
mismatches across all 2,835 and 28,350 rows.

Only `avg_duration` is used downstream: all five pivot tables in the workbook
declare it as their single data field, and no formula there references any
other column of these tables.

## Known issues

**Intermittent segfault at exit.** `run_analysis.jl` occasionally prints
`signal 11 (1): Segmentation fault` *after* logging `done`. It is a shutdown
race in the plotting stack, happens after every file is written and flushed,
and the exit code is still 0. Outputs are unaffected.

## Provenance

`astronomical_boot.jl` is the working script the analysis grew out of. The
edits that made it portable (environment, GLMakie to CairoMakie, the split of
`func.jl` into a pure-Julia and an R part) are marked `reproducibility fix` in
the source; the git history has the details.

## Citation and archive

The archived release of this code, together with the rest of the data for the
paper (including the inversion spreadsheet), is at
<https://doi.org/10.5281/zenodo.21038841>. That DOI resolves to the latest
version of the deposit. Please cite it alongside the paper when reusing this
code; `CITATION.cff` has the metadata.

## License

Code: GPL-3.0-or-later, see `LICENSE.md`. The ZB23 astronomical solutions are
redistributed by their authors at the URL above and are not included here.
