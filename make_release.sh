#!/usr/bin/env bash
# Assemble a clean, self-contained copy of the analysis for sharing / Zenodo.
#
#     ./make_release.sh [destination]     # default: ../hakatai-cyclostratigraphy
#
# Copies the code, the pinned environments, the reference tables and figures;
# verifies every copied file against the working tree; writes
# RELEASE_MANIFEST.txt (Julia version + sha256 of every shipped file) into the
# release; and zips the release to <destination>_code.zip next to this script.
# Fails if any reference output is missing: run run_analysis.jl first.
#
# Leaves behind the 3D scans, .blend files, manuscript drafts, the downloaded
# ZB23.RXX caches (run_analysis.jl re-downloads them) and the all-peaks tables.

set -euo pipefail
src="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dest="${1:-$(dirname "$src")/hakatai-cyclostratigraphy}"

code=(func.jl func_R.jl run_analysis.jl crosscheck_astrochron.jl
      astronomical_boot.jl filter.jl
      Project.toml Manifest.toml crosscheck/Project.toml crosscheck/Manifest.toml
      README.md LICENSE.md CITATION.cff CHECKSUMS.txt .gitignore make_release.sh)
# The duration tables the inversion workbook consumes, their summaries, and
# the two figures run_analysis.jl regenerates.
outputs=(out/ZB23.N64_filtered_durations_-1205_-1200_.csv
         out/ZB23.N64_filtered_durations_-1250_-1200_.csv
         out/summary_durations_-1205_-1200.csv
         out/summary_durations_-1250_-1200.csv
         imgs/duration_vs_time.png
         imgs/cycle_duration_means_both.png)
# optional: only present once the R cross-check has been run
optional=(out/crosscheck_astrochron_vs_julia.csv out/crosscheck_versions.txt
          imgs/crosscheck_astrochron_vs_julia.png)

for f in "${code[@]}" "${outputs[@]}"; do
    [ -f "$src/$f" ] || { echo "missing $f -- run run_analysis.jl first" >&2; exit 1; }
done

# start from an empty release (keep its git history, if any)
mkdir -p "$dest"
find "$dest" -mindepth 1 -not -path "$dest/.git" -not -path "$dest/.git/*" -delete

shipped=()
for f in "${code[@]}" "${outputs[@]}"; do
    mkdir -p "$dest/$(dirname "$f")"
    cp "$src/$f" "$dest/$f"
    shipped+=("$f")
done
for f in "${optional[@]}"; do
    if [ -f "$src/$f" ]; then
        mkdir -p "$dest/$(dirname "$f")"; cp "$src/$f" "$dest/$f"; shipped+=("$f")
    else
        echo "note: optional $f not present (R cross-check not run); skipped"
    fi
done

# verify the copy is exactly the working tree
for f in "${shipped[@]}"; do cmp "$src/$f" "$dest/$f"; done

{
    echo "Hakatai cyclostratigraphy release, $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "built with: $(julia --startup-file=no --version)"
    echo
    (cd "$dest" && sha256sum "${shipped[@]}")
} > "$dest/RELEASE_MANIFEST.txt"

zip="$src/$(basename "$dest")_code.zip"
python3 - "$dest" "$zip" <<'PY'
import os, sys, zipfile
dest, out = sys.argv[1], sys.argv[2]
name = os.path.basename(dest)
with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
    for root, dirs, files in os.walk(dest):
        dirs[:] = [d for d in dirs if d != ".git"]
        for f in sorted(files):
            p = os.path.join(root, f)
            z.write(p, os.path.join(name, os.path.relpath(p, dest)))
PY

echo "Wrote release to: $dest  ($(du -sh "$dest" --exclude=.git | cut -f1)), ${#shipped[@]} files verified"
echo "Wrote zip:        $zip"
echo
echo "Next:"
echo "  cd $dest && git add -A && git commit -m 'Hakatai cyclostratigraphy analysis v1.0.0' && git tag v1.0.0"
