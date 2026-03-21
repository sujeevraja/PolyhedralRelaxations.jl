#!/usr/bin/env bash
set -euo pipefail

# Run from the repo root
cd "$(dirname "$0")/.."

julia --project=. -e 'using Pkg; Pkg.test(coverage=true)'

julia --project=. -e '
    using Coverage
    coverage = process_folder("src")
    covered, total = get_summary(coverage)
    pct = round(100 * covered / total; digits=1)
    println("\n=== Coverage Summary ===")
    println("Lines covered: $covered / $total ($pct%)")
    clean_folder("src")
'
