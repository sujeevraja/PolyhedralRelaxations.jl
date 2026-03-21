#!/usr/bin/env bash
set -euo pipefail

# Run from the repo root
cd "$(dirname "$0")/.."

julia --project=. -e 'using Pkg; Pkg.test()'
