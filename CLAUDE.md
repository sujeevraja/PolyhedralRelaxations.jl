# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

**Run all tests:**
```
julia --project=. -e 'using Pkg; Pkg.test()'
```

**Run a single test file** (from the repo root):
```
julia --project=. -e 'include("test/runtests.jl")'
```
Or to run only a specific test file interactively:
```julia
using PolyhedralRelaxations, JuMP, HiGHS, Ipopt, Test
include("test/api_tests.jl")
```

**Format code** (uses `.JuliaFormatter.toml` config):
```
julia -e 'using JuliaFormatter; format(".")'
```

**Build/serve docs:**
```
julia --project=docs docs/make.jl
```

## Architecture

This is a Julia package that constructs MILP and LP polyhedral relaxations for:
1. **Univariate functions** — relaxations of `y = f(x)` over a bounded domain
2. **Bilinear terms** — relaxations of `z = x*y`
3. **Multilinear terms** — relaxations of `z = x1*x2*...*xn`

### Public API (`src/api.jl`)

Five exported functions:
- `construct_univariate_relaxation!(m, f, x, y, x_partition, milp; ...)` — adds univariate MILP or LP relaxation to a JuMP model
- `construct_bilinear_relaxation!(m, x, y, z, x_partition, y_partition; ...)` — adds bilinear relaxation; falls back to McCormick when both partitions have only 2 points
- `construct_multilinear_relaxation!(m, x, z, partitions; ...)` — adds multilinear relaxation; uses convex hull LP when all partitions have only 2 points, otherwise SOS2 MILP
- `add_multilinear_linking_constraints!(m, info, partitions; ...)` — adds Kim-Richard-Tawarmalani linking constraints tightening multiple multilinear relaxations
- `refine_partition!(partition, point; ...)` — refines a variable domain partition around a given point

All `construct_*` functions return a `FormulationInfo` struct containing the added variables and constraints (keyed by `Symbol`).

### Key Data Structures (`src/common.jl`)

- `FormulationInfo` — holds `variables`, `constraints`, `indices`, and `extra` dictionaries populated by the relaxation builders; returned by all public API functions
- `UnivariateFunctionData` — bundles function `f`, derivative `f_dash`, partition vector, and tolerance parameters for univariate relaxations

### Internal Source Files

| File | Purpose |
|------|---------|
| `src/univariate_lp_relaxation.jl` | LP relaxation via convex combination of triangle vertices (tangent/secant intersections) |
| `src/univariate_milp_relaxation.jl` | MILP incremental formulation (chain of triangles) |
| `src/bilinear_relaxation.jl` | McCormick relaxation + incremental MILP via chain of tetrahedra |
| `src/multilinear_relaxation.jl` | Convex hull LP and SOS2-based MILP for multilinear terms |
| `src/mutlilinear_linking.jl` | Linking constraints across multiple multilinear relaxations (note: filename has typo) |
| `src/partition_refinement.jl` | Adaptive partition refinement strategies (`:non_uniform`, `:bisect`, `:at_point`, `:bisect_all`) |
| `src/constants.jl` | Numerical constants (`EPS`, `INF`, refinement tolerances) |
| `src/logging.jl` | Custom logger with `set_logging_level!(:Debug/:Info/:Warn/:Error)` and `silence!()` |

### Relaxation Selection Logic

- Univariate: `milp=true` → incremental triangle-chain MILP; `milp=false` → lambda-form LP
- Bilinear: both partitions length 2 → McCormick; one partition longer → incremental tetrahedra MILP; both longer → error (use multilinear API instead)
- Multilinear: all partitions length 2 → convex hull LP; any partition longer → SOS2 MILP with binary variables

### Code Style

JuliaFormatter is enforced on PRs. Config (`.JuliaFormatter.toml`): 80-char margin, `always_for_in`, `always_use_return`, `short_to_long_function_def`, `remove_extra_newlines`. Internal helpers are prefixed with `_`.
