# PolyhedralRelaxations.jl Examples

This folder contains examples demonstrating how to use PolyhedralRelaxations.jl to construct MILP and LP relaxations for nonlinear programs.

## Setup

From the **repo root**, install the example dependencies (this also makes the local package available):

```
julia --project=examples -e 'using Pkg; Pkg.develop(PackageSpec(path=pwd())); Pkg.instantiate()'
```

## Running All Examples

```
julia --project=examples examples/run_examples.jl
```

## Running a Single Example

```
julia --project=examples examples/trig.jl
```

Or interactively from the Julia REPL:

```julia
using PolyhedralRelaxations, JuMP, HiGHS, Test
include("examples/trig.jl")
```

## Examples

### `trig.jl`

Solves the [trig MINLPLib instance](http://www.minlplib.org/trig.html) — a nonlinear, non-convex minimization problem with trigonometric functions. The example shows how to:

- Reformulate a nonlinear program by introducing auxiliary variables for each nonlinear term
- Construct MILP and LP relaxations using `construct_univariate_relaxation!`
- Use error tolerances to control relaxation tightness
- Verify that relaxation objectives provide valid lower bounds on the optimal value
