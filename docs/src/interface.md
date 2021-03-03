Using with JuMP
===============

This section provides an example of how to use the relaxations provided by PolyhdedralRelaxations in conjunction with JuMP.


The following examples illustrate the interfacing code for ``y = x^3`` on the partition ``[-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]``. 

```julia 
# Create MILP relaxation 
using PolyhderalRelaxations, JuMP, Cbc
cbc_optimizer = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

f, partition = x -> x^3, collect(-1.0:0.25:1.0)

# create the MILP relaxation of the univariate function.
milp = Model(cbc_optimizer)
@variable(milp, -1.0 <= x <= 1.0)
@variable(milp, y)
formulation_info_milp = construct_univariate_relaxation!(milp, f, x, y, partition, true)

# create the LP relaxation of the univariate function.
lp = Model(cbc_optimizer)
@variable(lp, -1.0 <= x_lp <= 1.0)
@variable(lp, y_lp)
formulation_info_lp = construct_univariate_relaxation!(milp, f, x_lp, y_lp, partition, false)
```

The `formulation_info_milp` and `formulation_info_lp` contains variable and constraint references for all the additional variables and constraints that are used to formulate the polyhedral relaxation. It definition is as follows:

```@docs 
FormulationInfo
```

The reader is referred to [manuscript](https://arxiv.org/abs/2005.13445) for details on the formulation.

For an exhaustive list of functions, the reader is referred to [PolyhedralRelaxations.jl Function Reference](@ref). 