Interfacing with JuMP
=====================

This section provides an example of how to interface the relaxations provided by PolyhdedralRelaxations with JuMP.

PolyhdedralRelaxations.jl was intentionally implemented without a direct dependency on JuMP.jl for easier maintenance. Nevertheless, interfacing the relaxations to a JuMP model is pretty straightforward. 

The following examples illustrate the interfacing code for ``y = x^3`` on the partition ``[-1.0, -0.75, -0.5, -0.25, 0.0, 0.25, 0.5, 0.75, 1.0]``. The MILP relaxation for the function consists constraints of type ``=`` and ``\leq`` where as the LP relaxation purely consists of constraints of type ``=``.

```julia 
using PolyhderalRelaxations, JuMP, SparseArrays, Cbc
cbc_optimizer = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

f, partition = x -> x^3, collect(-1.0:0.25:1.0)

# create the MILP relaxation of the univariate function.
milp_relaxation, milp_function_data = construct_milp_relaxation(f, partition)
milp = Model(cbc_optimizer)
lb, ub = get_variable_bounds(milp_relaxation)
num_variables = get_num_variables(milp_relaxation)
@variable(milp, lb[i] <= x[i = 1:num_variables] <= ub[i],
    binary = Bool(get_variable_type(milp)[i]),
    base_name = get_variable_names(milp)[i]
)
A, b = get_eq_constraint_matrices(milp_relaxation)
@constraint(milp, A * x .== b)
A, b = get_leq_constraint_matrices(milp_relaxation)
@constraint(milp, A * x .<= b)

# create LP relaxation of the univariate function using the convex hull formulation.
lp_relaxation, lp_function_data = construct_lp_relaxation(f, partition)
lp = Model(cbc_optimizer)
lb, ub = get_variable_bounds(lp_relaxation)
num_variables = get_num_variables(lp_relaxation)
@variable(lp, lb[i] <= y[i = 1:num_variables] <= ub[i])
A, b = get_eq_constraint_matrices(lp_relaxation)
@constraint(lp, A * y .== b)
```

## Helper Functions
The following helper functions can be used to query the properties of both the relaxations. 

```@docs 
get_variable_bounds
get_variable_names
has_geq_constraints
get_geq_constraint_matrices
has_eq_constraints
has_leq_constraints
get_eq_constraint_matrices
get_leq_constraint_matrices
get_num_variables
get_error_bound
get_variable_type
get_domain
get_domain_lb
get_domain_ub
get_partition
```

For an exhaustive list of functions, the reader is referred to [PolyhedralRelaxations.jl Function Reference](@ref). 