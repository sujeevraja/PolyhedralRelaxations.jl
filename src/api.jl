export construct_univariate_relaxation!

"""
    construct_univariate_relaxation!(m,x,y,f,x_partition;f_dash=x->ForwardDiff.derivative(f,x),error_tolerance=NaN64,length_tolerance=ϵ,derivative_tolerance=ϵ,num_additional_partitions=0)

Add MILP relaxation of `y=f(x)` to given JuMP model and return an object with
new variables and constraints.

# Mandatory Arguments
- `m::Jump.Model`: model to which relaxation is to be added.
- `f::Function`: function or oracle for which a polyhedral relaxation is
    required, usually non-linear.
- `x::Jump.VariableRef`: JuMP variable for domain of `f`.
- `y::JuMP.VariableRef`: JuMP variable for evaluation of `f`.
- `x_partition::Vector{<:Real}`: partition of the domain of `f`.
- `milp::Bool`: build MILP relaxation if true, LP relaxation otherwise. Note
    that the MILP relaxation uses the incremental formulation presented in the
    paper, but the LP relaxation uses a lambda form that builds a formulation
    as the convex combination of triangle vertices that are at the intersection
    of tangents to the curve.

# Optional Arguments
- `f_dash::Function`: function or oracle for derivative of `f`, defaults to 
    the `derivative` function from the `ForwardDiff` package.
- ` error_tolerance::Float64`: Maximum allowed vertical distance between over
    and under estimators of `f`, defaults to NaN64.
- `length_tolerance::Float64`: maximum length of a sub-interval in a partition,
    defaults to 0.001.
- `derivative_tolerance::Float64`: minimum absolute difference between
    derivaties at successive elements of a partition for them to be considered
    different, defaults to 0.001. If the difference of a partition sub-interval
    is smaller than this value, that sub-interval will be refined.
- `num_additional_partitions::Int64`: budget on number of sub-intervals in
    partition, defaults to 0. Note that if the number of partitions is `n` and
    the number of additional partitions is `m`, then the function will return a
    relaxation with at most `n+m` partitions.

Assume that:
- `f` is a bounded function of 1 variable.
- `x_partition` is a partition of the domain of `f` such that `f` is either
    convex or concave each sub-interval of the partition.
- `f_dash` is not equal at two consecutive elements of `x_partition`.

This function builds an incremental formulation, which is the formulation
presented in the paper.
"""
function construct_univariate_relaxation!(
    m::JuMP.Model,
    f::Function,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    x_partition::Vector{<:Real},
    milp::Bool;
    f_dash::Function = x -> ForwardDiff.derivative(f, x),
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = EPS,
    derivative_tolerance::Float64 = EPS,
    num_additional_partitions::Int64 = 0,
)::FormulationInfo
    univariate_function_data = UnivariateFunctionData(
        f,
        f_dash,
        x_partition,
        error_tolerance,
        length_tolerance,
        derivative_tolerance,
        num_additional_partitions,
        length(x_partition),
    )
    _validate(univariate_function_data)
    _validate(x, x_partition)
    _refine_partition!(univariate_function_data)
    func = milp ? _build_univariate_milp_relaxation! : _build_univariate_lp_relaxation!
    return func(m, x, y, univariate_function_data)
end

"""
    construct_bilinear_relaxation!(m,x,y,z,partition_x,partition_y)

Add polyhedral relaxation of `z = xy` to given JuMP model and return an object with
new variables and constraints.

# Mandatory Arguments
- `m::Jump.Model`: model to which relaxation is to be added.
- `x::Jump.VariableRef`: JuMP variable `x`.
- `y::JuMP.VariableRef`: JuMP variable `y`.
- `z::JuMP.VariableRef`: JuMP variable `z`.
- `partition_x::Vector{<:Real}`: partition of the domain of `x`.
- `partition_y::Vector{<:Real}`: partition of the domain of `y`.

This function builds an incremental formulation, and currently supports more than 
one partition only on one of the variables `x` or `y` and not on both. It will 
throw an error when more than one partitions are provided on both variables. 
When exactly one partition is input for both variables, it populates the model 
with the McCormick relaxation. The incremental formulation is similar to the triangle 
chain relaxation in the manuscript with the triangles replaced with quadrilaterals. 
"""
function construct_bilinear_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    z::JuMP.VariableRef,
    partition_x::Vector{<:Real},
    partition_y::Vector{<:Real}
)::FormulationInfo
    _validate(x, y, partition_x, partition_y)
    return build_bilinear_relaxation!(m, x, y, z, partition_x, partition_y)
end 
