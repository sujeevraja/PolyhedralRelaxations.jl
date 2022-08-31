export construct_univariate_relaxation!,
    construct_bilinear_relaxation!,
    construct_multilinear_relaxation!,
    add_multilinear_linking_constraints!,
    refine_partition!

"""
    construct_univariate_relaxation!(m,f,x,y,x_partition;f_dash=x->ForwardDiff.derivative(f,x),error_tolerance=NaN64,length_tolerance=1e-6,derivative_tolerance=1e-6,num_additional_partitions=0,variable_pre_base_name="",formulation_info=FormulationInfo())

Add MILP relaxation of `y=f(x)` to given JuMP model and return an object with
new variables and constraints.

# Mandatory Arguments
- `m::JuMP.Model`: model to which relaxation is to be added.
- `f::Function`: function or oracle for which a polyhedral relaxation is
    required, usually non-linear.
- `x::JuMP.VariableRef`: JuMP variable for domain of `f`.
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
    defaults to ``1 \\times 10^{-6}``.
- `derivative_tolerance::Float64`: minimum absolute difference between
    derivaties at successive elements of a partition for them to be considered
    different, defaults to ``1 \\times 10^{-6}``. If the difference of a partition sub-interval
    is smaller than this value, that sub-interval will be refined.
- `num_additional_partitions::Int64`: budget on number of sub-intervals in
    partition, defaults to 0. Note that if the number of partitions is `n` and
    the number of additional partitions is `m`, then the function will return a
    relaxation with at most `n+m` partitions.
- `variable_pre_base_name::AbstractString`: base_name that needs to be added to the auxiliary
    variables for meaningful LP files
- `formulation_info::FormulationInfo` : `FormulationInfo` for another univariate function on the 
    same variable that can be reused for the current relaxation. 

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
    variable_pre_base_name::AbstractString = "",
    formulation_info::FormulationInfo = FormulationInfo(),
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
    func =
        milp ? _build_univariate_milp_relaxation! :
        _build_univariate_lp_relaxation!
    return func(
        m,
        x,
        y,
        univariate_function_data,
        variable_pre_base_name,
        formulation_info,
    )
end

"""
    construct_bilinear_relaxation!(m,x,y,z,x_partition,y_partition)

Add polyhedral relaxation of `z = xy` to given JuMP model and return an object with
new variables and constraints.

# Mandatory Arguments
- `m::Jump.Model`: model to which relaxation is to be added.
- `x::Jump.VariableRef`: JuMP variable `x`.
- `y::JuMP.VariableRef`: JuMP variable `y`.
- `z::JuMP.VariableRef`: JuMP variable `z`.
- `x_partition::Vector{<:Real}`: partition of the domain of `x`.
- `y_partition::Vector{<:Real}`: partition of the domain of `y`.

# Optional Arguments
- `variable_pre_base_name::AbstractString`: base_name that needs to be added to the auxiliary
    variables for meaningful LP files
- `formulation_info::FormulationInfo` : `FormulationInfo` for another bilinear function where the 
    variables that is partitioned has to be same on both; to enable partition variable reuse
        

This function builds an incremental formulation, and currently supports more than 
one partition only on one of the variables `x` or `y` and not on both. It will 
throw an error when more than one partitions are provided on both variables. 
When exactly one partition is input for both variables, it populates the model 
with the McCormick relaxation. The incremental formulation is similar to the triangle 
chain relaxation in the manuscript with the triangles replaced with tetrahedrons. 
Adjacent tetrahedrons share an edge (instead of a point in the triangles case). 
"""
function construct_bilinear_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    z::JuMP.VariableRef,
    x_partition::Vector{<:Real},
    y_partition::Vector{<:Real};
    variable_pre_base_name::AbstractString = "",
    formulation_info::FormulationInfo = FormulationInfo(),
)::FormulationInfo
    _validate(x, y, x_partition, y_partition)
    if length(x_partition) == 2 && length(y_partition) == 2
        return _build_mccormick_relaxation!(m, x, y, z)
    end
    return _build_bilinear_milp_relaxation!(
        m,
        x,
        y,
        z,
        x_partition,
        y_partition,
        variable_pre_base_name,
        formulation_info,
    )
end

"""
    construct_multilinear_relaxation!(m,x,y,z,x_partition,y_partition)

Add polyhedral relaxation of `z = product(x)` to given JuMP model and return an object with
new variables and constraints.

# Mandatory Arguments
- `m::Jump.Model`: model to which relaxation is to be added.
- `x::Tuple`: JuMP variables in the multilinear term as a tuple.
- `z::JuMP.VariableRef`: JuMP variable `z`.
- `partitions::Dict{JuMP.VariableRef,Vector{T}} where {T<:Real}`: partition 
of the domain of the `x` variables.

# Optional Arguments
- `variable_pre_base_name::AbstractString`: base_name that needs 
to be added to the auxiliary variables for meaningful LP files

This function builds an lambda based SOS2 formulation for the piecewise polyhedral relaxation. 
Reference information:
    Kaarthik Sundar, Harsha Nagarajan, Jeff Linderoth, Site Wang, 
    Russell Bent, Piecewise Polyhedral Formulations for a Multilinear Term, 
    https://arxiv.org/abs/2001.00514 
"""
function construct_multilinear_relaxation!(
    m::JuMP.Model,
    x::Tuple,
    z::JuMP.VariableRef,
    partitions::Dict{JuMP.VariableRef,Vector{T}} where {T<:Real},
    variable_pre_base_name::AbstractString = "",
)::FormulationInfo
    _validate(x, partitions)
    lp = all([length(it) == 2 for it in values(partitions)])
    (lp) && (
        return _build_multilinear_convex_hull_relaxation!(
            m,
            x,
            z,
            partitions,
            variable_pre_base_name,
        )
    )
    return _build_multilinear_sos2_relaxation!(
        m,
        x,
        z,
        partitions,
        variable_pre_base_name,
    )
end

"""
    add_multilinear_linking_constraints!(m, info, partitions; max_degree_limit = nothing, 
        helper = Dict())::FormulationInfo

Add linking constraints for the different multilinear relaxations 
in the model using the lambda variables for each relaxation. 

Reference information:
    Jongeun Kim, Jean-Philippe P. Richard, Mohit Tawarmalani, 
    Piecewise Polyhedral Relaxations of Multilinear Optimization, 
    http://www.optimization-online.org/DB_HTML/2022/07/8974.html

# Mandatory Arguments
- `m::Jump.Model`: model to which relaxation is to be added.
- `info::Dict{T,FormulationInfo} where {T<:Any}`: dictionary keyed 
by tuple of variables involved in multilinear term whose value is 
the `FormulationInfo` struct returned by added the multilinear 
relaxation for that term
- `partitions::Dict{JuMP.VariableRef,Vector{T}} where {T<:Real}`: partition 
of the domain of the variables.

# Optional Keyword Arguments
- `max_degree_limit::Union{Nothing,T} where {T<:Int64}`: this is a 
control factor for limit the number of linking constraints added. 
Default value is nothing which will not limit the number of constraints added. 
- `helper::Dict` - default is the empty dictionary. This dictionary 
contains the common subterms that are shared between the different 
multilinear terms. When the function is invoked the first time, 
the `FormulationInfo` returned contains this dictionary in 
`extra[:common_subterm_data]`. When invoked the subsequent times, 
passing this dictionary can save a lot on the computation time required 
to generate these constraints.  
"""
function add_multilinear_linking_constraints!(
    m::JuMP.Model,
    info::Dict{T,FormulationInfo} where {T<:Any},
    partitions::Dict{JuMP.VariableRef,Vector{T}} where {T<:Real};
    max_degree_limit::Union{Nothing,T} where {T<:Int64} = nothing,
    helper = Dict(),
)::FormulationInfo
    if isempty(helper)
        is_needed =
            _check_if_linking_constraints_are_needed(info, max_degree_limit)
        (~is_needed.needed) && (return FormulationInfo())
        linking_constraint_helper = is_needed.linking_constraint_helper
    end

    formulation_info = _build_linking_constraints!(
        m,
        info,
        partitions,
        linking_constraint_helper,
    )

    formulation_info.extra[:common_subterm_data] = linking_constraint_helper
    return formulation_info
end

"""
    refine_partition!(
        partition::Vector{<:Real}, 
        point::T where {T<:Real};
        refinement_type::Symol = :non_uniform, 
        refinement_ratio::Float64 = REFINEMENT_RATIO, 
        refinement_width_tol::Float64 = REFINEMENT_WIDTH_TOL,
        refinement_added_point_width_tolerance::Float64 = REFINEMENT_ADDED_POINT_WIDTH_TOL,
        refine_largest::Bool = true
    )

This function refines a variable domain (refered to as a `partition`) using 
the `point` that is contained in one of the partitions. If the point is not 
contained in the variable domain, the function throws an error message. 

# Mandatory Arguments
- `partition::Vector{<:Real}`: domain of the variable to be refined.
- `point::T where {T<:Real}`: a point contained in the domain.

# Optional Keyword Arguments
- `refinement_type::Symbol`: this variable chooses the type of refinement
to be done. The choices are `:non_uniform`, `:bisect`, `:at_point`, and `:bisect_all`. 
The default is `:non_uniform`
- `refinement_ratio::Float64 = 0.1`: parameter to perform refinement (do not change unless
you know what you are doing). This parameter is applicable only for the `non_uniform` choic
of refinement scheme.  
- `refinement_width_tol::Float64 = 1E-2`: a width of the partition beyond which it is not refined.
Also, if `point` is within `refinement_width_tolerance` of  
- `refinement_added_point_width_tolerance::Float64 = 1E-3` - if the refinement points 
are within `refinement_added_point_width_tolerance` of the left or the right of the partition 
containing the `point`, the corresponding refinements are not performed. 
- `refine_largest::Bool = true`: bisects the largest sub-interval if the width of 
the partition where the `point` lies is less than `refinement_width_tol` i.e., the 
sub-interval containing the point is too small for refinement
"""

function refine_partition!(
    partition::Vector{<:Real},
    point::T where {T<:Real};
    refinement_type::Symbol = :non_uniform,
    refinement_ratio::Float64 = REFINEMENT_RATIO,
    refinement_width_tol::Float64 = REFINEMENT_WIDTH_TOL,
    refinement_added_point_width_tolerance::Float64 = REFINEMENT_ADDED_POINT_WIDTH_TOL,
    refine_largest::Bool = true,
)::RefinementInfo
    if (point < partition[1] || point > partition[end])
        error("$point is not contained in the variable domain")
        return RefinementInfo()
    end
    return _refine_partition!(
        partition,
        point,
        refinement_type,
        refinement_ratio,
        refinement_width_tol,
        refinement_added_point_width_tolerance,
        refine_largest,
    )
end
