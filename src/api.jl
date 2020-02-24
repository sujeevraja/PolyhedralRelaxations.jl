export
    construct_lp_relaxation,
    construct_milp_relaxation,
    get_variable_bounds,
    get_variable_names,
    has_geq_constraints,
    get_geq_constraint_matrices,
    get_error_bound, has_eq_constraints,
    has_leq_constraints,
    get_eq_constraint_matrices,
    get_leq_constraint_matrices,
    get_num_variables,
    has_eq_constraints,
    get_variable_type,
    get_function,
    get_derivative,
    get_domain_lb,
    get_domain_ub,
    get_domain,
    get_partition


function construct_milp_relaxation(
    f::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_binary_variables::Int64 = 0,
)::Pair{MILPRelaxation,FunctionData}
    return construct_milp_relaxation(
        f,
        x -> ForwardDiff.derivative(f, x),
        base_partition,
        error_tolerance = error_tolerance,
        length_tolerance = length_tolerance,
        derivative_tolerance = derivative_tolerance,
        num_additional_binary_variables = num_additional_binary_variables,
    )
end

function construct_milp_relaxation(
    f::Function,
    f_dash::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_binary_variables::Int64 = 0,
)::Pair{MILPRelaxation,FunctionData}
    function_data = FunctionData(
        f,
        f_dash,
        base_partition,
        copy(base_partition),
        error_tolerance,
        length_tolerance,
        derivative_tolerance,
        num_additional_binary_variables,
    )
    validate(function_data)
    refine_partition!(function_data)
    return build_milp_relaxation(function_data)
end

function construct_lp_relaxation(
    f::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_binary_variables::Int64 = 0,
)::Pair{LPRelaxation,FunctionData}
    return construct_lp_relaxation(
        f,
        x -> ForwardDiff.derivative(f, x),
        base_partition,
        error_tolerance = error_tolerance,
        length_tolerance = length_tolerance,
        derivative_tolerance = derivative_tolerance,
        num_additional_binary_variables = num_additional_binary_variables,
    )
end

function construct_lp_relaxation(
    f::Function,
    f_dash::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_binary_variables::Int64 = 0,
)::Pair{LPRelaxation,FunctionData}
    function_data = FunctionData(
        f,
        f_dash,
        base_partition,
        copy(base_partition),
        error_tolerance,
        length_tolerance,
        derivative_tolerance,
        num_additional_binary_variables,
    )
    validate(function_data)
    refine_partition!(function_data)
    return build_lp_relaxation(function_data)
end

"Getters for AbstractFormulation"
@inline get_variable_bounds(formulation::AbstractFormulation) =
    formulation.lower_bounds, formulation.upper_bounds
@inline get_variable_names(formulation::AbstractFormulation) = formulation.variable_names
@inline has_geq_constraints(formulation::AbstractFormulation) = false
@inline get_geq_constraint_matrices(formulation::AbstractFormulation) =
    Memento.error(_LOGGER, "both the LP and the MILP relaxation have no >= constraints")
@inline get_error_bound(formulation::AbstractFormulation) = formulation.error_bound

"Getters for LPRelaxation"
@inline has_eq_constraints(lp::LPRelaxation) = true
@inline has_leq_constraints(lp::LPRelaxation) = false
@inline get_eq_constraint_matrices(lp::LPRelaxation) = lp.A, lp.b
@inline get_leq_constraint_matrices(lp::LPRelaxation) =
    Memento.error(_LOGGER, "the LP relaxation has no <= constraints")
@inline get_num_variables(lp::LPRelaxation) = length(lp.λ_indices) + 2

"Getters for MILPRelaxation"
@inline has_eq_constraints(milp::MILPRelaxation) = true
@inline has_leq_constraints(milp::MILPRelaxation) = true
@inline get_eq_constraint_matrices(milp::MILPRelaxation) = milp.A_eq, milp.b_eq
@inline get_leq_constraint_matrices(milp::MILPRelaxation) = milp.A_leq, milp.b_leq
@inline get_variable_type(milp::MILPRelaxation) = milp.binary
@inline get_num_binary_variables(milp::MILPRelaxation) = length(milp.z_indices)
@inline get_num_variables(milp::MILPRelaxation) =
    length(milp.δ_1_indices) + length(milp.δ_2_indices) + length(milp.z_indices) + 2
@inline get_binary_indices(milp::MILPRelaxation) = findnz(milp.binary)[1]

"Getters for FunctionData"
@inline get_function(function_data::FunctionData) = function_data.f
@inline get_derivative(function_data::FunctionData) = function_data.f_dash
@inline get_domain_lb(function_data::FunctionData) = function_data.partition[1]
@inline get_domain_ub(function_data::FunctionData) = function_data.partition[end]
@inline get_domain(function_data::FunctionData) =
    get_domain_lb(function_data), get_domain_ub(function_data)
@inline get_partition(function_data::FunctionData) = function_data.partition
