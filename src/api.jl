export construct_lp_relaxation,
    construct_milp_relaxation,
    construct_univariate_milp_relaxation,
    construct_univariate_lp_relaxation,
    get_variable_bounds,
    get_variable_names,
    has_geq_constraints,
    get_geq_constraint_matrices,
    get_error_bound,
    has_eq_constraints,
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

"""
    construct_univariate_milp_relaxation(m, x, y, f, x_partition; error_tolerance = NaN64, length_tolerance = ϵ, derivative_tolerance = ϵ, num_additional_partitions = 0)

This function is the entry point for constructing a MILP relaxation for the univariate function ``y = f(x)``. The mandatory inputs to the functions are the the JuMP model, the JuMP variables x and y, the function (oracle), the partition on the ``x`` variable's domain. All other arguments are optional. The keyword argument `error_tolerance` is used to obtain an MILP relaxation by adding more partitions to the `base_partition`. Partitions are added till the vertical distance between the over-estimator and under-estimator of the relaxation is less than the `error_tolerance` for any value in the domain. `length_tolerance` is similar in the sense that it prohibits adding more partitions to an interval that is less than this value and `derivative_tolerance` is used to check if derivaties of the function as successive discretization points are not equal to each other. Finally, the `num_additional_partitions` can be used generate an LP relaxation of the function with a budget on the number of partitions. Note that if the number of partitions is ``n`` and the number of additional partitions is ``m``, then the function will return a relaxation with at most ``n+m`` partitions. 
"""
function construct_univariate_milp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    f::Function,
    x_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_partitions::Int64 = 0,
)::FormulationInfo
    return construct_univariate_milp_relaxation!(
        m,
        x,
        y,
        f,
        x -> ForwardDiff.derivative(f, x),
        x_partition,
        error_tolerance = error_tolerance,
        length_tolerance = length_tolerance,
        derivative_tolerance = derivative_tolerance,
        num_additional_partitions = num_additional_partitions,
    )
end

"""
    construct_univariate_milp_relaxation(m, x, y, f,f_dash, x_partition; error_tolerance = NaN64, length_tolerance = ϵ, derivative_tolerance = ϵ, num_additional_partitions = 0)

This function is the entry point for constructing a MILP relaxation for the univariate function ``y = f(x)``. The mandatory inputs to the functions are the the JuMP model, the JuMP variables x and y, the function (oracle), the derivative of f (oracle), the partition on the ``x`` variable's domain. All other arguments are optional. The keyword argument `error_tolerance` is used to obtain an MILP relaxation by adding more partitions to the `base_partition`. Partitions are added till the vertical distance between the over-estimator and under-estimator of the relaxation is less than the `error_tolerance` for any value in the domain. `length_tolerance` is similar in the sense that it prohibits adding more partitions to an interval that is less than this value and `derivative_tolerance` is used to check if derivaties of the function as successive discretization points are not equal to each other. Finally, the `num_additional_partitions` can be used generate an LP relaxation of the function with a budget on the number of partitions. Note that if the number of partitions is ``n`` and the number of additional partitions is ``m``, then the function will return a relaxation with at most ``n+m`` partitions. 
"""
function construct_univariate_milp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    f::Function,
    f_dash::Function,
    x_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
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
    )
    validate(univariate_function_data)
    refine_partition!(univariate_function_data)
    return build_milp_relaxation(univariate_function_data)
end

"""
    construct_univariate_lp_relaxation(m, x, y, f, x_partition; error_tolerance = NaN64, length_tolerance = ϵ, derivative_tolerance = ϵ, num_additional_partitions = 0)

This function is the entry point for constructing a LP relaxation for the univariate function ``y = f(x)``. The mandatory inputs to the functions are the JuMP model, the JuMP variables x and y, the function (oracle), the partition on ``x``'s domain. All other arguments are optional. The keyword argument `error_tolerance` is used to obtain an LP relaxation by adding more partitions to the `base_partition`. Partitions are added till the vertical distance between the over-estimator and under-estimator of the relaxation is less than the `error_tolerance` for any value in the domain. `length_tolerance` is similar in the sense that it prohibits adding more partitions to an interval that is less than this value and `derivative_tolerance` is used to check if derivaties of the function as successive discretization points are not equal to each other. Finally, the `num_additional_partitions` can be used generate an LP relaxation of the function with a budget on the number of partitions. Note that if the number of partitions is ``n`` and the number of additional partitions is ``m``, then the function will return a relaxation with at most ``n+m`` partitions. 
"""
function construct_univariate_lp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    f::Function,
    x_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_partitions::Int64 = 0,
)::FormulationInfo
    return construct_univariate_lp_relaxation!(
        m,
        x,
        y,
        f,
        x -> ForwardDiff.derivative(f, x),
        x_partition,
        error_tolerance = error_tolerance,
        length_tolerance = length_tolerance,
        derivative_tolerance = derivative_tolerance,
        num_additional_partitions = num_additional_partitions,
    )
end

"""
    construct_univariate_lp_relaxation(m, x, y, f, f_dash, x_partition; error_tolerance = NaN64, length_tolerance = ϵ, derivative_tolerance = ϵ, num_additional_partitions = 0)

This function is the entry point for constructing a LP relaxation for the univariate function ``y = f(x)``. The mandatory inputs to the functions are the JuMP model, the JuMP variables x and y, the function (oracle), the function's derivative (oracle), the partition on ``x``'s domain. All other arguments are optional. The keyword argument `error_tolerance` is used to obtain an LP relaxation by adding more partitions to the `base_partition`. Partitions are added till the vertical distance between the over-estimator and under-estimator of the relaxation is less than the `error_tolerance` for any value in the domain. `length_tolerance` is similar in the sense that it prohibits adding more partitions to an interval that is less than this value and `derivative_tolerance` is used to check if derivaties of the function as successive discretization points are not equal to each other. Finally, the `num_additional_partitions` can be used generate an LP relaxation of the function with a budget on the number of partitions. Note that if the number of partitions is ``n`` and the number of additional partitions is ``m``, then the function will return a relaxation with at most ``n+m`` partitions. 
"""
function construct_univariate_lp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    f::Function,
    f_dash::Function,
    x_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
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
    )
    validate(univariate_function_data)
    validate(x, x_partition)
    refine_partition!(univariate_function_data)
    return build_univariate_lp_relaxation!(m, x, y, univariate_function_data)
end

"get variable bounds from formulation"
@inline get_variable_bounds(formulation::AbstractFormulation) =
    formulation.lower_bounds, formulation.upper_bounds

"get variable names from formulation"
@inline get_variable_names(formulation::AbstractFormulation) = formulation.variable_names

"check if formulation has ``\\geq`` constraints"
@inline has_geq_constraints(formulation::AbstractFormulation) = false

"get the ``\\geq`` constraints from the formulation"
@inline get_geq_constraint_matrices(formulation::AbstractFormulation) =
    Memento.error(_LOGGER, "both the LP and the MILP relaxation have no >= constraints")

"get the maximum vertical distance between the over- and under-estimators across any point in the domain"
@inline get_error_bound(formulation::AbstractFormulation) = formulation.error_bound

"check if the LP formulation has ``=`` constraints"
@inline has_eq_constraints(lp::LPRelaxation) = true

"check if the LP formulation has ``\\leq`` constraints"
@inline has_leq_constraints(lp::LPRelaxation) = false

"get the ``=`` constraints from the LP formulation"
@inline get_eq_constraint_matrices(lp::LPRelaxation) = lp.A, lp.b

"get the ``\\leq`` constraints from the LP formulation"
@inline get_leq_constraint_matrices(lp::LPRelaxation) =
    Memento.error(_LOGGER, "the LP relaxation has no <= constraints")

"get number of variables from the LP formulation"
@inline get_num_variables(lp::LPRelaxation) = length(lp.λ_indices) + 2

"check if the MILP formulation has ``=`` constraints"
@inline has_eq_constraints(milp::MILPRelaxation) = true

"check if the MILP formulation has ``\\leq`` constraints"
@inline has_leq_constraints(milp::MILPRelaxation) = true

"get the ``=`` constraints from the MILP formulation"
@inline get_eq_constraint_matrices(milp::MILPRelaxation) = milp.A_eq, milp.b_eq

"get the ``\\leq`` constraints from the MILP formulation"
@inline get_leq_constraint_matrices(milp::MILPRelaxation) = milp.A_leq, milp.b_leq

"get the variable types (binary/continuous) from the MILP formulation"
@inline get_variable_type(milp::MILPRelaxation) = milp.binary

"get number of binary variables from the MILP formulation"
@inline get_num_binary_variables(milp::MILPRelaxation) = length(milp.z_indices)

"get number of variables from the MILP formulation"
@inline get_num_variables(milp::MILPRelaxation) =
    length(milp.δ_1_indices) + length(milp.δ_2_indices) + length(milp.z_indices) + 2

"get the binary variable indices from the MILP formulation"
@inline get_binary_indices(milp::MILPRelaxation) = findnz(milp.binary)[1]

"get the function from `univariate_function_data`"
@inline get_function(univariate_function_data::UnivariateFunctionData) =
    univariate_function_data.f

"get the derivative from `univariate_function_data`"
@inline get_derivative(univariate_function_data::UnivariateFunctionData) =
    univariate_function_data.f_dash

"get the domain's lower bound from `univariate_function_data`"
@inline get_domain_lb(univariate_function_data::UnivariateFunctionData) =
    univariate_function_data.partition[1]

"get the domain's upper bound from `univariate_function_data`"
@inline get_domain_ub(univariate_function_data::UnivariateFunctionData) =
    univariate_function_data.partition[end]

"get the function's domain from `univariate_function_data`"
@inline get_domain(univariate_function_data::UnivariateFunctionData) =
    get_domain_lb(univariate_function_data), get_domain_ub(univariate_function_data)

"get the partition from `univariate_function_data`"
@inline get_partition(univariate_function_data::UnivariateFunctionData) =
    univariate_function_data.partition
