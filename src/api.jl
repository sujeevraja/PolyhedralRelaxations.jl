export construct_lp_relaxation,
    construct_milp_relaxation,
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
    construct_milp_relaxation(f, base_partition; error_tolerance = NaN64, length_tolerance = ϵ, derivative_tolerance = ϵ, num_additional_partitions = 0)

This function is the entry point for constructing a MILP relaxation for the constraint ``y = f(x)``. The mandatory inputs to the functions are the function, the base partition. All other arguments are optional. The keyword argument `error_tolerance` is used to obtain an MILP relaxation by adding more partitions to the `base_partition`. Partitions are added till the vertical distance between the over-estimator and under-estimator of the relaxation is less than the `error_tolerance` for any value in the domain. `length_tolerance` is similar in the sense that it prohibits adding more partitions to an interval that is less than this value and `derivative_tolerance` is used to check if derivaties of the function as successive discretization points are not equal to each other. Finally, the `num_additional_partitions` can be used generate an LP relaxation of the function with a budget on the number of partitions. Note that if the number of partitions is ``n`` and the number of additional partitions is ``m``, then the function will return a relaxation with at most ``n+m`` partitions. 
"""
function construct_milp_relaxation(
    f::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_partitions::Int64 = 0,
)::Pair{MILPRelaxation,FunctionData}
    return construct_milp_relaxation(
        f,
        x -> ForwardDiff.derivative(f, x),
        base_partition,
        error_tolerance = error_tolerance,
        length_tolerance = length_tolerance,
        derivative_tolerance = derivative_tolerance,
        num_additional_partitions = num_additional_partitions,
    )
end

"""
    construct_milp_relaxation(f, f_dash, base_partition; error_tolerance = NaN64, length_tolerance = ϵ, derivative_tolerance = ϵ, num_additional_partitions = 0)
    
This function is the entry point for constructing a MILP relaxation for the constraint ``y = f(x)``. The mandatory inputs to the functions are the function, its derivative, and the base partition. All other arguments are optional. The keyword argument `error_tolerance` is used to obtain an MILP relaxation by adding more partitions to the `base_partition`. Partitions are added till the vertical distance between the over-estimator and under-estimator of the relaxation is less than the `error_tolerance` for any value in the domain. `length_tolerance` is similar in the sense that it prohibits adding more partitions to an interval that is less than this value and `derivative_tolerance` is used to check if derivaties of the function as successive discretization points are not equal to each other. Finally, the `num_additional_partitions` can be used generate an LP relaxation of the function with a budget on the number of partitions. Note that if the number of partitions is ``n`` and the number of additional partitions is ``m``, then the function will return a relaxation with at most ``n+m`` partitions. 
"""
function construct_milp_relaxation(
    f::Function,
    f_dash::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_partitions::Int64 = 0,
)::Pair{MILPRelaxation,FunctionData}
    function_data = FunctionData(
        f,
        f_dash,
        base_partition,
        copy(base_partition),
        error_tolerance,
        length_tolerance,
        derivative_tolerance,
        num_additional_partitions,
    )
    validate(function_data)
    refine_partition!(function_data)
    return build_milp_relaxation(function_data)
end

"""
    construct_lp_relaxation(f, base_partition; error_tolerance = NaN64, length_tolerance = ϵ, derivative_tolerance = ϵ, num_additional_partitions = 0)

This function is the entry point for constructing a LP relaxation for the constraint ``y = f(x)``. The mandatory inputs to the functions are the function, the base partition. All other arguments are optional. The keyword argument `error_tolerance` is used to obtain an LP relaxation by adding more partitions to the `base_partition`. Partitions are added till the vertical distance between the over-estimator and under-estimator of the relaxation is less than the `error_tolerance` for any value in the domain. `length_tolerance` is similar in the sense that it prohibits adding more partitions to an interval that is less than this value and `derivative_tolerance` is used to check if derivaties of the function as successive discretization points are not equal to each other. Finally, the `num_additional_partitions` can be used generate an LP relaxation of the function with a budget on the number of partitions. Note that if the number of partitions is ``n`` and the number of additional partitions is ``m``, then the function will return a relaxation with at most ``n+m`` partitions. 
"""
function construct_lp_relaxation(
    f::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_partitions::Int64 = 0,
)::Pair{LPRelaxation,FunctionData}
    return construct_lp_relaxation(
        f,
        x -> ForwardDiff.derivative(f, x),
        base_partition,
        error_tolerance = error_tolerance,
        length_tolerance = length_tolerance,
        derivative_tolerance = derivative_tolerance,
        num_additional_partitions = num_additional_partitions,
    )
end

"""
    construct_lp_relaxation(f, f_dash, base_partition; error_tolerance = NaN64, length_tolerance = ϵ, derivative_tolerance = ϵ, num_additional_partitions = 0)

This function is the entry point for constructing a LP relaxation for the constraint ``y = f(x)``. The mandatory inputs to the functions are the function, the derivate, and the base partition. All other arguments are optional. The keyword argument `error_tolerance` is used to obtain an LP relaxation by adding more partitions to the `base_partition`. Partitions are added till the vertical distance between the over-estimator and under-estimator of the relaxation is less than the `error_tolerance` for any value in the domain. `length_tolerance` is similar in the sense that it prohibits adding more partitions to an interval that is less than this value and `derivative_tolerance` is used to check if derivaties of the function as successive discretization points are not equal to each other. Finally, the `num_additional_partitions` can be used generate an LP relaxation of the function with a budget on the number of partitions. Note that if the number of partitions is ``n`` and the number of additional partitions is ``m``, then the function will return a relaxation with at most ``n+m`` partitions. 
"""
function construct_lp_relaxation(
    f::Function,
    f_dash::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64 = NaN64,
    length_tolerance::Float64 = ϵ,
    derivative_tolerance::Float64 = ϵ,
    num_additional_partitions::Int64 = 0,
)::Pair{LPRelaxation,FunctionData}
    function_data = FunctionData(
        f,
        f_dash,
        base_partition,
        copy(base_partition),
        error_tolerance,
        length_tolerance,
        derivative_tolerance,
        num_additional_partitions,
    )
    validate(function_data)
    refine_partition!(function_data)
    return build_lp_relaxation(function_data)
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

"get the function from `function_data`"
@inline get_function(function_data::FunctionData) = function_data.f

"get the derivative from `function_data`"
@inline get_derivative(function_data::FunctionData) = function_data.f_dash

"get the domain's lower bound from `function_data`"
@inline get_domain_lb(function_data::FunctionData) = function_data.partition[1]

"get the domain's upper bound from `function_data`"
@inline get_domain_ub(function_data::FunctionData) = function_data.partition[end]

"get the function's domain from `function_data`"
@inline get_domain(function_data::FunctionData) =
    get_domain_lb(function_data), get_domain_ub(function_data)

"get the partition from `function_data`"
@inline get_partition(function_data::FunctionData) = function_data.partition
