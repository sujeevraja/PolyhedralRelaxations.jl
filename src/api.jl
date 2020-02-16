function construct_milp_relaxation(
    f::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Real=NaN64,
    length_tolerance::Real=ϵ,
    derivative_tolerance::Real=ϵ,
        num_additional_binary_variables::Int=0)::Pair{FormulationData, FunctionData}
    return construct_milp_relaxation(
        f,
        x -> ForwardDiff.derivative(f, x),
        base_partition,
        error_tolerance=error_tolerance,
        length_tolerance=length_tolerance,
        derivative_tolerance=derivative_tolerance,
        num_additional_binary_variables=num_additional_binary_variables)
end

function construct_milp_relaxation(
    f::Function,
    d::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Real=NaN64,
    length_tolerance::Real=ϵ,
    derivative_tolerance::Real=ϵ,
        num_additional_binary_variables::Int=0)::Pair{FormulationData,FunctionData}
    function_data = FunctionData(
        f,
        d,
        base_partition,
        copy(base_partition),
        error_tolerance,
        length_tolerance,
        derivative_tolerance,
        num_additional_binary_variables)
    validate(function_data)
    refine_partition!(function_data)
    return build_formulation(function_data)
end

function construct_convex_hull_relaxation(
    f::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Real=NaN64,
    length_tolerance::Real=ϵ,
    derivative_tolerance::Real=ϵ,
        num_additional_binary_variables::Int=0)::Pair{ConvexHullFormulation,FunctionData}
    return construct_convex_hull_relaxation(
        f,
        x -> ForwardDiff.derivative(f, x),
        base_partition,
        error_tolerance=error_tolerance,
        length_tolerance=length_tolerance,
        derivative_tolerance=derivative_tolerance,
        num_additional_binary_variables=num_additional_binary_variables)
end

function construct_convex_hull_relaxation(
    f::Function,
    d::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Real=NaN64,
    length_tolerance::Real=ϵ,
    derivative_tolerance::Real=ϵ,
        num_additional_binary_variables::Int=0)::Pair{ConvexHullFormulation,FunctionData}
    function_data = FunctionData(
        f,
        d,
        base_partition,
        copy(base_partition),
        error_tolerance,
        length_tolerance,
        derivative_tolerance,
        num_additional_binary_variables)
    validate(function_data)
    refine_partition!(function_data)
    return build_convex_hull_formulation(function_data)
end
