export 
    construct_lp_relaxation, 
    construct_milp_relaxation

function construct_milp_relaxation(
    f::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64=NaN64,
    length_tolerance::Float64=ϵ,
    derivative_tolerance::Float64=ϵ,
        num_additional_binary_variables::Int64=0)::Pair{MILPRelaxation, FunctionData}
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
    f_dash::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64=NaN64,
    length_tolerance::Float64=ϵ,
    derivative_tolerance::Float64=ϵ,
        num_additional_binary_variables::Int64=0)::Pair{MILPRelaxation,FunctionData}
    function_data = FunctionData(
        f,
        f_dash,
        base_partition,
        copy(base_partition),
        error_tolerance,
        length_tolerance,
        derivative_tolerance,
        num_additional_binary_variables)
    validate(function_data)
    refine_partition!(function_data)
    return build_milp_relaxation(function_data)
end

function construct_lp_relaxation(
    f::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64=NaN64,
    length_tolerance::Float64=ϵ,
    derivative_tolerance::Float64=ϵ,
        num_additional_binary_variables::Int64=0)::Pair{LPRelaxation,FunctionData}
    return construct_lp_relaxation(
        f,
        x -> ForwardDiff.derivative(f, x),
        base_partition,
        error_tolerance=error_tolerance,
        length_tolerance=length_tolerance,
        derivative_tolerance=derivative_tolerance,
        num_additional_binary_variables=num_additional_binary_variables)
end

function construct_lp_relaxation(
    f::Function,
    f_dash::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Float64=NaN64,
    length_tolerance::Float64=ϵ,
    derivative_tolerance::Float64=ϵ,
        num_additional_binary_variables::Int64=0)::Pair{LPRelaxation,FunctionData}
    function_data = FunctionData(
        f,
        f_dash,
        base_partition,
        copy(base_partition),
        error_tolerance,
        length_tolerance,
        derivative_tolerance,
        num_additional_binary_variables)
    validate(function_data)
    refine_partition!(function_data)
    return build_lp_relaxation(function_data)
end
