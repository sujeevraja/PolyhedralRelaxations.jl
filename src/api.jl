function construct_milp_relaxation(f::Function, partition::Vector{<:Real};
        error_tolerance::Float64=Inf,
        num_additional_binary_variables::Int=-1)::FormulationData
    return construct_milp_relaxation(
        f,
        x -> ForwardDiff.derivative(f, x),
        partition,
        error_tolerance=error_tolerance,
        num_additional_binary_variables=num_additional_binary_variables)
end

function construct_milp_relaxation(f::Function, d::Function,
        partition::Vector{<:Real}; error_tolerance::Float64=Inf,
        num_additional_binary_variables::Int=0)::FormulationData
    function_data = refine_partition(FunctionData(f, d, partition),
        error_tolerance=error_tolerance,
        num_additional_binary_variables=num_additional_binary_variables)
    return build_formulation(function_data)
end
