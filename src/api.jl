function construct_milp_relaxation(f::Function, partition::Vector{<:Real};
        error_tolerance::Float64=Inf,
        num_additional_binary_variables::Int=0)::FormulationData
    return construct_milp_relaxation(
        f,
        x -> ForwardDiff.derivative(f, x),
        partition,
        error_tolerance=error_tolerance,
        num_additional_binary_variables=num_additional_binary_variables)
end

function construct_milp_relaxation(f::Function, d::Function,
        partition::Vector{<:Real}; error_tolerance=Inf,
        num_additional_binary_variables::Int=0)::FormulationData
    return build_formulation(FunctionData(f, d, partition))
end
