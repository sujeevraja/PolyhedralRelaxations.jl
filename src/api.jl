function construct_milp_relaxation(
    f::Function,
    base_partition::Vector{<:Real};
    error_tolerance::Real=NaN64,  # maximum allowed error bound
    length_tolerance::Real=EPS, # maximum allowed distance between successive partition points
    derivative_tolerance::Real=EPS,  # maximum difference between successive derivative values
        num_additional_binary_variables::Int=0)::FormulationData  # maximum additional added partition intervals
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
    length_tolerance::Real=EPS,
    derivative_tolerance::Real=EPS,
        num_additional_binary_variables::Int=0)::FormulationData
    function_data = FunctionData(f, d, base_partition, copy(base_partition), error_tolerance,
        length_tolerance, derivative_tolerance, num_additional_binary_variables)
    refine_partition!(function_data)
    return build_formulation(function_data)
end
