function construct_relaxation(uf::UnivariateFunction; 
    partition_points::Vector{Float64}=[], 
    error_tolerance::Float64=Inf,
    binary_variable_budget::Int=0)

    if length(partition_points) == 0 && isinf(error_tolerance) && binary_variable_budget == 0
        return build_model(uf)
    end 

end 