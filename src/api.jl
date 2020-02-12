function construct_milp_relaxation(f::Function, partition::Vector{Real})::Model
    return construct_milp_relaxation(f, x -> ForwardDiff(f, x), partition)
end

function construct_milp_relaxation(f::Function, d::Function, partition::Vector{Real})::Model
    return build_model(FunctionData(f, d, partition))
end

# function construct_relaxation(function_data::FunctionData;
#     partition_points::Vector{Float64}=[],
#     error_tolerance::Float64=Inf,
#     binary_variable_budget::Int=0)

#     if (length(partition_points) == 0 && isinf(error_tolerance) &&
#         binary_variable_budget == 0)
#         return build_model(uf)
#     end
# end
