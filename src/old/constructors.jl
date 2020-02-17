"Constructor for the struct MILPRelaxation"
function MILPRelaxation(
    function_data::FunctionData,
    milp_variable_indices::MILPVariableIndices,
    eq_constraint_data::ConstraintData,
    leq_constraint_data::ConstraintData,
    f_min::Real, f_max::Real)::MILPRelaxation

    num_variables = get_num_variables(milp_variable_indices)
    A_eq, b_eq = get_constraint_matrix(eq_constraint_data, num_variables)
    A_leq, b_leq = get_constraint_matrix(leq_constraint_data, num_variables)

    lower_bounds::Vector{Real} = zeros(num_variables)
    upper_bounds::Vector{Real} = ones(num_variables)

    lower_bounds[milp_variable_indices.x_index] = function_data.partition[1]
    upper_bounds[milp_variable_indices.x_index] = function_data.partition[end]
    lower_bounds[milp_variable_indices.y_index] = f_min
    upper_bounds[milp_variable_indices.y_index] = f_max
    binary = sparsevec(milp_variable_indices.z_indices,
        ones(length(milp_variable_indices.z_indices)),
        num_variables)

    variable_names::Vector{String} = ["" for _ in 1:num_variables]
    variable_names[milp_variable_indices.x_index] = "x"
    variable_names[milp_variable_indices.y_index] = "y"
    for i in 1:length(milp_variable_indices.δ_1_indices)
        variable_names[milp_variable_indices.δ_1_indices[i]] = "delta_1_$i"
        variable_names[milp_variable_indices.δ_2_indices[i]] = "delta_2_$i"
        variable_names[milp_variable_indices.z_indices[i]] = "z_$i"
    end

    return MILPRelaxation(
        A_eq, b_eq, eq_constraint_data.num_constraints,
        A_leq, b_leq, leq_constraint_data.num_constraints,
        milp_variable_indices.x_index,
        milp_variable_indices.y_index,
        milp_variable_indices.δ_1_indices,
        milp_variable_indices.δ_2_indices,
        milp_variable_indices.z_indices,
        lower_bounds,
        upper_bounds,
        binary,
        variable_names)
end

"Construct for the struct LPRelaxation"
function LPRelaxation(
    function_data::FunctionData,
    constraint_data::ConstraintData,
    lp_variable_indices::LPVariableIndices,
    f_min::Real, f_max::Real)::LPRelaxation 

    num_variables = length(function_data.partition) + 3
    A, b = get_constraint_matrix(constraint_data, num_variables)

    lower_bounds = zeros(num_variables)
    upper_bounds = ones(num_variables)
    lower_bounds[lp_variable_indices.x_index] = function_data.partition[1]
    upper_bounds[lp_variable_indices.x_index] = function_data.partition[end]
    lower_bounds[lp_variable_indices.y_index] = f_min
    upper_bounds[lp_variable_indices.y_index] = f_max

    variable_names = ["" for _ in 1:num_variables]
    variable_names[lp_variable_indices.x_index] = "x"
    variable_names[lp_variable_indices.y_index] = "y"
    for i in 1:length(lp_variable_indices.λ_indices)
        variable_names[lp_variable_indices.λ_indices[i]] = "lambda_$i"
    end

    return LPRelaxation(
        A, b,
        constraint_data.num_constraints,
        lp_variable_indices.x_index,
        lp_variable_indices.y_index,
        lp_variable_indices.λ_indices,
        lower_bounds,
        upper_bounds,
        variable_names)
end

"""
    MILPVariableIndices(num_partition_points::Int64)::MILPVariableIndices

Constructor for the struct MILPVariableIndices; the only input it takes is the number of partition points.
Return the collection of column indices of all variables in the MILP constraint matrix.

If there are k partition points, there are k-1 intervals, with each interval corresponding to a
triangle. As we need one delta_1 variable, delta_2 variable and z variable for each triangle, we
need k-1 of each of these variables in total. For instance, if k=3, we need 2 delta_1 variables with
indices 3,4. As the collect() function includes both endpoints, we should collect only up to
num_vars-1. Then, we will get the correct count for each variable set.
"""
function MILPVariableIndices(num_partition_points::Int64)::MILPVariableIndices
    x_index, y_index = 1, 2
    num_vars = num_partition_points - 1

    start = 3
    δ_1_indices = collect(start:(start+num_vars-1))

    start = δ_1_indices[end]+1
    δ_2_indices = collect(start:(start+num_vars-1))

    start = δ_2_indices[end]+1
    z_indices = collect(start:(start+num_vars-1))

    return MILPVariableIndices(x_index, y_index, δ_1_indices, δ_2_indices, z_indices)
end

"""
    LPVariableIndices(num_partition_points::Int64)::LPVariableIndices

Constructor for the LPVariableIndices struct. The only input it takes is the number of partition points.
"""
function LPVariableIndices(num_partition_points::Int64)::LPVariableIndices
    num_λ_variables = num_partition_points + 1
    return LPVariableIndices(1, 2, collect(3:3+num_λ_variables-1))
end
