"""
    get_num_variables(milp_variable_indices::MILPVariableIndices)::Int64 

    Returns the number of variables in the MILP relaxation
"""
function get_num_variables(milp_variable_indices::MILPVariableIndices)::Int64
    return length(milp_variable_indices.δ_1_indices) + 
        length(milp_variable_indices.δ_2_indices) +
        length(milp_variable_indices.z_indices) + 2
end 

"""
    get_num_variables(lp_variable_indices::LPVariableIndices)::Int64 

    Returns the number of variables in the LP relaxation 
"""
get_num_variables(lp_variable_indices::LPVariableIndices)::Int64 = 
    length(lp_variable_indices.λ_indices) + 2

"""
    get_constraint_matrix(constraint_data::ConstraintData, num_columns::Int64)::ConstraintMatrix

    Creates the Pair(A, b) from the constraint_data and number of columns (number of variables).
    Here A is a sparse matrix and b is a dense vector. 
"""
function get_constraint_matrix(constraint_data::ConstraintData, num_columns::Int64)::ConstraintMatrix
    A = sparse(constraint_data.constraint_row_indices,
        constraint_data.constraint_column_indices,
        constraint_data.constraint_coefficients,
        constraint_data.num_constraints, num_columns)
    b = sparsevec(constraint_data.rhs_row_indices,
        constraint_data.rhs_values,
        constraint_data.num_constraints)
    return Pair(A,Vector(b))
end