struct FunctionData
    f::Function
    d::Function
    base_partition::Vector{Real}
    partition::Vector{Real}  # refinement of `base_partition` based on specified limits
    error_tolerance::Real  # maximum allowed error bound
    length_tolerance::Real  # maximum allowed distance between successive partition points
    derivative_tolerance::Real  # maximum difference between successive derivative values
    num_additional_binary_variables::Int64  # maximum number of additional partition intervals
end

"Getters for FunctionData"
@inline get_function(uf::FunctionData)::Function = uf.f
@inline get_derivative(uf::FunctionData)::Function = uf.d
@inline get_domain_lb(uf::FunctionData)::Real = uf.partition[1]
@inline get_domain_ub(uf::FunctionData)::Real = uf.partition[end]
@inline get_domain(uf::FunctionData)::Pair{Real,Real} = get_domain_lb(uf), get_domain_ub(uf)
@inline get_partition(uf::FunctionData)::Vector{<:Real} = uf.partition

mutable struct ConstraintData
    constraint_row_indices::Vector{Int64}
    constraint_column_indices::Vector{Int64}
    constraint_coefficients::Vector{Real}
    rhs_row_indices::Vector{Int64}
    rhs_values::Vector{Real}
    num_constraints::Int64
end

function ConstraintData()::ConstraintData
    row_indices = Int64[]
    col_indices = Int64[]
    coefs = Real[]
    rhs_row_indices = Int64[]
    rhs_values = Real[]
    num_constraints = 0
    return ConstraintData(row_indices, col_indices, coefs, rhs_row_indices,
        rhs_values, num_constraints)
end

const ConstraintMatrix = Pair{SparseMatrixCSC{Real,Int64},Vector{Real}}

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

"""
Column indices of variables in constraint matrix of polyhedral relaxation.
Indices of δ_1^i, δ_2^i and z_i start from 1.
"""
mutable struct IndexData
    x_index::Int64
    y_index::Int64
    δ_1_indices::Vector{Int64}
    δ_2_indices::Vector{Int64}
    z_indices::Vector{Int64}
end

"""
    IndexData(num_points)

Return the collection of column indices of all variables in the constraint matrix of the polyhedral
relaxation.

If there are k partition points, there are k-1 intervals, with each interval corresponding to a
triangle. As we need one delta_1 variable, delta_2 variable and z variable for each triangle, we
need k-1 of each of these variables in total. For instance, if k=3, we need 2 delta_1 variables with
indices 3,4. As the collect() function includes both endpoints, we should collect only up to
num_vars-1. Then, we will get the correct count for each variable set.
"""
function IndexData(num_points::Int64)::IndexData
    x_index, y_index = 1, 2
    num_vars = num_points - 1

    start = 3
    δ_1_indices = collect(start:(start+num_vars-1))

    start = δ_1_indices[end]+1
    δ_2_indices = collect(start:(start+num_vars-1))

    start = δ_2_indices[end]+1
    z_indices = collect(start:(start+num_vars-1))

    return IndexData(x_index, y_index, δ_1_indices, δ_2_indices, z_indices)
end

"Getters for IndexData"
function get_num_variables(index_data::IndexData)::Int64
    return length(index_data.δ_1_indices) + length(index_data.δ_2_indices) +
        length(index_data.z_indices) + 2
end

"""
Constraint coefficients and right-hand-side of MIP relaxation.

Variables are ordered as: x,y,δ_1^i,δ_2^i,z_i.

All constraints are either equality or less-than-or-equal-to constraints. Row indices of equality
constraints are stored in `equality_row_indices`.
"""
struct FormulationData
    A_eq::SparseMatrixCSC{Real,Int64}
    b_eq::Vector{Real}
    num_eq_constraints::Int64
    A_leq::SparseMatrixCSC{Real,Int64}
    b_leq::Vector{Real}
    num_leq_constraints::Int64
    x_index::Int64
    y_index::Int64
    δ_1_indices::Vector{Int64}
    δ_2_indices::Vector{Int64}
    z_indices::Vector{Int64}
    lower_bounds::Vector{Real}
    upper_bounds::Vector{Real}
    binary::SparseVector{Int64}
    variable_names::Vector{String}
end

function FormulationData(
    function_data::FunctionData,
    index_data::IndexData,
    eq_constraint_data::ConstraintData,
    leq_constraint_data::ConstraintData,
    f_min::Real,
        f_max::Real)::FormulationData
    num_variables = get_num_variables(index_data)
    A_eq, b_eq = get_constraint_matrix(eq_constraint_data, num_variables)
    A_leq, b_leq = get_constraint_matrix(leq_constraint_data, num_variables)

    lower_bounds::Vector{Real} = zeros(num_variables)
    upper_bounds::Vector{Real} = ones(num_variables)

    lower_bounds[index_data.x_index] = function_data.partition[1]
    upper_bounds[index_data.x_index] = function_data.partition[end]
    lower_bounds[index_data.y_index] = f_min
    upper_bounds[index_data.y_index] = f_max
    binary = sparsevec(
        index_data.z_indices,
        ones(length(index_data.z_indices)),
        num_variables)

    variable_names::Vector{String} = ["" for _ in 1:num_variables]
    variable_names[index_data.x_index] = "x"
    variable_names[index_data.y_index] = "y"
    for i in 1:length(index_data.δ_1_indices)
        variable_names[index_data.δ_1_indices[i]] = "delta_1_$i"
        variable_names[index_data.δ_2_indices[i]] = "delta_2_$i"
        variable_names[index_data.z_indices[i]] = "z_$i"
    end

    return FormulationData(
        A_eq, b_eq, eq_constraint_data.num_constraints,
        A_leq, b_leq, leq_constraint_data.num_constraints,
        index_data.x_index,
        index_data.y_index,
        index_data.δ_1_indices,
        index_data.δ_2_indices,
        index_data.z_indices,
        lower_bounds,
        upper_bounds,
        binary,
        variable_names)
end

function get_num_variables(formulation_data::FormulationData)
    return length(formulation_data.δ_1_indices) + length(formulation_data.δ_2_indices) +
        length(formulation_data.z_indices) + 2
end

    struct ConvexHullVariableIndices
        x_index
        y_index
        λ_indices
    end

function ConvexHullVariableIndices(num_λ_variables::Int64)::ConvexHullVariableIndices
    return ConvexHullVariableIndices(1, 2, collect(3:3+num_λ_variables-1))
end

struct ConvexHullFormulation
    A::SparseMatrixCSC{Real,Int64}
    b::Vector{Real}
    num_constraints::Int64
    x_index::Int64
    y_index::Int64
    λ_indices::Vector{Int64}
    lower_bounds::Vector{Real}
    upper_bounds::Vector{Real}
    variable_names::Vector{String}
end

function ConvexHullFormulation(
    function_data::FunctionData,
    constraint_data::ConstraintData,
    index_data::ConvexHullVariableIndices,
    f_min::Real,
        f_max::Real)::ConvexHullFormulation
    num_variables = length(function_data.partition) + 3
    A, b = get_constraint_matrix(constraint_data, num_variables)

    lower_bounds = zeros(num_variables)
    upper_bounds = ones(num_variables)
    lower_bounds[index_data.x_index] = function_data.partition[1]
    upper_bounds[index_data.x_index] = function_data.partition[end]
    lower_bounds[index_data.y_index] = f_min
    upper_bounds[index_data.y_index] = f_max

    variable_names = ["" for _ in 1:num_variables]
    variable_names[index_data.x_index] = "x"
    variable_names[index_data.y_index] = "y"
    for i in 1:length(index_data.λ_indices)
        variable_names[index_data.λ_indices[i]] = "lambda_$i"
    end

    return ConvexHullFormulation(
        A,
        b,
        constraint_data.num_constraints,
        index_data.x_index,
        index_data.y_index,
        index_data.λ_indices,
        lower_bounds,
        upper_bounds,
        variable_names)
end

const Vertex = Pair{Real,Real}
