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
    equality_row_indices::Set{Int64}
    num_constraints::Int64
end

function ConstraintData()::ConstraintData
    row_indices = Int64[]
    col_indices = Int64[]
    coefs = Real[]
    rhs_row_indices = Int64[]
    rhs_values = Real[]
    equality_row_indices = Set{Int64}()
    num_constraints = 0
    return ConstraintData(row_indices, col_indices, coefs, rhs_row_indices,
        rhs_values, equality_row_indices, num_constraints)
end

"""
Column indices of variables in constraint matrix of polyhedral relaxation.
Indices of delta_1^i, delta_2^i and z_i start from 1.
"""
mutable struct IndexData
    x_index::Int64
    y_index::Int64
    delta_1_indices::Vector{Int64}
    delta_2_indices::Vector{Int64}
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
    delta_1_indices = collect(start:(start+num_vars-1))

    start = delta_1_indices[end]+1
    delta_2_indices = collect(start:(start+num_vars-1))

    start = delta_2_indices[end]+1
    z_indices = collect(start:(start+num_vars-1))

    return IndexData(x_index, y_index, delta_1_indices, delta_2_indices, z_indices)
end

"Getters for IndexData"
function get_num_variables(index_data::IndexData)::Int64
    return length(index_data.delta_1_indices) + length(index_data.delta_2_indices) +
        length(index_data.z_indices) + 2
end

"""
Constraint coefficients and right-hand-side of MIP relaxation.

Variables are ordered as: x, y, delta_1^i, delta_2^i, z_i.

All constraints are either equality or less-than-or-equal-to constraints. Row indices of equality
constraints are stored in `equality_row_indices`.
"""
struct FormulationData
    A::SparseArrays.SparseMatrixCSC{Real,Int64}
    b::SparseArrays.SparseVector{Real,Int64}
    x_index::Int64
    y_index::Int64
    delta_1_indices::Vector{Int64}
    delta_2_indices::Vector{Int64}
    z_indices::Vector{Int64}
    lower_bounds::Vector{Real}
    upper_bounds::Vector{Real}
    binary::SparseArrays.SparseVector{Int64}
    equality_row_indices::SparseArrays.SparseVector{Int64}
    num_constraints::Int64
end

function FormulationData(
    function_data::FunctionData,
    constraint_data::ConstraintData,
        index_data::IndexData)::FormulationData
    A = SparseArrays.sparse(constraint_data.constraint_row_indices,
        constraint_data.constraint_column_indices,
        constraint_data.constraint_coefficients)
    b = SparseArrays.sparsevec(constraint_data.rhs_row_indices, constraint_data.rhs_values)

    num_variables = get_num_variables(index_data)
    lower_bounds::Vector{Real} = zeros(num_variables)
    upper_bounds::Vector{Real} = ones(num_variables)

    lower_bounds[index_data.x_index] = function_data.partition[1]
    upper_bounds[index_data.x_index] = function_data.partition[end]
    lower_bounds[index_data.y_index] = -Inf
    upper_bounds[index_data.y_index] = Inf
    binary = SparseArrays.sparsevec(index_data.z_indices, ones(length(index_data.z_indices)))

    eq_row_index_vec = Vector{Real}()
    append!(eq_row_index_vec, constraint_data.equality_row_indices)
    equality_row_indices = SparseArrays.sparsevec(eq_row_index_vec, ones(length(eq_row_index_vec)))

    return FormulationData(A, b,
        index_data.x_index,
        index_data.y_index,
        index_data.delta_1_indices,
        index_data.delta_2_indices,
        index_data.z_indices,
        lower_bounds,
        upper_bounds,
        binary,
        equality_row_indices,
        constraint_data.num_constraints)
end

function get_num_variables(formulation_data::FormulationData)
    return length(formulation_data.delta_1_indices) + length(formulation_data.delta_2_indices) +
        length(formulation_data.z_indices) + 2
end

const Vertex = Pair{Real,Real}
