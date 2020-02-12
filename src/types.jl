struct FunctionData
    f::Function
    f_dash::Function
    partition::Vector{Real}
end

"Constructor without definition of derivative"
FunctionData(f::Function, partition::Vector{Real}) =
    FunctionData(f, x -> ForwardDiff.derivative(f, x), partition)

"Getters for UnivariateFunction"
@inline get_function(uf::FunctionData)::Function = uf.f
@inline get_derivative(uf::FunctionData)::Function = uf.f_dash
@inline get_domain_lb(uf::FunctionData)::Real = uf.partition[1]
@inline get_domain_ub(uf::FunctionData)::Real = uf.partition[end]
@inline get_domain(uf::FunctionData)::Pair{Real,Real} =
    get_domain_lb(uf), get_domain_ub(uf)
@inline get_partition(uf::FunctionData)::Vector{<:Real} = uf.partition

const Vertex = Pair{Real,Real}

"""
Constraint coefficients and right-hand-side of MIP relaxation.

Variables are ordered as: x, y, delta_1^i, delta_2^i, z_i.

All constraints are either equality or less-than-or-equal-to constraints.
Row indices of equality constraints are stored in `equality_row_indices`.
"""
struct Model
    A::SparseArrays.SparseMatrixCSC{Real,Int64}
    b::SparseArrays.SparseVector{Real,Int64}
    x_index::Int64
    y_index::Int64
    delta_1_indices::Vector{Int64}
    delta_2_indices::Vector{Int64}
    z_indices::Vector{Int64}
    equality_row_indices::Set{Int64}
    num_constraints::Int64
end

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
Indices to recover variable values from model. Indices of delta_1^i, delta_2^i
and z_i start from 1.
"""
mutable struct IndexData
    x_index::Int64
    y_index::Int64
    delta_1_indices::Vector{Int64}
    delta_2_indices::Vector{Int64}
    z_indices::Vector{Int64}
end

function IndexData(num_points::Int64)::IndexData
    x_index = 1
    y_index = 2

    # If there are k partition points, there are k-1 intervals, with each
    # interval corresponding to a triangle. As we need one delta_1 variable,
    # delta_2 variable and z variable for each triangle, we need k-1 of
    # each of these variables in total. For instance, if k=3, we need 2 delta_1
    # variables with indices 3,4. As the collect() function includes both
    # endpoints, we need to set num_vars to k-2. Then, we will get the correct
    # count for each variable set.

    num_vars = num_points - 2

    start = 3
    delta_1_indices = collect(start:(start+num_vars))

    start = delta_1_indices[end]+1
    delta_2_indices = collect(start:(start+num_vars))

    start = delta_2_indices[end]+1
    z_indices = collect(start:(start+num_vars))

    return IndexData(x_index, y_index, delta_1_indices, delta_2_indices,
        z_indices)
end
