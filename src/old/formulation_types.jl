"""
MILP relaxation struct which contains 
Constraint coefficients and right-hand-side of MIP relaxation.

Variables are ordered as: x,y,δ_1^i,δ_2^i,z_i.

All constraints are either equality or less-than-or-equal-to constraints. 
"""
struct MILPRelaxation <: AbstractFormulation 
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

"""
LP relaxation struct which contains 
Constraint coefficients and right-hand-side of LP relaxation.

Variables are ordered as: x,y,λ_i.

All constraints are equality constraints. 
"""
struct LPRelaxation <: AbstractFormulation
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

"""
Column indices of variables in constraint matrix of MILP relaxation.
Indices of δ_1^i, δ_2^i and z_i start from 3.
"""
struct MILPVariableIndices <: AbstractVariableIndices
    x_index::Int64
    y_index::Int64
    δ_1_indices::Vector{Int64}
    δ_2_indices::Vector{Int64}
    z_indices::Vector{Int64}
end

"""
Column indices of variables in constraint matrix of LP relaxation.
"""
struct LPVariableIndices <: AbstractVariableIndices
    x_index
    y_index
    λ_indices
end
