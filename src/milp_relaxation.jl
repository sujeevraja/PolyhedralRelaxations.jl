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

"Getters for the LP relaxation struct"
@inline has_eq_constraints(milp::MILPRelaxation)::Bool = true
@inline has_leq_constraints(milp::MILPRelaxation)::Bool = true
@inline get_equality_constraint_matrices(milp::MILPRelaxation)::Tuple{SparseMatrixCSC{<:Real,Int64}, Vector{<:Real}} = milp.A_eq, milp.b_eq 
@inline get_leq_constraint_matrices(milp::MILPRelaxation)::Tuple{SparseMatrixCSC{<:Real,Int64}, Vector{<:Real}} = milp.A_leq, milp.b_leq 
@inline get_variable_type(milp::MILPRelaxation)::SparseVector{Int64} = milp.binary
@inline get_num_variables(milp::MILPRelaxation)::Int64 = 
    length(milp.δ_1_indices) + length(milp.δ_2_indices) + length(milp.z_indices) + 2

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

"Getters for MILPVariableIndices"
@inline get_num_variables(milp_variable_indices::MILPVariableIndices)::Int64 = 
    length(milp_variable_indices.δ_1_indices) + length(milp_variable_indices.δ_2_indices) + 
    length(milp_variable_indices.z_indices) + 2

"""
    MILPRelaxation(
        function_data::FunctionData,
        milp_variable_indices::MILPVariableIndices,
        eq_constraint_data::ConstraintData,
        leq_constraint_data::ConstraintData,
        f_min::Real, f_max::Real)::MILPRelaxation
Constructor for the struct MILPRelaxation
"""
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
    add_vertex_constraints!(constraint_data, milp_variable_indices, secant_vertices, tangent_vertices)

Add vertex constraints to `constraint_data` using variable indices from `milp_variable_indices`.

These constraints link the x and y coordinate variables to the δ variables. The lists
`secant_vertices` and `tangent_vertices` are used to compute coefficients of δ variables.
"""
function add_vertex_constraints!(
    constraint_data::ConstraintData,
    milp_variable_indices::MILPVariableIndices,
    secant_vertices::Vector{Vertex},
        tangent_vertices::Vector{Vertex})
    indices = [milp_variable_indices.x_index, milp_variable_indices.y_index]
    num_vars = length(secant_vertices) - 1

    for c in [1,2]  # c is the coordinate index (1 for x, 2 for y).
        row = constraint_data.num_constraints+1

        # Add coordinate variable to constraint.
        add_coeff!(constraint_data, row, indices[c], 1)

        for i in 1:num_vars
            # Add δ_1 variable to constraint.
            column = milp_variable_indices.δ_1_indices[i]
            value =  secant_vertices[i][c] - tangent_vertices[i][c]
            add_coeff!(constraint_data, row, column, value)

            # Add δ_2 variable to constraint.
            column = milp_variable_indices.δ_2_indices[i]
            value =  secant_vertices[i][c] - secant_vertices[i+1][c]
            add_coeff!(constraint_data, row, column, value)
        end

        add_rhs!(constraint_data, row, secant_vertices[1][c])
        constraint_data.num_constraints += 1
    end

    Memento.debug(_LOGGER, "built vertex constraints.")
end

"""
    add_first_δ_constraint!(constraint_data, milp_variable_indices)

Add the constraint "δ_1^1 + δ_2^1 <= 1 to `constraint_data` using variable indices from
`milp_variable_indices`.
"""
function add_first_δ_constraint!(constraint_data::ConstraintData, milp_variable_indices::MILPVariableIndices)
    row = constraint_data.num_constraints + 1
    add_coeff!(constraint_data, row, milp_variable_indices.δ_1_indices[1], 1)
    add_coeff!(constraint_data, row, milp_variable_indices.δ_2_indices[1], 1)
    add_rhs!(constraint_data, row, 1)
    constraint_data.num_constraints += 1
    Memento.debug(_LOGGER, "built first δ constraint.")
end

"""
    add_linking_constraint!(constraint_data, milp_variable_indices, num_vars)

Add the constraint families

    δ_1^i + δ_2^i - z_{i-1} <= 0
    δ_2^{i-1} >= z_{i-1}

to `constraint_data` using variable indices from `milp_variable_indices`. The number of each of these
constraints corresponds to the number of triangles specified by `num_triangles`.
"""
function add_linking_constraints!(
    constraint_data::ConstraintData,
    milp_variable_indices::MILPVariableIndices,
        num_triangles::Int64)
    for i in 2:num_triangles
        constraint_data.num_constraints += 1
        row = constraint_data.num_constraints

        # Add δ_1^i + δ_2^i - z_{i-1} <= 0 constraint.
        add_coeff!(constraint_data, row, milp_variable_indices.δ_1_indices[i], 1)
        add_coeff!(constraint_data, row, milp_variable_indices.δ_2_indices[i], 1)
        add_coeff!(constraint_data, row, milp_variable_indices.z_indices[i-1], -1)

        # Add z_{i-1} - δ_2^{i-1} <= 0 constraint.
        constraint_data.num_constraints += 1
        row = constraint_data.num_constraints
        add_coeff!(constraint_data, row, milp_variable_indices.z_indices[i-1], 1)
        add_coeff!(constraint_data, row, milp_variable_indices.δ_2_indices[i-1], -1)
    end
    Memento.debug(_LOGGER, "added linking constraints.")
end

"""
    build_milp_relaxation(function_data)

Return a MILPRelaxation object with constraint and RHS information of the MILP formulation of the
polyhedral relaxation.
"""
function build_milp_relaxation(function_data::FunctionData)::Pair{MILPRelaxation, FunctionData}
    num_points = length(function_data.partition)
    milp_variable_indices = MILPVariableIndices(num_points)
    secant_vertices, tangent_vertices = collect_vertices(function_data)
    f_min, f_max = Inf, -Inf
    for sv in secant_vertices
        f_val = sv[2]
        f_min = min(f_min, f_val)
        f_max = max(f_max, f_val)
    end
    for tv in tangent_vertices
        f_val = tv[2]
        f_min = min(f_min, f_val)
        f_max = max(f_max, f_val)
    end

    eq_constraint_data = ConstraintData()
    add_vertex_constraints!(eq_constraint_data, milp_variable_indices, secant_vertices, tangent_vertices)

    leq_constraint_data = ConstraintData()
    add_first_δ_constraint!(leq_constraint_data, milp_variable_indices)
    add_linking_constraints!(leq_constraint_data, milp_variable_indices, num_points-1)

    milp_relaxation = MILPRelaxation(function_data, milp_variable_indices, eq_constraint_data,
        leq_constraint_data, f_min, f_max)

    return Pair(milp_relaxation, function_data)
end