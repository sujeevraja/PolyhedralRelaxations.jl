"""
LP relaxation struct which contains
Constraint coefficients and right-hand-side of LP relaxation.

Variables are ordered as: x,y,λ_i.

All constraints are equality constraints.
"""
struct LPRelaxation <: AbstractFormulation
    A::SparseMatrixCSC{Float64,Int64}
    b::Vector{Float64}
    num_constraints::Int64
    x_index::Int64
    y_index::Int64
    λ_indices::Vector{Int64}
    lower_bounds::Vector{Float64}
    upper_bounds::Vector{Float64}
    variable_names::Vector{String}
    error_bound::Float64
end

"""
Column indices of variables in constraint matrix of LP relaxation.
"""
struct LPVariableIndices <: AbstractVariableIndices
    x_index
    y_index
    λ_indices
end

"""
    LPRelaxation(
        function_data::FunctionData,
        constraint_data::ConstraintData,
        lp_variable_indices::LPVariableIndices,
        f_min::Float64,
        f_max::Float64)::LPRelaxation

Constructor for the struct LPRelaxation
"""
function LPRelaxation(
    function_data::FunctionData,
    constraint_data::ConstraintData,
    lp_variable_indices::LPVariableIndices,
    f_min::Float64,
    f_max::Float64,
)::LPRelaxation
    num_variables = length(function_data.partition) + 3
    A, b = get_constraint_matrix(constraint_data, num_variables)

    lower_bounds = zeros(num_variables)
    upper_bounds = ones(num_variables)
    lower_bounds[lp_variable_indices.x_index] = function_data.partition[1]
    upper_bounds[lp_variable_indices.x_index] = function_data.partition[end]
    lower_bounds[lp_variable_indices.y_index] = f_min
    upper_bounds[lp_variable_indices.y_index] = f_max

    variable_names = ["" for _ = 1:num_variables]
    variable_names[lp_variable_indices.x_index] = "x"
    variable_names[lp_variable_indices.y_index] = "y"
    for i = 1:length(lp_variable_indices.λ_indices)
        variable_names[lp_variable_indices.λ_indices[i]] = "lambda_$i"
    end

    return LPRelaxation(
        A,
        b,
        constraint_data.num_constraints,
        lp_variable_indices.x_index,
        lp_variable_indices.y_index,
        lp_variable_indices.λ_indices,
        lower_bounds,
        upper_bounds,
        variable_names,
        get_max_error_bound(function_data),
    )
end

"""
    LPVariableIndices(num_partition_points::Int64)::LPVariableIndices

Constructor for the LPVariableIndices struct. The only input it takes is the number of partition
points.
"""
function LPVariableIndices(num_partition_points::Int64)::LPVariableIndices
    num_λ_variables = num_partition_points + 1
    return LPVariableIndices(1, 2, collect(3:3+num_λ_variables-1))
end

"""
    get_convex_hull_vertices(function_data::FunctionData)::Vector{Vertex}

Returns all the vertices that are a part of the LP relaxation given the function data
"""
function get_lp_relaxation_vertices(function_data::FunctionData)::Vector{Vertex}
    secant_vertices, tangent_vertices = collect_vertices(function_data)
    vertices = Vertex[]
    push!(vertices, secant_vertices[1])
    append!(vertices, tangent_vertices)
    push!(vertices, secant_vertices[end])
    return vertices
end

"""
    build_lp_relaxation(function_data::FunctionData)::Pair{LPRelaxation,FunctionData}

Function to build the LP relaxation given the function data
"""
function build_lp_relaxation(function_data::FunctionData)::Pair{LPRelaxation,FunctionData}
    vertices = get_lp_relaxation_vertices(function_data)
    f_min, f_max = Inf, -Inf
    for v in vertices
        f_val = v[2]
        f_min = min(f_min, f_val)
        f_max = max(f_max, f_val)
    end

    lp_variable_indices = LPVariableIndices(length(function_data.partition))

    constraint_data = ConstraintData()

    # Add convex combination constraint coefficients.
    x_row, y_row, simplex_row = 1, 2, 3
    add_coeff!(constraint_data, x_row, lp_variable_indices.x_index, 1.0)
    add_coeff!(constraint_data, y_row, lp_variable_indices.y_index, 1.0)
    for i = 1:length(vertices)
        add_coeff!(
            constraint_data,
            x_row,
            lp_variable_indices.λ_indices[i],
            float(-vertices[i][1]),
        )
        add_coeff!(
            constraint_data,
            y_row,
            lp_variable_indices.λ_indices[i],
            float(-vertices[i][2]),
        )
        add_coeff!(constraint_data, simplex_row, lp_variable_indices.λ_indices[i], 1.0)
    end
    add_rhs!(constraint_data, x_row, 0.0)
    add_rhs!(constraint_data, y_row, 0.0)
    add_rhs!(constraint_data, simplex_row, 1.0)
    constraint_data.num_constraints = 3

    validate(constraint_data)

    # Build and return LP relaxation with function data.
    lp_relaxation =
        LPRelaxation(function_data, constraint_data, lp_variable_indices, f_min, f_max)
    return Pair(lp_relaxation, function_data)
end
