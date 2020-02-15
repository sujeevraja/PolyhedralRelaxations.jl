"""
    collect_vertices(function_data)

Return a pair of lists with secant vertices as the first element and tangent vertices as the second
element.

Each element in the secant vertex list is a pair (x,y) where
* x is an element of `function_data.partition`,
* y is the value of the given univariate function at x.

Each element in the tangent vertex list is also a pair (x,y). Each position i of the list contains
the vertex formed by intersection of tangents of the curve `y=function_data.f(x)` at
`secant_vertices[i]` and `secant_vertices[i+1]`.

In terms of notation in the paper, `secant_vertices[i]` is the vertex v_i, `secant_vertices[i+1]`
is the vertex v_{i+1} and `tangent_vertices[i]` is the vertex v_{i,i+1}.
"""
function collect_vertices(function_data::FunctionData)::Pair{Vector{Vertex},Vector{Vertex}}
    secant_vertices, tangent_vertices = Vertex[], Vertex[]
    for x in function_data.partition
        push!(secant_vertices, Pair(x, function_data.f(x)))
        if length(secant_vertices) >= 2
            tv = get_tangent_vertex(secant_vertices[end-1], secant_vertices[end], function_data.d)
            push!(tangent_vertices, tv)
        end
    end
    return Pair(secant_vertices, tangent_vertices)
end

"""
    get_tangent_vertex(prev_secant_vertex, next_secant_vertex, derivative)

Return (x,y) coordinates of the intersection of tangents drawn at `prev_secant_vertex` and
`next_secant_vertex`.
"""
function get_tangent_vertex(
    prev_secant_vertex::Vertex,
    next_secant_vertex::Vertex,
        derivative::Function)::Vertex
    x_prev, f_prev = prev_secant_vertex
    x_next, f_next = next_secant_vertex
    d_prev, d_next = derivative(x_prev), derivative(x_next)
    x_t = (f_next - f_prev + (d_prev*x_prev) - (d_next*x_next)) / (d_prev - d_next)
    y_t = f_prev + (d_prev*(x_t - x_prev))
    return Pair(x_t, y_t)
end

"""
    build_formulation(function_data)

Return a FormulationData object with constraint and RHS information of the MILP formulation of the
polyhedral relaxation.
"""
function build_formulation(function_data::FunctionData)::FormulationData
    num_points = length(function_data.partition)
    index_data = IndexData(num_points)
    secant_vertices, tangent_vertices = collect_vertices(function_data)
    constraint_data = ConstraintData()
    add_vertex_constraints!(constraint_data, index_data, secant_vertices, tangent_vertices)
    add_first_delta_constraint!(constraint_data, index_data)
    add_linking_constraints!(constraint_data, index_data, num_points-1)
    return FormulationData(constraint_data, index_data)
end

"""
    add_vertex_constraints!(constraint_data, index_data, secant_vertices,
        tangent_vertices)

Add vertex constraints to `constraint_data` using variable indices from `index_data`.

These constraints link the x and y coordinate variables to the delta variables. The lists
`secant_vertices` and `tangent_vertices` are used to compute coefficients of delta variables.
"""
function add_vertex_constraints!(
    constraint_data::ConstraintData,
    index_data::IndexData,
    secant_vertices::Vector{Vertex},
        tangent_vertices::Vector{Vertex})
    indices = [index_data.x_index, index_data.y_index]
    num_vars = length(secant_vertices) - 1

    for c in [1,2]  # c is the coordinate index (1 for x, 2 for y).
        row = constraint_data.num_constraints+1

        # Add coordinate variable to constraint.
        add_coeff!(constraint_data, row, indices[c], 1)

        for i in 1:num_vars
            # Add delta_1 variable to constraint.
            column = index_data.delta_1_indices[i]
            value =  secant_vertices[i][c] - tangent_vertices[i][c]
            add_coeff!(constraint_data, row, column, value)

            # Add delta_2 variable to constraint.
            column = index_data.delta_2_indices[i]
            value =  secant_vertices[i][c] - secant_vertices[i+1][c]
            add_coeff!(constraint_data, row, column, value)
        end

        # Complete the constraint.
        push!(constraint_data.equality_row_indices, row)
        add_rhs!(constraint_data, row, secant_vertices[1][c])
        constraint_data.num_constraints += 1
    end

    Memento.info(_LOGGER, "built vertex constraints.")
end

"""
    add_first_delta_constraint!(constraint_data, index_data)

Add the constraint "delta_1^1 + delta_2^1 <= 1 to `constraint_data` using variable indices from
`index_data`.
"""
function add_first_delta_constraint!(constraint_data::ConstraintData, index_data::IndexData)
    row = constraint_data.num_constraints + 1
    add_coeff!(constraint_data, row, index_data.delta_1_indices[1], 1)
    add_coeff!(constraint_data, row, index_data.delta_2_indices[1], 1)
    add_rhs!(constraint_data, row, 1)
    constraint_data.num_constraints += 1
    Memento.info(_LOGGER, "built first delta constraint.")
end

"""
    add_linking_constraint!(constraint_data, index_data, num_vars)

Add the constraint families

    delta_1^i + delta_2^i - z_{i-1} <= 0
    delta_2^{i-1} >= z_{i-1}

to `constraint_data` using variable indices from `index_data`. The number of each of these
constraints corresponds to the number of triangles specified by `num_triangles`.
"""
function add_linking_constraints!(
    constraint_data::ConstraintData,
    index_data::IndexData,
        num_triangles::Int64)
    for i in 2:num_triangles
        constraint_data.num_constraints += 1
        row = constraint_data.num_constraints

        # Add delta_1^i + delta_2^i - z_{i-1} <= 0 constraint.
        add_coeff!(constraint_data, row, index_data.delta_1_indices[i], 1)
        add_coeff!(constraint_data, row, index_data.delta_2_indices[i], 1)
        add_coeff!(constraint_data, row, index_data.z_indices[i-1], -1)
        add_rhs!(constraint_data, row, 0)

        # Add z_{i-1} - delta_2^{i-1} <= 0 constraint.
        constraint_data.num_constraints += 1
        row = constraint_data.num_constraints
        add_coeff!(constraint_data, row, index_data.z_indices[i-1], 1)
        add_coeff!(constraint_data, row, index_data.delta_2_indices[i-1], -1)
        add_rhs!(constraint_data, row, 0)
    end
    Memento.info(_LOGGER, "added linking constraints.")
end

"""
    add_coeff!(constraint_data, row, col, value)

Add the coefficient `value` of the variable with index `col` to the constraint with index `row` to
`constraint_data`.
"""
function add_coeff!(constraint_data::ConstraintData, row::Int64, col::Int64, value::Real)
    push!(constraint_data.constraint_row_indices, row)
    push!(constraint_data.constraint_column_indices, col)
    push!(constraint_data.constraint_coefficients, value)
end

"""
    add_rhs!(constraint_data, row, value)

Add the right-hand-side `value` for row `row` to `constraint_data`.
"""
function add_rhs!(constraint_data::ConstraintData, row::Int64, value::Real)
    push!(constraint_data.rhs_row_indices, row)
    push!(constraint_data.rhs_values, value)
end

"""
    main()

Generate model data for the polyhedral relaxation of a univariate function.
"""
function main()
    f = x -> x^3
    partition = Vector{Real}(collect(-1.0:1.0:1.0))
    return construct_milp_relaxation(f, partition)
end
