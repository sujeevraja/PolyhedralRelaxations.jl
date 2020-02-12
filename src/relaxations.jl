"""
    collect_vertices(function_data)

Return a pair of lists with secant vertices as the first element and tangent
vertices as the second element. All vertices

Each element in the secant vertex list is a pair (x,y) where
* x is an element of `function_data.partition`,
* y is the value of the given univariate function at x.

Each element in the tangent vertex list is also a pair (x,y). Each position i
of the list contains the vertex formed by intersection of tangents of the curve
`y=function_data.f(x)` at `secant_vertices[i]` and `secant_vertices[i+1]`.

In terms of notation in the paper, `secant_vertices[i]` is the vertex v_i,
`secant_vertices[i+1]` is the vertex v_{i+1} and `tangent_vertices[i]` is
the vertex v_{i,i+1}.
"""
function collect_vertices(function_data::FunctionData)::Pair{Vector{Vertex},Vector{Vertex}}
    validate_partition(function_data.partition)

    secant_vertices = Vertex[]
    tangent_vertices = Vertex[]

    # Add first secant vertex.
    x_prev = function_data.partition[1]
    push!(secant_vertices, get_secant_vertex(x_prev, function_data.f))

    # Add remaining secant vertices and all tangent vertices.
    for i in 2:length(function_data.partition)
        x = function_data.partition[i]
        push!(secant_vertices, get_secant_vertex(x, function_data.f))
        push!(tangent_vertices, get_tangent_vertex(secant_vertices[end-1],
            secant_vertices[end], function_data.d))
    end

    return Pair(secant_vertices, tangent_vertices)
end

function validate_partition(partition::Vector{<:Real})
    # Check partition size.
    num_points = length(partition)
    if num_points < 2
        Memento.error(_LOGGER, "partition must have at least 2 points")
    end

    # Check finiteness of first point.
    if !isfinite(partition[1])
        Memento.error(_LOGGER, "all partition points must be finite")
    end

    for i in 2:num_points
        x_prev = partition[i-1]
        x = partition[i]

        # Check finiteness of current point.
        if !isfinite(x)
            Memento.error(_LOGGER, "all partition points must be finite")
        end

        # Check ascending order of partition.
        if x <= x_prev
            Memento.error(_LOGGER,
                "partition must be sorted, violation for $x, $x_prev")
        end

        # Check length of partition interval.
        if (x - x_prev) <= EPS
            Memento.erro(_LOGGER,
                "partition length from $x_prev to $x shorter than $EPS")
        end
    end

    Memento.info(_LOGGER, "partition points are valid.")
end

"""
    get_secant_vertex(x, f)

Return the pair (x,f(x)) after verifying finiteness of f(x).
"""
function get_secant_vertex(x::Real, f::Function)::Vertex
    fx = f(x)
    if !isfinite(fx)
        Memento.error(_LOGGER,
            "boundedness necessary, function unbounded at $x")
    end

    return Pair(x,fx)
end

"""
    get_tangent_vertex(prev_secant_vertex, next_secant_vertex, derivative)

Return (x,y) coordinates of the intersection of tangents drawn at
`prev_secant_vertex` and `next_secant_vertex`.
"""
function get_tangent_vertex(prev_secant_vertex::Vertex, next_secant_vertex::Vertex, derivative::Function)::Vertex
    x_prev = prev_secant_vertex[1]
    d_prev = derivative(x_prev)
    if !isfinite(d_prev)
        Memento.error(_LOGGER, "derivative unbounded at $x_prev")
    end

    x_next = next_secant_vertex[1]
    d_next = derivative(x_next)
    if !isfinite(d_next)
        Memento.error(_LOGGER, "derivative unbounded at $x_next")
    end

    if isapprox(d_prev, d_next, atol=EPS)
        Memento.error(_LOGGER,
            "$x_prev,$x_next have derivatives with difference less than $EPS")
    end

    f_prev = prev_secant_vertex[2]
    f_next = next_secant_vertex[2]
    x_t = (f_next - f_prev + (d_prev*x_prev) - (d_next*x_next)) / (d_prev - d_next)
    y_t = f_prev + (d_prev*(x_t - x_prev))
    return Pair(x_t, y_t)
end

"""
    build_formulation(function_data)

Return a FormulationData object with constraint and RHS information of the MILP
formulation of the polyhedral relaxation.
"""
function build_formulation(function_data::FunctionData)::FormulationData
    Memento.info(_LOGGER, "starting to build formulation data...")
    secant_vertices, tangent_vertices = collect_vertices(function_data)
    Memento.info(_LOGGER, "got $(length(secant_vertices)) secant vertices.")
    Memento.info(_LOGGER, "got $(length(tangent_vertices)) tangent vertices.")

    # Indices to recover variable values from model. Indices of delta_1^i,
    # delta_2^i and z_i start from 1.
    num_points = length(secant_vertices)
    index_data = IndexData(num_points)
    Memento.info(_LOGGER, "number of partition points: $num_points")
    Memento.info(_LOGGER, "x index: $(index_data.x_index)")
    Memento.info(_LOGGER, "y index: $(index_data.y_index)")

    i_start = index_data.delta_1_indices[1]
    i_end = index_data.delta_1_indices[end]
    Memento.info(_LOGGER, "delta_1_indices: $i_start to $i_end")

    i_start = index_data.delta_2_indices[1]
    i_end = index_data.delta_2_indices[end]
    Memento.info(_LOGGER, "delta_2_indices: $i_start to $i_end")

    i_start = index_data.z_indices[1]
    i_end = index_data.z_indices[end]
    Memento.info(_LOGGER, "z_indices: $i_start to $i_end")

    constraint_data = ConstraintData()
    add_vertex_constraints!(constraint_data, index_data, secant_vertices, tangent_vertices)
    add_first_delta_constraint!(constraint_data, index_data)
    add_linking_constraints!(constraint_data, index_data, num_points-1)

    # Store constraint data into a Model object and return it.
    A = SparseArrays.sparse(constraint_data.constraint_row_indices,
        constraint_data.constraint_column_indices,
        constraint_data.constraint_coefficients)
    b = SparseArrays.sparsevec(constraint_data.rhs_row_indices, constraint_data.rhs_values)
    Memento.info(_LOGGER, "built formulation data.")

    return FormulationData(A, b,
        index_data.x_index,
        index_data.y_index,
        index_data.delta_1_indices,
        index_data.delta_2_indices,
        index_data.z_indices,
        constraint_data.equality_row_indices,
        constraint_data.num_constraints)
end

"""
    add_vertex_constraints!(constraint_data, index_data, secant_vertices,
        tangent_vertices)

Add vertex constraints to `constraint_data` using variable indices from
`index_data`.

These constraints link the x and y coordinate variables to the delta variables.
The lists `secant_vertices` and `tangent_vertices` are used to compute
coefficients of delta variables.
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

Add the constraint "delta_1^1 + delta_2^1 <= 1 to `constraint_data` using
variable indices from `index_data`.
"""
function add_first_delta_constraint!(
        constraint_data::ConstraintData,
        index_data::IndexData)
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

to `constraint_data` using variable indices from `index_data`. The number of
each of these constraints corresponds to the number of triangles specified by
`num_triangles`.
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

Add the coefficient `value` of the variable with index `col` to the constraint
with index `row` to `constraint_data`.
"""
function add_coeff!(
        constraint_data::ConstraintData,
        row::Int64,
        col::Int64,
        value::Real)
    push!(constraint_data.constraint_row_indices, row)
    push!(constraint_data.constraint_column_indices, col)
    push!(constraint_data.constraint_coefficients, value)
end

"""
    add_rhs!(constraint_data, row, value)

Add the right-hand-side `value` for row `row` to `constraint_data`.
"""
function add_rhs!(
        constraint_data::ConstraintData,
        row::Int64,
        value::Real)
    push!(constraint_data.rhs_row_indices, row)
    push!(constraint_data.rhs_values, value)
end

"""
    main()

Generate model data for the polyhedral relaxation of a univariate function.
"""
function main()
    function_data = FunctionData(
        x -> x^3,  # f
        x -> 3*(x^2),  # f'
        Vector{Real}(collect(-1.0:1.0:1.0)))
    return build_formulation(function_data)
end
