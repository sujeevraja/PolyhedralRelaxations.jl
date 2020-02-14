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
    # Add first secant vertex.
    secant_vertices = Vertex[]
    x_prev = function_data.partition[1]
    push!(secant_vertices, get_secant_vertex(x_prev, function_data.f))

    # Add remaining secant vertices and all tangent vertices.
    tangent_vertices = Vertex[]
    for i in 2:length(function_data.partition)
        x = function_data.partition[i]
        push!(secant_vertices, get_secant_vertex(x, function_data.f))
        push!(tangent_vertices, get_tangent_vertex(secant_vertices[end-1],
            secant_vertices[end], function_data.d))
    end

    return Pair(secant_vertices, tangent_vertices)
end

function validate_partition(partition::Vector{<:Real}, tolerance::Real)
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
            Memento.error(_LOGGER, "partition must be sorted, violation for $x, $x_prev")
        end

        # Check length of partition interval.
        if (x - x_prev) <= tolerance
            Memento.error(_LOGGER, "$x_prev and $x closer than $(function_data.length_tolerance)")
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
        Memento.error(_LOGGER, "boundedness necessary, function unbounded at $x")
    end
    return Pair(x,fx)
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

    f_prev = prev_secant_vertex[2]
    f_next = next_secant_vertex[2]
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
    Memento.info(_LOGGER, "starting to build formulation data...")
    validate_partition(function_data.partition, function_data.length_tolerance)
    secant_vertices, tangent_vertices = collect_vertices(function_data)
    Memento.info(_LOGGER, "got $(length(secant_vertices)) secant vertices.")
    Memento.info(_LOGGER, "got $(length(tangent_vertices)) tangent vertices.")

    # Indices to recover variable values from model. Indices of delta_1^i, delta_2^i and z_i start
    # from 1.
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

function refine_partition!(function_data::FunctionData)
    # Don't refine the partition if no additional constraints are specified.
    if isnan(function_data.error_tolerance) && function_data.num_additional_binary_variables <= 0
        return
    end

    error_queue = get_error_queue(function_data)
    num_partitions_added = 0
    partition = function_data.partition
    while true
        if !is_refinement_feasible(function_data, error_queue)
            Memento.info(_LOGGER, "stopping refinement")
            return
        end

        start, max_error = peek(error_queue)
        x_start = partition[start]
        x_end = partition[start+1]
        Memento.info(_LOGGER, "max error: $max_error between $x_start, $x_end")

        # Errors of partition intervals in `error_queue` are indexed by positions of interval
        # starts in `refined_partition`. As we will be inserting `x_new` into `refined_partition`
        # between positions `start` and `start+1`, the positions of interval-starts after `x_new`
        # will all increase by 1 after the insertions. Upade the queue with this new indexing.
        num_starts = length(partition)
        for i in num_starts:-1:start+1
            error_queue[i] = error_queue[i-1]
        end

        # Add errors of the new partition intervals to the queue.
        x_new = 0.5 * (x_start + x_end)
        error_queue[start] = get_error_bound(function_data.d, x_start, x_new)
        error_queue[start+1] = get_error_bound(function_data.d, x_new, x_end)

        # Add `x_new` to `refined_partition`.
        splice!(partition, start+1:start, x_new)
        num_partitions_added += 1
    end
end

function is_refinement_feasible(function_data::FunctionData, error_queue::PriorityQueue)::Bool
    # Check if error bound is below tolerance.
    start, max_error = peek(error_queue)
    if !isnan(function_data.error_tolerance) && max_error <= function_data.error_tolerance
        Memento.info(_LOGGER, "error: $max_error less than limit: $(function_data.error_tolerance)")
        return false
    end

    # Check if the maximum allowed number of binary variables have been added.
    if function_data.num_additional_binary_variables > 0
        num_added = length(function_data.partition) - length(function_data.base_partition)
        if num_added >= function_data.num_additional_binary_variables
            Memento.info(_LOGGER, "number of new binary variables: $num_added")
            Memento.info(_LOGGER, "budget: $(function_data.num_additional_binary_variables)")
            return false
        end
    end

    # Check if new partition lengths are smaller than allowed.
    x_start, x_end = function_data.partition[start], function_data.partition[start+1]
    if x_end - x_start <= (2 * function_data.length_tolerance)
        Memento.info(_LOGGER, "start: $x_start, end: $x_end")
        Memento.info(_LOGGER, "new interval length will be too small")
        return false
    end

    # Check if derivative differences at endpoints of new partitions are smaller than allowed.
    d_start, d_end = function_data.d(x_start), function_data.d(x_end)
    x_new = 0.5 * (x_start + x_end)
    d_new = function_data.d(x_new)
    if (abs(d_new - d_start) <= function_data.derivative_tolerance ||
        abs(d_end - d_new) <= function_data.derivative_tolerance)
        Memento.info(_LOGGER, "d_start: $d_start, d_new: $d_new, d_end: $d_end")
        Memento.info(_LOGGER, "adjacent derivative difference will be too small")
        return false
    end

    return true
end

"""
    get_error_queue(function_data)

Build a max-priority-queue holding errors of each partition interval in `function_data.partition`.
Keys of the queue are indices of starting positions of the partition interval in
`function_data.partition`. Priorities are error bounds of partition intervals. The queue is built as
a max-queue for easy access to the  maximum error.
"""
function get_error_queue(function_data::FunctionData)::PriorityQueue
    pq = PriorityQueue{Int64,Float64}(Base.Order.Reverse)
    for i in 1:length(function_data.partition)-1
        lb = function_data.partition[i]
        ub = function_data.partition[i+1]
        pq[i] = get_error_bound(function_data.d, lb, ub)
    end
    return pq
end

"""
    get_error_bound(derivative, lb, ub)

Get error bound of a function with derivative `derivative` in the closed interval `[lb,ub]`.
"""
function get_error_bound(derivative::Function, lb::Real, ub::Real)
    return (ub-lb) * abs(derivative(ub) - derivative(lb)) / 4.0
end

"""
    main()

Generate model data for the polyhedral relaxation of a univariate function.
"""
function main()
    f = x -> x^3
    partition = Vector{Real}(collect(-1.0:1.0:1.0))
    return construct_milp_relaxation(f, partition, error_tolerance=0.01,
        num_additional_binary_variables=10)
end
