export get_function, get_derivative, get_domain_lb, get_domain_ub, get_domain, get_partition

"Abstract formulation class"
abstract type AbstractFormulation end

"Abstract variable index class"
abstract type AbstractVariableIndices end

"ConstraintMatrix is contains the pair (A, b) where A is a sparse matrix"
const ConstraintMatrix = Pair{SparseMatrixCSC{<:Real,Int64},Vector{<:Real}}

"Vertex is a pair (x, y)"
const Vertex = Pair{<:Real,<:Real}

"""
    Struct to hold the function data with the inputs provided by the user.

It takes in the function and its derivative (as lambda expressions), the base partition
that the user provides, the partition (refinement of the base partition), 3 tolerance values
(1) error tolerance denotes the strength of the relaxation (the closer to zero the stronger
    the relaxation),
(2) length tolerance (maximum difference between successive derivative values), and
(3) derivative tolerance denotes the tolerance for checking equality of derivative values at
    subsequent partition points, and finally, the maximum number of additional partition intervals.
"""
struct FunctionData
    f::Function
    f_dash::Function
    base_partition::Vector{<:Real}
    partition::Vector{<:Real}  # refinement of `base_partition` based on specified limits
    error_tolerance::Float64  # maximum allowed error bound
    length_tolerance::Float64  # maximum allowed distance between successive partition points
    derivative_tolerance::Float64  # maximum difference between successive derivative values
    num_additional_binary_variables::Int64  # maximum number of additional partition intervals
end

"Getters for FunctionData struct"
@inline get_function(uf::FunctionData) = uf.f
@inline get_derivative(uf::FunctionData) = uf.f_dash
@inline get_domain_lb(uf::FunctionData) = uf.partition[1]
@inline get_domain_ub(uf::FunctionData) = uf.partition[end]
@inline get_domain(uf::FunctionData) = get_domain_lb(uf), get_domain_ub(uf)
@inline get_partition(uf::FunctionData) = uf.partition

"""
    Struct to hold the constraint data.
    It contains info to create the sparse constraint matrices.
"""
mutable struct ConstraintData
    constraint_row_indices::Vector{Int64}
    constraint_column_indices::Vector{Int64}
    constraint_coefficients::Vector{Float64}
    rhs_row_indices::Vector{Int64}
    rhs_values::Vector{Float64}
    num_constraints::Int64
end

"Empty Constructor for ConstraintData"
function ConstraintData()::ConstraintData
    row_indices = Int64[]
    col_indices = Int64[]
    coefs = Float64[]
    rhs_row_indices = Int64[]
    rhs_values = Float64[]
    num_constraints = 0
    return ConstraintData(
        row_indices,
        col_indices,
        coefs,
        rhs_row_indices,
        rhs_values,
        num_constraints,
    )
end

"Helper function to construct ConstraintMatrix from ConstraintData and number of variables"
function get_constraint_matrix(
    constraint_data::ConstraintData,
    num_variables::Int64,
)::ConstraintMatrix
    A = sparse(
        constraint_data.constraint_row_indices,
        constraint_data.constraint_column_indices,
        constraint_data.constraint_coefficients,
        constraint_data.num_constraints,
        num_variables,
    )
    b = sparsevec(
        constraint_data.rhs_row_indices,
        constraint_data.rhs_values,
        constraint_data.num_constraints,
    )
    return Pair(A, Vector(b))
end

"""
    get_tangent_vertex(prev_secant_vertex, next_secant_vertex, derivative)

Return (x,y) coordinates of the intersection of tangents drawn at `prev_secant_vertex` and
`next_secant_vertex`.
"""
function get_tangent_vertex(
    prev_secant_vertex::Vertex,
    next_secant_vertex::Vertex,
    derivative::Function,
)::Vertex
    x_prev, f_prev = prev_secant_vertex
    x_next, f_next = next_secant_vertex
    d_prev, d_next = derivative(x_prev), derivative(x_next)
    x_t = (f_next - f_prev + (d_prev * x_prev) - (d_next * x_next)) / (d_prev - d_next)
    y_t = f_prev + (d_prev * (x_t - x_prev))
    return Pair(x_t, y_t)
end

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
            tv = get_tangent_vertex(
                secant_vertices[end-1],
                secant_vertices[end],
                function_data.f_dash,
            )
            push!(tangent_vertices, tv)
        end
    end
    return Pair(secant_vertices, tangent_vertices)
end

"""
    add_coeff!(constraint_data, row, col, value)

Add the coefficient `value` of the variable with index `col` to the constraint with index `row` to
`constraint_data`.
"""
function add_coeff!(constraint_data::ConstraintData, row::Int64, col::Int64, value::Float64)
    push!(constraint_data.constraint_row_indices, row)
    push!(constraint_data.constraint_column_indices, col)
    push!(constraint_data.constraint_coefficients, value)
end

"""
    add_rhs!(constraint_data, row, value)

Add the right-hand-side `value` for row `row` to `constraint_data`.
"""
function add_rhs!(constraint_data::ConstraintData, row::Int64, value::Float64)
    push!(constraint_data.rhs_row_indices, row)
    push!(constraint_data.rhs_values, value)
end

"ConstraintData validator - checks for correctness"
function validate(constraint_data::ConstraintData)
    num_entries = length(constraint_data.constraint_row_indices)
    @assert num_entries == length(constraint_data.constraint_column_indices)
    @assert num_entries == length(constraint_data.constraint_coefficients)
    @assert length(constraint_data.rhs_row_indices) == length(constraint_data.rhs_values)
end

"Input data point validator"
function validate_point(function_data::FunctionData, x::Float64)
    if !isfinite(x) || abs(x) >= ∞
        Memento.error(_LOGGER, "all partition points must be finite")
    end

    fx = function_data.f(x)
    if abs(fx) >= ∞
        Memento.error(_LOGGER, "absolute function value at $x larger than $∞")
    end

    dx = function_data.f_dash(x)
    if abs(dx) >= ∞
        Memento.error(_LOGGER, "absolute derivative value at $x larger than $∞")
    end
end

"Input data validator"
function validate(function_data::FunctionData)
    if length(function_data.partition) < 2
        Memento.error(_LOGGER, "partition must have at least 2 points")
    end

    x_prev::Float64 = NaN64
    d_prev::Float64 = NaN64
    for x in function_data.partition
        validate_point(function_data, x)
        dx = function_data.f_dash(x)
        if !isnan(x_prev)
            if x <= x_prev
                Memento.error(
                    _LOGGER,
                    "partition must be sorted, violation for $x, $x_prev",
                )
            end
            if x - x_prev <= function_data.length_tolerance
                Memento.error(
                    _LOGGER,
                    string(
                        "$x_prev and $x difference less than ",
                        "$(function_data.length_tolerance)",
                    ),
                )
            end
            if abs(dx - d_prev) <= function_data.derivative_tolerance
                Memento.error(
                    _LOGGER,
                    string(
                        "difference of derivatives at $x and $x_prev less than ",
                        "$(function_data.derivative_tolerance)",
                    ),
                )
            end
        end
        x_prev = x
        d_prev = dx
    end

    Memento.debug(_LOGGER, "input data valid.")
end

"Partition refinement schemes (interval bisection)"
function refine_partition!(function_data::FunctionData)
    # Don't refine the partition if no additional constraints are specified.
    if isnan(function_data.error_tolerance) &&
       function_data.num_additional_binary_variables <= 0
        return
    end

    error_queue = get_error_queue(function_data)
    num_partitions_added = 0
    partition = function_data.partition
    while true
        if !is_refinement_feasible(function_data, error_queue)
            Memento.debug(_LOGGER, "stopping refinement")
            return
        end

        start, max_error = peek(error_queue)
        x_start = partition[start]
        x_end = partition[start+1]
        Memento.debug(_LOGGER, "max error: $max_error between $x_start, $x_end")

        # Errors of partition intervals in `error_queue` are indexed by positions of interval
        # starts in `refined_partition`. As we will be inserting `x_new` into `refined_partition`
        # between positions `start` and `start+1`, the positions of interval-starts after `x_new`
        # will all increase by 1 after the insertions. Upade the queue with this new indexing.
        num_starts = length(partition)
        for i = num_starts:-1:start+1
            error_queue[i] = error_queue[i-1]
        end

        # Add errors of the new partition intervals to the queue.
        x_new = 0.5 * (x_start + x_end)
        error_queue[start] = get_error_bound(function_data.f_dash, x_start, x_new)
        error_queue[start+1] = get_error_bound(function_data.f_dash, x_new, x_end)

        # Add `x_new` to `refined_partition`.
        splice!(partition, start+1:start, x_new)
        num_partitions_added += 1
    end
end

"This function checks if the refinement is feasible"
function is_refinement_feasible(
    function_data::FunctionData,
    error_queue::PriorityQueue,
)::Bool
    # Check if error bound is below tolerance.
    start, max_error = peek(error_queue)
    if !isnan(function_data.error_tolerance) && max_error <= function_data.error_tolerance
        Memento.debug(
            _LOGGER,
            "error: $max_error less than limit: $(function_data.error_tolerance)",
        )
        return false
    end

    # Check if the maximum allowed number of binary variables have been added.
    if function_data.num_additional_binary_variables > 0
        num_added = length(function_data.partition) - length(function_data.base_partition)
        if num_added >= function_data.num_additional_binary_variables
            Memento.debug(_LOGGER, "number of new binary variables: $num_added")
            Memento.debug(
                _LOGGER,
                "budget: $(function_data.num_additional_binary_variables)",
            )
            return false
        end
    end

    # Check if new partition lengths are smaller than allowed.
    x_start, x_end = function_data.partition[start], function_data.partition[start+1]
    if x_end - x_start <= (2 * function_data.length_tolerance)
        Memento.debug(_LOGGER, "start: $x_start, end: $x_end")
        Memento.debug(_LOGGER, "new interval length will be too small")
        return false
    end

    # Check if function and derivative are well-defined at the new point.
    x_new = 0.5 * (x_start + x_end)
    validate_point(function_data, x_new)

    # Check if derivative differences at endpoints of new partitions are smaller than allowed.
    d_start, d_end = function_data.f_dash(x_start), function_data.f_dash(x_end)
    d_new = function_data.f_dash(x_new)
    if (
        abs(d_new - d_start) <= function_data.derivative_tolerance ||
        abs(d_end - d_new) <= function_data.derivative_tolerance
    )
        Memento.debug(_LOGGER, "d_start: $d_start, d_new: $d_new, d_end: $d_end")
        Memento.debug(_LOGGER, "adjacent derivative difference will be too small")
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
    for i = 1:length(function_data.partition)-1
        lb = function_data.partition[i]
        ub = function_data.partition[i+1]
        pq[i] = get_error_bound(function_data.f_dash, lb, ub)
    end
    return pq
end

"""
    get_error_bound(derivative, lb, ub)

Get error bound of a function with derivative `derivative` in the closed interval `[lb,ub]`.
"""
function get_error_bound(derivative::Function, lb::Float64, ub::Float64)
    return (ub - lb) * abs(derivative(ub) - derivative(lb)) / 4.0
end
