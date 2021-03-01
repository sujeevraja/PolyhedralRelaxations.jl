"Vertex is a pair ``(x, y)``"
const Vertex = Pair{<:Real,<:Real}

"""
The struct `FormulationInfo` holds two dictionaries, one for variable
references and the other for constraint references.
"""
struct FormulationInfo
    variables::Dict{Symbol,Any}
    constraints::Dict{Symbol,Any}
end

"Empty contructor for struct `FormulationInfo`"
FormulationInfo()::FormulationInfo = FormulationInfo(Dict{Symbol,Any}(), Dict{Symbol,Any}())

"""
The struct `UnivariateFunctionData` holds the inputs provided by the user. It
takes in the function and its derivative (as lambda expressions), the partition
on the independent variable that the user provides, and the following 3
tolerance values:

* error tolerance denotes the strength of the relaxation (the closer to zero
    the stronger the relaxation)
* length tolerance (maximum difference between successive derivative values)
* derivative tolerance denotes the tolerance for checking equality of
    derivative values at subsequent partition points, and finally, the maximum
    number of additional partition intervals.
"""
struct UnivariateFunctionData
    f::Function
    f_dash::Function
    partition::Vector{<:Real}
    error_tolerance::Float64
    length_tolerance::Float64
    derivative_tolerance::Float64
    num_additional_partitions::Int64
end

"""
    get_tangent_vertex(prev_secant_vertex, next_secant_vertex, derivative)

Return ``(x,y)`` coordinates of the intersection of tangents drawn at
`prev_secant_vertex` and `next_secant_vertex`.
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
    collect_vertices(univariate_function_data)

Return a pair of lists with secant vertices as the first element and tangent
vertices as the second element.

Each element in the secant vertex list is a pair (x,y) where
* x is an element of `univariate_function_data.partition`,
* y is the value of the given univariate function at x.
All secant vertices lie on the curve.

Each element in the tangent vertex list is also a pair (x,y). Each position i
of the list contains the vertex formed by intersection of tangents of the curve
`y=univariate_function_data.f(x)` at `secant_vertices[i]` and
`secant_vertices[i+1]`. No tangent vertex will lie on the curve (except for the
trivial case where the curve is linear, and all triangles are flat lines).
"""
function collect_vertices(
    univariate_function_data::UnivariateFunctionData,
)::Pair{Vector{Vertex},Vector{Vertex}}
    secant_vertices, tangent_vertices = Vertex[], Vertex[]
    for x in univariate_function_data.partition
        push!(secant_vertices, Pair(x, univariate_function_data.f(x)))
        if length(secant_vertices) >= 2
            tv = get_tangent_vertex(
                secant_vertices[end - 1],
                secant_vertices[end],
                univariate_function_data.f_dash,
            )
            push!(tangent_vertices, tv)
        end
    end
    return Pair(secant_vertices, tangent_vertices)
end

"""
    validate(x, partition)

Variable bounds and partition consistency checker
"""
function validate(x::JuMP.VariableRef, partition::Vector{<:Real})
    x_lb = isinf(JuMP.lower_bound(x)) ? -Inf : JuMP.lower_bound(x)
    x_ub = isinf(JuMP.upper_bound(x)) ? Inf : JuMP.upper_bound(x)
    lb = partition[1]
    ub = partition[end]
    if isfinite(x_lb) && !isapprox(x_lb, lb, rtol=1e-8)
        Memento.error(_LOGGER, "partition lower bound and variable lower bound not equal")
    end 
    if isfinite(x_ub) && !isapprox(x_ub, ub, rtol=1e-8)
        Memento.error(_LOGGER, "partition upper bound and variable upper bound not equal")
    end 
    (isinf(x_lb)) && (JuMP.set_lower_bound(x, lb))
    (isinf(x_ub)) && (JuMP.set_upper_bound(x, ub))
end

"""
    validate_point(univariate_function_data, x)

Input data point validator
"""
function validate_point(univariate_function_data::UnivariateFunctionData, x::Float64)
    if !isfinite(x) || abs(x) >= INF 
        Memento.error(_LOGGER, "all partition points must be finite")
    end

    fx = univariate_function_data.f(x)
    if abs(fx) >= INF
        Memento.error(_LOGGER, "absolute function value at $x larger than $∞")
    end

    dx = univariate_function_data.f_dash(x)
    if abs(dx) >= INF
        Memento.error(_LOGGER, "absolute derivative value at $x larger than $∞")
    end
end

"""
    validate(univariate_function_data)

Input univariate function data validator
"""
function validate(univariate_function_data::UnivariateFunctionData)
    if length(univariate_function_data.partition) < 2
        Memento.error(_LOGGER, "partition must have at least 2 points")
    end

    x_prev::Float64 = NaN64
    d_prev::Float64 = NaN64
    for x in univariate_function_data.partition
        validate_point(univariate_function_data, x)
        dx = univariate_function_data.f_dash(x)
        if !isnan(x_prev)
            if x <= x_prev
                Memento.error(
                    _LOGGER,
                    "partition must be sorted, violation for $x, $x_prev",
                )
            end
            if x - x_prev <= univariate_function_data.length_tolerance
                Memento.error(
                    _LOGGER,
                    string(
                        "$x_prev and $x difference less than ",
                        "$(univariate_function_data.length_tolerance)",
                    ),
                )
            end
            if abs(dx - d_prev) <= univariate_function_data.derivative_tolerance
                Memento.error(
                    _LOGGER,
                    string(
                        "difference of derivatives at $x and $x_prev less than ",
                        "$(univariate_function_data.derivative_tolerance)",
                    ),
                )
            end
        end
        x_prev = x
        d_prev = dx
    end

    Memento.debug(_LOGGER, "input data valid.")
end

"""
    refine_partition!(univariate_function_data)

Partition refinement schemes (interval bisection)
"""
function refine_partition!(univariate_function_data::UnivariateFunctionData)
    # Don't refine the partition if no additional constraints are specified.
    if isnan(univariate_function_data.error_tolerance) &&
       univariate_function_data.num_additional_partitions <= 0
        return
    end

    error_queue = get_error_queue(univariate_function_data)
    num_partitions_added = 0
    partition = univariate_function_data.partition
    while true
        if !is_refinement_feasible(univariate_function_data, error_queue)
            Memento.debug(_LOGGER, "stopping refinement")
            return
        end

        start, max_error = peek(error_queue)
        x_start = partition[start]
        x_end = partition[start + 1]
        Memento.debug(_LOGGER, "max error: $max_error between $x_start, $x_end")

        """
        Errors of partition intervals in `error_queue` are indexed by positions
        of interval starts in `refined_partition`. As we will be inserting
        `x_new` into `refined_partition` between positions `start` and
        `start+1`, the positions of interval-starts after `x_new` will all
        increase by 1 after the insertions. Upade the queue with this new
        indexing.
        """
        num_starts = length(partition)
        for i = num_starts:-1:start + 1
            error_queue[i] = error_queue[i - 1]
        end

        # Add errors of the new partition intervals to the queue.
        x_new = 0.5 * (x_start + x_end)
        error_queue[start] =
            get_error_bound(univariate_function_data.f_dash, x_start, x_new)
        error_queue[start + 1] =
            get_error_bound(univariate_function_data.f_dash, x_new, x_end)

        # Add `x_new` to `refined_partition`.
        splice!(partition, start + 1:start, x_new)
        num_partitions_added += 1
    end
end

"""
    is_refinement_feasible(univariate_function_data, error_queue)

This function checks if the refinement is feasible
"""
function is_refinement_feasible(
    univariate_function_data::UnivariateFunctionData,
    error_queue::PriorityQueue,
)::Bool
    # Check if error bound is below tolerance.
    start, max_error = peek(error_queue)
    if !isnan(univariate_function_data.error_tolerance) &&
       max_error <= univariate_function_data.error_tolerance
        Memento.debug(
            _LOGGER,
            "error: $max_error less than limit: $(univariate_function_data.error_tolerance)",
        )
        return false
    end

    # Check if the maximum allowed number of binary variables have been added.
    if univariate_function_data.num_additional_partitions > 0
        num_added =
            length(univariate_function_data.partition) -
            length(univariate_function_data.base_partition)
        if num_added >= univariate_function_data.num_additional_partitions
            Memento.debug(_LOGGER, "number of new binary variables: $num_added")
            Memento.debug(
                _LOGGER,
                "budget: $(univariate_function_data.num_additional_partitions)",
            )
            return false
        end
    end

    # Check if new partition lengths are smaller than allowed.
    x_start, x_end = univariate_function_data.partition[start],
    univariate_function_data.partition[start + 1]
    if x_end - x_start <= (2 * univariate_function_data.length_tolerance)
        Memento.debug(_LOGGER, "start: $x_start, end: $x_end")
        Memento.debug(_LOGGER, "new interval length will be too small")
        return false
    end

    # Check if function and derivative are well-defined at the new point.
    x_new = 0.5 * (x_start + x_end)
    validate_point(univariate_function_data, x_new)

    # Check if derivative differences at endpoints of new partitions are smaller than allowed.
    d_start, d_end =
        univariate_function_data.f_dash(x_start), univariate_function_data.f_dash(x_end)
    d_new = univariate_function_data.f_dash(x_new)
    if (
        abs(d_new - d_start) <= univariate_function_data.derivative_tolerance ||
        abs(d_end - d_new) <= univariate_function_data.derivative_tolerance
    )
        Memento.debug(_LOGGER, "d_start: $d_start, d_new: $d_new, d_end: $d_end")
        Memento.debug(_LOGGER, "adjacent derivative difference will be too small")
        return false
    end

    return true
end

"""
    get_error_queue(univariate_function_data)

Build a max-priority-queue holding errors of each partition interval in
`univariate_function_data.partition`. Keys of the queue are indices of starting
positions of the partition interval in `univariate_function_data.partition`.
Priorities are error bounds of partition intervals. The queue is built as
a max-queue for easy access to the  maximum error.
"""
function get_error_queue(univariate_function_data::UnivariateFunctionData)::PriorityQueue
    pq = PriorityQueue{Int64,Float64}(Base.Order.Reverse)
    for i = 1:length(univariate_function_data.partition) - 1
        lb = univariate_function_data.partition[i]
        ub = univariate_function_data.partition[i + 1]
        pq[i] = get_error_bound(univariate_function_data.f_dash, lb, ub)
    end
    return pq
end

"""
    get_error_bound(derivative, lb, ub)

Get error bound of a function with derivative `derivative` in the closed
interval `[lb,ub]`.
"""
function get_error_bound(derivative::Function, lb::Float64, ub::Float64)
    return (ub - lb) * abs(derivative(ub) - derivative(lb)) / 4.0
end

"""
    get_max_error_bound(univariate_function_data)

Compute and return maximum value of error bound among all partition intervals.
"""
function get_max_error_bound(univariate_function_data::UnivariateFunctionData)::Float64
    max_err = -Inf
    for i in length(univariate_function_data.partition) - 1
        err = get_error_bound(
            univariate_function_data.f_dash,
            univariate_function_data.partition[i],
            univariate_function_data.partition[i + 1],
        )
        max_err = max(err, max_err)
    end
    return max_err
end
