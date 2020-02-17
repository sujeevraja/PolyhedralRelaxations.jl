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

    # Check if function and derivative are well-defined at the new point.
    x_new = 0.5 * (x_start + x_end)
    validate_point(function_data, x_new)

    # Check if derivative differences at endpoints of new partitions are smaller than allowed.
    d_start, d_end = function_data.d(x_start), function_data.d(x_end)
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
