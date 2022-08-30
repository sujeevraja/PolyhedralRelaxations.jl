export RefinementInfo

struct RefinementInfo
    num_additional_points::Int64
    refined_largest_partition::Bool
end

RefinementInfo() = RefinementInfo(0, false)

"""
    _bisect_all!(partition::Vector{<:Real}, refinement_width_tol::Float64)::RefinementInfo

This scheme bisects every partition provided the width of the partition is greater
than `refinement_width_tol` value.
"""
function _bisect_all!(
    partition::Vector{<:Real},
    refinement_width_tol::Float64,
)::RefinementInfo
    additional_points = []
    for i in 1:(length(partition)-1)
        width = partition[i+1] - partition[i]
        (width < refinement_width_tol) && (continue)
        push!(additional_points, partition[i] + 0.5 * width)
    end
    sort!(append!(partition, additional_points))
    return RefinementInfo(length(additional_points), false)
end

"""
    _at_point!(partition::Vector{<:Real}, point::T where {T<:Real})::RefinementInfo 

This scheme adds the point to the partition. 
"""
function _add_point!(
    partition::Vector{<:Real},
    point::T where {T<:Real},
)::RefinementInfo
    sort!(append!(partition, [point]))
    return RefinementInfo(1, false)
end

""" 
    _bisect!(partition::Vector{<:Real}, lower::Float64, upper::Float64)::RefinementInfo

This scheme bisects the partition defined by `lower` and `upper` values.
"""
function _bisect!(
    partition::Vector{<:Real},
    lower::Float64,
    upper::Float64,
)::RefinementInfo
    halfway = (lower + upper) / 2.0
    sort!(append!(partition, [halfway]))
    return RefinementInfo(1, false)
end

"""
    _non_uniform!(partition, point, lower, upper, added_point_width_tol, ratio)::RefinementInfo

This scheme adds a small partition around the `point` using the `ratio` provided while 
satisfying the `added_point_width_tolerance`
"""
function _non_uniform!(
    partition::Vector{<:Real},
    point::T where {T<:Real},
    lower::Float64,
    upper::Float64,
    added_point_width_tol::Float64,
    ratio::Float64,
)
    # potential points .
    left = point - ratio * (point - lower)
    right = point + ratio * (upper - point)
    if abs(left - lower) > added_point_width_tol
        # Left point can be added
        if abs(upper - right) > added_point_width_tol
            # Both left and right points can be added
            sort!(append!(partition, [left, right]))
            return RefinementInfo(2, false)
        else
            # Only the left point can be added
            sort!(append!(partition, [left]))
            return RefinementInfo(1, false)
        end
    elseif abs(upper - right) > added_point_width_tol
        # Only the right point can be added
        sort!(append!(partition, [left]))
        return RefinementInfo(1, false)
    else
        error("Issue encountered when trying to add points in a subinterval")
    end
    return RefinementInfo
end

"""
    _refine_partition!(
        partition::Vector{<:Real}, 
        point::T where {T<:Real},
        refinement_type::REFINEMENT_TYPE,
        refinement_ratio::Float64, 
        refinement_width_tol::Float64,
        refinement_added_point_width_tolerance::Float64,
        refine_largest::Bool)::RefinementInfo

Internal helper function to help in partition refinement         
"""
function _refine_partition!(
    partition::Vector{<:Real},
    point::T where {T<:Real},
    refinement_type::Symbol,
    refinement_ratio::Float64,
    refinement_width_tol::Float64,
    refinement_added_point_width_tolerance::Float64,
    refine_largest::Bool,
)::RefinementInfo
    if (refinement_type == :bisect_all)
        return _bisect_all!(partition, refinement_width_tol)
    end

    # Flag indicating if the largest subinterval needs to be refined. This will
    # only happen when the subinterval containing the point is too small for futher
    # refinement
    adapt = false
    # Flag indicating whether the subinterval containing the point has been found
    found = false
    # Tracking the largest subinterval found in the partition 
    width = 0
    largest = [-Inf, Inf]

    # Iterating over all subintervals in the partition
    for index in 1:length(partition)-1
        # Grabbing the lower and upper value of the subinterval
        lower = partition[index]
        upper = partition[index+1]
        # Updating the largest subinterval found
        if upper - lower > width
            # New largest subinterval found
            width = upper - lower
            largest = [lower, upper]
        end

        # Checking if the point is effectively equal to one of the endpoints
        # of the subinterval or the point is too close to add a point. 
        # If so, do not further refine the partition
        (index == 1 && abs(lower - point) <= EPS_ZERO) &&
            (return RefinementInfo())
        (abs(upper - point) <= EPS_ZERO) && (return RefinementInfo())

        # Checking if the point is in the subinterval [lower, upper] if the subinterval
        # containing the point has not yet been found
        if !found && lower < point && point < upper
            # Subinterval containing the point has been found, so we no longer
            # need to continue searching for the subinterval containing the point. 
            # If the subinterval containing the point is too small, we still 
            # will need to find the largest subinterval in the partition. 
            found = true
            # Checking if the subinterval containing the point is too small for further
            # refinement. 
            if abs(upper - lower) < refinement_width_tol
                # Subinterval containing the point is too small, so refine the
                # largest subinterval in the partition instead.
                adapt = true
            else
                # Subinterval is sufficiently large for further refinement. 
                if refinement_type == :at_point
                    return _at_point!(partition, point)
                elseif refinement_type == :bisect
                    return _bisect!(partition, lower, upper)
                else
                    return _non_uniform!(
                        partition,
                        point,
                        lower,
                        upper,
                        refinement_added_point_width_tolerance,
                        refinement_ratio,
                    )
                end
            end
        end
    end

    # The point may have been found but the subinterval containing it may have been too
    # small for further refinement. If this is the case, refine the largest subinterval
    # using bisection instead.
    if adapt && refine_largest
        # Checking the largest subinterval is sufficiently large for further refinement
        if width >= refinement_width_tol
            # Largest subinterval is sufficiently large for further refinement.
            halfway = sum(largest) / 2.0
            sort!(append!(partition, [halfway]))
            return (RefinementInfo(1, false))
        end
    end
end
