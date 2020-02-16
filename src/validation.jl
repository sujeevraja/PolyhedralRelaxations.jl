
function validate_point(function_data::FunctionData, x::Real)
    if !isfinite(x) || abs(x) >= ∞
        Memento.error(_LOGGER, "all partition points must be finite")
    end

    fx = function_data.f(x)
    if abs(fx) >= ∞
        Memento.error(_LOGGER, "absolute function value at $x larger than $∞")
    end

    dx = function_data.d(x)
    if abs(dx) >= ∞
        Memento.error(_LOGGER, "absolute derivative value at $x larger than $∞")
    end
end

function validate(function_data::FunctionData)
    if length(function_data.partition) < 2
        Memento.error(_LOGGER, "partition must have at least 2 points")
    end

    x_prev::Float64 = NaN64
    d_prev::Float64 = NaN64
    for x in function_data.partition
        validate_point(function_data, x)
        dx = function_data.d(x)
        if !isnan(x_prev)
            if x <= x_prev
                Memento.error(_LOGGER, "partition must be sorted, violation for $x, $x_prev")
            end
            if x - x_prev <= function_data.length_tolerance
                Memento.error(_LOGGER, string(
                    "$x_prev and $x difference less than ",
                    "$(function_data.length_tolerance)"))
            end
            if abs(dx - d_prev) <= function_data.derivative_tolerance
                Memento.error(_LOGGER, string(
                    "difference of derivatives at $x and $x_prev less than ",
                    "$(function_data.derivative_tolerance)"))
            end
        end
        x_prev = x
        d_prev = dx
    end

    Memento.info(_LOGGER, "input data valid.")
end
