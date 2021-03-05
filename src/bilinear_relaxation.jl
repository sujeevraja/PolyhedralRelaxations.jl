"""
    function _variable_domain(var)

Computes the valid domain of a given JuMP variable taking into account bounds
and the varaible's implicit bounds (e.g. binary).
"""
function _variable_domain(var::JuMP.VariableRef)
    lb = -Inf
    if JuMP.has_lower_bound(var)
        lb = JuMP.lower_bound(var)
    end
    if JuMP.is_binary(var)
        lb = max(lb, 0.0)
    end

    ub = Inf
    if JuMP.has_upper_bound(var)
        ub = JuMP.upper_bound(var)
    end
    if JuMP.is_binary(var)
        ub = min(ub, 1.0)
    end

    return (lower_bound=lb, upper_bound=ub)
end

"""
    function _build_mccormick_relaxation!(m, x, y, z)

McCormick relaxation of binlinear term 
```
z >= JuMP.lower_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.lower_bound(x)*JuMP.lower_bound(y)
z >= JuMP.upper_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.upper_bound(x)*JuMP.upper_bound(y)
z <= JuMP.lower_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.lower_bound(x)*JuMP.upper_bound(y)
z <= JuMP.upper_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.upper_bound(x)*JuMP.lower_bound(y)
```
"""
function _build_mccormick_relaxation!(m::JuMP.Model, x::JuMP.VariableRef, y::JuMP.VariableRef, z::JuMP.VariableRef)::FormulationInfo
    x_lb, x_ub = _variable_domain(x)
    y_lb, y_ub = _variable_domain(y)

    formulation_info = FormulationInfo()

    formulation_info.constraints[:lb_1] = @constraint(m, z >= x_lb*y + y_lb*x - x_lb*y_lb)
    formulation_info.constraints[:lb_2] = @constraint(m, z >= x_ub*y + y_ub*x - x_ub*y_ub)
    formulation_info.constraints[:ub_1] = @constraint(m, z <= x_lb*y + y_ub*x - x_lb*y_ub)
    formulation_info.constraints[:ub_2] = @constraint(m, z <= x_ub*y + y_lb*x - x_ub*y_lb)

    return formulation_info
end


"""
    _build_bilinear_relaxation!(m, x, y, z, partition_x, partition_y)

Build incremental formulation for ``z = xy`` given partition data.
"""
function _build_bilinear_milp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    z::JuMP.VariableRef,
    partition_x::Vector{<:Real},
    partition_y::Vector{<:Real}
)::FormulationInfo

    return FormulationInfo()
end 