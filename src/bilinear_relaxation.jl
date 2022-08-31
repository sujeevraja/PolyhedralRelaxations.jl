"""
    _build_mccormick_relaxation!(m, x, y, z)

McCormick relaxation of binlinear term 
```
z >= JuMP.lower_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.lower_bound(x)*JuMP.lower_bound(y)
z >= JuMP.upper_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.upper_bound(x)*JuMP.upper_bound(y)
z <= JuMP.lower_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.lower_bound(x)*JuMP.upper_bound(y)
z <= JuMP.upper_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.upper_bound(x)*JuMP.lower_bound(y)
```
"""
function _build_mccormick_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    z::JuMP.VariableRef,
)::FormulationInfo
    x_lb, x_ub = _variable_domain(x)
    y_lb, y_ub = _variable_domain(y)
    formulation_info = FormulationInfo()
    JuMP.@constraint(m, z >= x_lb * y + y_lb * x - x_lb * y_lb)
    JuMP.@constraint(m, z >= x_ub * y + y_ub * x - x_ub * y_ub)
    JuMP.@constraint(m, z <= x_lb * y + y_ub * x - x_lb * y_ub)
    JuMP.@constraint(m, z <= x_ub * y + y_lb * x - x_ub * y_lb)
    return formulation_info
end

"""
    _check_partition_variable_consistency(formulation_info, num_vars) 
Checks consistency between provided variables and partition sizes 
"""
function _check_partition_variable_consistency(
    formulation_info::FormulationInfo,
    num_vars,
)::Bool
    var = formulation_info.variables
    if haskey(var, :bin)
        (length(var[:bin]) == num_vars) && (return true)
    end
    return false
end

"""
    _build_bilinear_relaxation!(m, x, y, z, x_partition, y_partition, pre_base_name)

Build incremental formulation for ``z = xy`` given partition data.
"""
function _build_bilinear_milp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    z::JuMP.VariableRef,
    x_partition::Vector{<:Real},
    y_partition::Vector{<:Real},
    variable_pre_base_name::AbstractString,
    reuse::FormulationInfo,
)::FormulationInfo
    origin_vs, non_origin_vs =
        _collect_bilinear_vertices(x_partition, y_partition)
    formulation_info = FormulationInfo()
    num_vars = max(length(x_partition), length(y_partition)) - 1

    reuse_variables = reuse.variables

    is_consistent = _check_partition_variable_consistency(reuse, num_vars)

    # add variables
    delta_1 =
        formulation_info.variables[:delta_1] = JuMP.@variable(
            m,
            [1:num_vars],
            lower_bound = 0.0,
            upper_bound = 1.0,
            base_name = variable_pre_base_name * "delta_1"
        )
    delta_2 =
        formulation_info.variables[:delta_2] = JuMP.@variable(
            m,
            [1:num_vars],
            lower_bound = 0.0,
            upper_bound = 1.0,
            base_name = variable_pre_base_name * "delta_2"
        )
    delta_3 =
        formulation_info.variables[:delta_3] = JuMP.@variable(
            m,
            [1:num_vars],
            lower_bound = 0.0,
            upper_bound = 1.0,
            base_name = variable_pre_base_name * "delta_3"
        )
    z_bin =
        (is_consistent) ? reuse_variables[:bin] :
        JuMP.@variable(
            m,
            [1:num_vars],
            binary = true,
            base_name = variable_pre_base_name * "z"
        )
    formulation_info.variables[:bin] = z_bin

    # add x constraints
    JuMP.@constraint(
        m,
        x ==
        origin_vs[1][1] + sum(
            delta_1[i] * (non_origin_vs[i][1] - origin_vs[i][1]) +
            delta_2[i] * (non_origin_vs[i+1][1] - origin_vs[i][1]) +
            delta_3[i] * (origin_vs[i+1][1] - origin_vs[i][1]) for
            i in 1:num_vars
        )
    )

    # add y constraints
    JuMP.@constraint(
        m,
        y ==
        origin_vs[1][2] + sum(
            delta_1[i] * (non_origin_vs[i][2] - origin_vs[i][2]) +
            delta_2[i] * (non_origin_vs[i+1][2] - origin_vs[i][2]) +
            delta_3[i] * (origin_vs[i+1][2] - origin_vs[i][2]) for
            i in 1:num_vars
        )
    )

    # add z constraints
    JuMP.@constraint(
        m,
        z ==
        origin_vs[1][3] + sum(
            delta_1[i] * (non_origin_vs[i][3] - origin_vs[i][3]) +
            delta_2[i] * (non_origin_vs[i+1][3] - origin_vs[i][3]) +
            delta_3[i] * (origin_vs[i+1][3] - origin_vs[i][3]) for
            i in 1:num_vars
        )
    )

    # add first delta constraint
    JuMP.@constraint(m, delta_1[1] + delta_2[1] + delta_3[1] <= 1)

    # add linking constraints between delta_1, delta_2 and z
    JuMP.@constraint(
        m,
        [i = 2:num_vars],
        delta_1[i] + delta_2[i] + delta_3[i] <= z_bin[i-1]
    )
    JuMP.@constraint(m, [i = 2:num_vars], z_bin[i-1] <= delta_3[i-1])

    return formulation_info
end
