"""
    _check_consistency(formulation_info, num_vars) 
Checks consistency between provided variables and partition sizes 
"""
function _check_consistency(formulation_info::FormulationInfo, num_vars)::Bool
    var = formulation_info.variables
    if haskey(var, :z)
        (length(var[:z]) == num_vars) && (return true)
    end
    return false
end

"""
    _build_univariate_milp_relaxation!(m,x,y,function_data,pre_base_name)

Return a MILPRelaxation object with constraint and RHS information of the MILP
formulation of the polyhedral relaxation.
"""
function _build_univariate_milp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    univariate_function_data::UnivariateFunctionData,
    variable_pre_base_name::AbstractString,
    reuse::FormulationInfo,
)::FormulationInfo
    sec_vs, tan_vs = _collect_vertices(univariate_function_data)
    formulation_info = FormulationInfo()

    # create variables
    num_vars = length(univariate_function_data.partition) - 1

    reuse_variables = reuse.variables

    is_consistent = _check_consistency(reuse, num_vars)

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
    z =
        (is_consistent) ? reuse_variables[:z] :
        JuMP.@variable(
            m,
            [1:num_vars],
            binary = true,
            base_name = variable_pre_base_name * "_z"
        )
    formulation_info.variables[:z] = z

    # add x constraints
    JuMP.@constraint(
        m,
        x ==
        sec_vs[1][1] + sum(
            delta_1[i] * (tan_vs[i][1] - sec_vs[i][1]) +
            delta_2[i] * (sec_vs[i+1][1] - sec_vs[i][1]) for i in 1:num_vars
        )
    )

    # add y constraints
    JuMP.@constraint(
        m,
        y ==
        sec_vs[1][2] + sum(
            delta_1[i] * (tan_vs[i][2] - sec_vs[i][2]) +
            delta_2[i] * (sec_vs[i+1][2] - sec_vs[i][2]) for i in 1:num_vars
        )
    )

    # add first delta constraint
    JuMP.@constraint(m, delta_1[1] + delta_2[1] <= 1)
    # add linking constraints between delta_1, delta_2 and z
    JuMP.@constraint(m, [i = 2:num_vars], delta_1[i] + delta_2[i] <= z[i-1])
    JuMP.@constraint(m, [i = 2:num_vars], z[i-1] <= delta_2[i-1])

    return formulation_info
end
