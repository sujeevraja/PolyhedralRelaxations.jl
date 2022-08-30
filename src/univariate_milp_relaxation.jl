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
    pre_base_name::AbstractString,
)::FormulationInfo
    sec_vs, tan_vs = _collect_vertices(univariate_function_data)
    formulation_info = FormulationInfo()

    # create variables
    num_vars = length(univariate_function_data.partition) - 1

    delta_1 =
        formulation_info.variables[:delta_1] = JuMP.@variable(
            m,
            [1:num_vars],
            lower_bound = 0.0,
            upper_bound = 1.0,
            base_name = pre_base_name * "delta_1"
        )
    delta_2 =
        formulation_info.variables[:delta_2] = JuMP.@variable(
            m,
            [1:num_vars],
            lower_bound = 0.0,
            upper_bound = 1.0,
            base_name = pre_base_name * "delta_2"
        )
    z =
        formulation_info.variables[:z] = JuMP.@variable(
            m,
            [1:num_vars],
            binary = true,
            base_name = pre_base_name * "z"
        )

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
