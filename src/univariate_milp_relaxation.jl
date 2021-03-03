"""
    _build_univariate_milp_relaxation!(m,x,y,function_data)

Return a MILPRelaxation object with constraint and RHS information of the MILP
formulation of the polyhedral relaxation.
"""
function _build_univariate_milp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    univariate_function_data::UnivariateFunctionData,
)::FormulationInfo
    sec_vs, tan_vs = _collect_vertices(univariate_function_data)
    formulation_info = FormulationInfo()

    # create variables
    num_vars = length(univariate_function_data.partition) - 1
    @variable(m, 0 <= delta_1[1:num_vars] <= 1)
    @variable(m, 0 <= delta_2[1:num_vars] <= 1)
    @variable(m, 0 <= z[1:num_vars] <= 1, Bin)
    formulation_info.variables[:delta_1] = delta_1
    formulation_info.variables[:delta_2] = delta_2
    formulation_info.variables[:z] = z

    # add x constraints
    formulation_info.constraints[:x] = @constraint(
        m,
        x ==
        sec_vs[1][1] + sum(
            delta_1[i] * (tan_vs[i][1] - sec_vs[i][1]) +
            delta_2[i] * (sec_vs[i+1][1] - sec_vs[i][1]) for i = 1:num_vars
        )
    )

    # add y constraints
    formulation_info.constraints[:y] = @constraint(
        m,
        y ==
        sec_vs[1][2] + sum(
            delta_1[i] * (tan_vs[i][2] - sec_vs[i][2]) +
            delta_2[i] * (sec_vs[i+1][2] - sec_vs[i][2]) for i = 1:num_vars
        )
    )

    # add first delta constraint
    formulation_info.constraints[:first_delta] =
        @constraint(m, delta_1[1] + delta_2[1] <= 1)

    # add linking constraints between delta_1, delta_2 and z
    @constraint(m, zcon1[i = 2:num_vars], delta_1[i] + delta_2[i] <= z[i-1])
    formulation_info.constraints[:below_z] = zcon1
    @constraint(m, z_con2[i = 2:num_vars], z[i-1] <= delta_2[i-1])
    formulation_info.constraints[:above_z] = z_con2

    return formulation_info
end
