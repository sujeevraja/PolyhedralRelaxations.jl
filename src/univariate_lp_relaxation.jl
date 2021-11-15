"""
    _get_lp_relaxation_vertices(univariate_function_data::UnivariateFunctionData)::Vector{Vertex2d}

Return vertices of the LP relaxation of the given function.
"""
function _get_lp_relaxation_vertices(
    univariate_function_data::UnivariateFunctionData,
)::Vector{Vertex2d}
    sec_vs, tan_vs = _collect_vertices(univariate_function_data)
    vertices = Vertex2d[]
    push!(vertices, sec_vs[1])
    append!(vertices, tan_vs)
    push!(vertices, sec_vs[end])
    return vertices
end

"""
    _build_univariate_lp_relaxation!(m, x, y, univariate_function_data, 
        variable_pre_base_name, constraint_pre_base_name)

Build LP relaxation for ``y=f(x)`` given the univariate function data.
"""
function _build_univariate_lp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    univariate_function_data::UnivariateFunctionData,
    variable_pre_base_name::AbstractString,
    constraint_pre_base_name::AbstractString,
    ::FormulationInfo
)::FormulationInfo
    vertices = _get_lp_relaxation_vertices(univariate_function_data)
    num_vars = length(vertices)
    formulation_info = FormulationInfo()

    # add variables 
    lambda = @variable(m, [1:num_vars], 
            lower_bound = 0.0, upper_bound = 1.0, 
            base_name = variable_pre_base_name * "lambda")
    formulation_info.variables[:lambda] = lambda

    # add constraints 
    formulation_info.constraints[:sum_lambda] = 
        @constraint(m, sum(lambda) == 1,
            base_name = constraint_pre_base_name * "sum_lambda") 

    formulation_info.constraints[:x] = 
        @constraint(m, x == sum(lambda[i] * vertices[i][1] for i = 1:num_vars), 
            base_name = constraint_pre_base_name * "x")

    formulation_info.constraints[:y] =
        @constraint(m, y == sum(lambda[i] * vertices[i][2] for i = 1:num_vars),
            base_name = constraint_pre_base_name * "y")

    return formulation_info
end

"""
    _build_univariate_on_off_relaxation!(m, x, y, z, active_when_z_is_one, 
        variable_pre_base_name, constraint_pre_base_name)

Build on-off relaxation for ``y=f(x)`` given the univariate function data.
"""
function _build_univariate_on_off_relaxation(
    m::JuMP.Model, 
    x::JuMP.VariableRef, 
    y::JuMP.VariableRef, 
    z::JuMP.VariableRef, 
    univariate_function_data::UnivariateFunctionData,
    active_when_z_is_one::Bool,
    variable_pre_base_name::AbstractString,
    constraint_pre_base_name::AbstractString
)::FormulationInfo
    vertices = _get_lp_relaxation_vertices(univariate_function_data)
    num_vars = length(vertices)
    formulation_info = FormulationInfo()

    # add variables 
    lambda = @variable(m, [1:num_vars], 
            lower_bound = 0.0, upper_bound = 1.0, 
            base_name = variable_pre_base_name * "lambda")
    formulation_info.variables[:lambda] = lambda

    # add constraints 
    formulation_info.constraints[:sum_lambda] = 
        if active_when_z_is_one == true
            @constraint(m, sum(lambda) == z,
                base_name = constraint_pre_base_name * "sum_lambda") 
        else 
            @constraint(m, sum(lambda) == 1 - z,
                base_name = constraint_pre_base_name * "sum_lambda") 
        end 

    # `x` constraints are slightly different because `x` has to be within its limit when `y=f(x)` is not active 
    x_lb, x_ub = _variable_domain(x)
    formulation_info.constraints[:x_lb] = 
        if active_when_z_is_one == true
            @constraint(m, x >= 
                sum(lambda[i] * vertices[i][1] for i = 1:num_vars) + (1 - z) * x_lb, 
                base_name = constraint_pre_base_name * "x_lb")
        else 
            @constraint(m, x >= 
                sum(lambda[i] * vertices[i][1] for i = 1:num_vars) + z * x_lb, 
                base_name = constraint_pre_base_name * "x_lb")
        end
    
    formulation_info.constraints[:x_ub] = 
        if active_when_z_is_one == true
            @constraint(m, x <= 
                sum(lambda[i] * vertices[i][1] for i = 1:num_vars) + (1 - z) * x_ub, 
                base_name = constraint_pre_base_name * "x_ub")
        else 
            @constraint(m, x <= 
                sum(lambda[i] * vertices[i][1] for i = 1:num_vars) + z * x_ub, 
                base_name = constraint_pre_base_name * "x_ub")
        end 
    
    formulation_info.constraints[:y] =
        @constraint(m, y == sum(lambda[i] * vertices[i][2] for i = 1:num_vars),
            base_name = constraint_pre_base_name * "y")

    return formulation_info
end 