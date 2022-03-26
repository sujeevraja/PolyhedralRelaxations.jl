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
    _build_univariate_lp_relaxation!(m, x, y, univariate_function_data, pre_base_name)

Build LP relaxation for ``y=f(x)`` given the univariate function data.
"""
function _build_univariate_lp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    univariate_function_data::UnivariateFunctionData,
    pre_base_name::AbstractString,
)::FormulationInfo
    vertices = _get_lp_relaxation_vertices(univariate_function_data)
    num_vars = length(vertices)
    formulation_info = FormulationInfo()

    # add variables 
    lambda =
        formulation_info.variables[:lambda] = @variable(
            m,
            [1:num_vars],
            lower_bound = 0.0,
            upper_bound = 1.0,
            base_name = pre_base_name * "_lambda"
        )
    formulation_info.variables[:lambda] = lambda

    # add constraints 
    formulation_info.constraints[:sum_lambda] = @constraint(m, sum(lambda) == 1)
    formulation_info.constraints[:x] =
        @constraint(m, x == sum(lambda[i] * vertices[i][1] for i in 1:num_vars))
    formulation_info.constraints[:y] =
        @constraint(m, y == sum(lambda[i] * vertices[i][2] for i in 1:num_vars))

    return formulation_info
end
