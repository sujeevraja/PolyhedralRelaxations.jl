"""
    _get_lp_relaxation_vertices(univariate_function_data::UnivariateFunctionData)::Vector{Vertex}

Return vertices of the LP relaxation of the given function.
"""
function _get_lp_relaxation_vertices(
    univariate_function_data::UnivariateFunctionData,
)::Vector{Vertex}
    sec_vs, tan_vs = collect_vertices(univariate_function_data)
    vertices = Vertex[]
    push!(vertices, sec_vs[1])
    append!(vertices, tan_vs)
    push!(vertices, sec_vs[end])
    return vertices
end

"""
    build_univariate_lp_relaxation(univariate_function_data)

Build LP relaxation for ``y=f(x)`` given the univariate function data.
"""
function build_univariate_lp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    univariate_function_data::UnivariateFunctionData,
)::FormulationInfo
    vertices = _get_lp_relaxation_vertices(univariate_function_data)
    num_vars = length(vertices)
    formulation_info = FormulationInfo()

    # add variables 
    @variable(m, 0 <= lambda[1:num_vars] <= 1)
    formulation_info.variables[:lambda] = lambda

    # add constraints 
    formulation_info.constraints[:sum_lambda] = @constraint(m, sum(lambda) == 1)
    formulation_info.constraints[:x] = @constraint(m,
        x == sum(lambda[i] * vertices[i][1] for i = 1:num_vars))
    formulation_info.constraints[:y] = @constraint(m,
        y == sum(lambda[i] * vertices[i][2] for i = 1:num_vars))

    return formulation_info
end
