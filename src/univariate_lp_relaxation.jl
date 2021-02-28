"""
    get_lp_relaxation_vertices(univariate_function_data::UnivariateFunctionData)::Vector{Vertex}

Return all the vertices that are a part of the LP relaxation given the
function data.
"""
function get_lp_relaxation_vertices(
    univariate_function_data::UnivariateFunctionData,
)::Vector{Vertex}
    secant_vertices, tangent_vertices = collect_vertices(univariate_function_data)
    vertices = Vertex[]
    push!(vertices, secant_vertices[1])
    append!(vertices, tangent_vertices)
    push!(vertices, secant_vertices[end])
    return vertices
end

"""
    build_lp_relaxation(univariate_function_data)

Build LP relaxation for ``y=f(x)`` given the univariate function data.
"""
function build_univariate_lp_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    univariate_function_data::UnivariateFunctionData,
)::Pair{LPRelaxation,UnivariateFunctionData}
    vertices = get_lp_relaxation_vertices(univariate_function_data)
    f_min, f_max = Inf, -Inf
    for v in vertices
        f_val = v[2]
        f_min = min(f_min, f_val)
        f_max = max(f_max, f_val)
    end
    num_vertices = length(vertices)

    lp_formulation_info = FormulationInfo()

    # add variables 
    num_lambda_variables = length(univariate_function_data) + 1
    @assert num_vertices == num_lambda_variables
    lambda =
        lp_formulation_info.var[:lambda] =
            @variable(m, 0 <= [i = 1:num_lambda_variables] <= 1)

    # add constraints 
    lp_formulation_info.con[:sum_lambda] = @constraint(m, sum(lambda) == 1)
    lp_formulation_info.con[:x] =
        @constraint(m, x == sum(lambda[i] * vertices[i][1] for i = 1:num_vertices))
    lp_formulation_info.con[:y] =
        @constraint(m, y == sum(lambda[i] * vertices[i][2] for i = 1:num_vertices))

    return lp_formulation_info
end
