"""
    get_lp_relaxation_vertices(univariate_function_data::UnivariateFunctionData)::Vector{Vertex}

Returns all the vertices that are a part of the LP relaxation given the function data
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

Function to build the LP relaxation for ``y=f(x)`` given the univariate function data
"""
function build_univariate_lp_relaxation(
    m::JuMP.Model, x::JuMP.VariableRef, y::JuMP.VariableRef,
    univariate_function_data::UnivariateFunctionData,
)::Pair{LPRelaxation,UnivariateFunctionData}
    vertices = get_lp_relaxation_vertices(univariate_function_data)
    f_min, f_max = Inf, -Inf
    for v in vertices
        f_val = v[2]
        f_min = min(f_min, f_val)
        f_max = max(f_max, f_val)
    end

    lp_formulation_info = FormulationInfo()
    
    # Add the LP relaxation formulation 
    return lp_formulation_info
end
