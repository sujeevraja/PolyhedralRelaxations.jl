"""
    _check_if_linking_constraints_are_needed(
        info::Vector{FormulationInfo}, 
        degree_limit::Union{Nothing,T} where {T<:Int64}
    )::NamedTuple

This function checks to see if linking constraints are 
required given the vector of each multilinear terms' 
`FormulationInfo` struct 
"""
function _check_if_linking_constraints_are_needed(
    info::Vector{FormulationInfo},
    degree_limit::Union{Nothing,T} where {T<:Int64},
)::NamedTuple
    # check maximum degree of multilinear terms and return if <= 2
    terms = map(x -> keys(x.indices), info)
    max_degree = map(x -> x |> length, terms) |> maximum
    (max_degree <= 2) && (return (linking_info = nothing, needed = false))
    
    # check if terms share more than two variables 
    if isnothing(degree_limit) || degree_limit < actual_max_degree
        max_degree = degree_limit
    end
    all_variables_in_terms = sort(union(terms...), by = x -> x.index.value)
    linking_info = Dict(
        combination => filter(r -> issubset(combination, r), terms) for
        degree in 2:(max_degree-1) for combination in
        Combinatorics.combinations(all_variables_in_terms, degree)
    )
    filter!(r -> length(r.second) >= 2, link_info)
    (isempty(linking_info)) && (return (linking_info = nothing, needed = false))
    return (linking_info = linking_info, needed = true)
end
