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
    terms = map(it -> (keys(it.indices) |> collect)[1], info)
    max_degree = map(it -> length(it), terms) |> maximum
    (max_degree <= 2) && (return (linking_info = nothing, needed = false))

    # check if terms share more than two variables 
    if !isnothing(degree_limit) && degree_limit < max_degree
        max_degree = degree_limit
    end
    all_variables_in_terms = sort(union(terms...), by = x -> x.index.value)
    linking_constraint_helper = Dict(
        combination => filter(r -> issubset(combination, r), terms) for
        degree in 2:(max_degree-1) for combination in
        Combinatorics.combinations(all_variables_in_terms, degree)
    )
    filter!(r -> length(r.second) >= 2, linking_constraint_helper)
    if isempty(linking_constraint_helper)
        return (linking_constraint_helper = nothing, needed = false)
    end
    return (
        linking_constraint_helper = linking_constraint_helper,
        needed = true,
    )
end

function _build_linking_constraints!(
    m::JuMP.Model,
    info::Vector{FormulationInfo},
    partitions::Dict{JuMP.VariableRef,Vector{T}} where {T<:Real},
    helper::Dict
)::FormulationInfo

    return FormulationInfo()
end
