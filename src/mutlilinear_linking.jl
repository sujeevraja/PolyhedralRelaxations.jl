"""
    _check_if_linking_constraints_are_needed(
        info::Dict{T,FormulationInfo} where {T<:Any}, 
        degree_limit::Union{Nothing,T} where {T<:Int64}
    )::NamedTuple

This function checks to see if linking constraints are 
required given the vector of each multilinear terms' 
`FormulationInfo` struct 
"""
function _check_if_linking_constraints_are_needed(
    info::Dict{T,FormulationInfo} where {T<:Any},
    degree_limit::Union{Nothing,T} where {T<:Int64},
)::NamedTuple
    # check maximum degree of multilinear terms and return if <= 2
    terms = keys(info) |> collect
    max_degree = terms .|> (it -> length(it)) |> maximum
    (max_degree <= 2) && (return (linking_info = nothing, needed = false))

    # check if terms share more than two variables 
    if !isnothing(degree_limit) && degree_limit < max_degree
        max_degree = degree_limit
    end
    all_variables_in_terms = sort(union(terms...), by = x -> x.index.value)
    linking_constraint_helper = Dict(
        tuple(combination...) =>
            filter(r -> issubset(combination, r), terms) for
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
    info::Dict{T,FormulationInfo} where {T<:Any},
    partitions::Dict{JuMP.VariableRef,Vector{T}} where {T<:Real},
    helper::Dict,
)::FormulationInfo
    formulation_info = FormulationInfo()
    # linking variable declaration
    formulation_info.variables[:mu] =
        mu = JuMP.@variable(m, [i in keys(helper)])

    # Add linking constraints to the JuMP model
    for (common_subterm, multilinear_terms) in helper,
        multilinear_term in multilinear_terms

        var_positions =
            common_subterm .|> (it -> findfirst(==(it), multilinear_term))

        indices_iterator =
            info[multilinear_term].indices[multilinear_term][:indices_iterator]
        cartesian_indices = CartesianIndices(indices_iterator |> collect)
        num_lambda = length(indices_iterator)
        lambda = info[multilinear_term].variables[:lambda]

        fcn_value = []
        for j in 1:num_lambda
            cartesian_index = Tuple(cartesian_indices[j])
            pt = [
                partitions[multilinear_term[k]][cartesian_index[k]] for
                k in var_positions
            ]
            push!(fcn_value, prod(pt))
        end

        JuMP.@constraint(
            m,
            mu[common_subterm] ==
            sum(lambda[j] * fcn_value[j] for j in 1:num_lambda)
        )
    end
    return formulation_info
end
