"""
    _build_multilinear_convex_hull_relaxation!(m, x, z, partitions, pre_base_name)

Build a convex hull formulation LP for  for ``z = prod(x)`` given partition data. 
"""
function _build_multilinear_convex_hull_relaxation!(
    m::JuMP.Model,
    x::Tuple,
    z::JuMP.VariableRef,
    partitions::Dict{JuMP.VariableRef,Vector{T}} where {T<:Real},
    pre_base_name::AbstractString,
)::FormulationInfo
    formulation_info = FormulationInfo()
    dim = length(x)
    index_range = 1:2
    index_ranges = [index_range for _ in 1:dim]
    indices_iterator = Iterators.product(index_ranges...)

    cartesian_indices = CartesianIndices(indices_iterator |> collect)
    formulation_info.indices[x] = Dict(
        :index_ranges => index_ranges,
        :indices_iterator => indices_iterator,
    )
    num_lambda = length(indices_iterator)

    extreme_points = []
    for i in 1:num_lambda
        cartesian_index = Tuple(cartesian_indices[i])
        pt = [partitions[var][cartesian_index[k]] for (k, var) in enumerate(x)]
        push!(pt, prod(pt))
        push!(extreme_points, pt)
    end

    lambda =
        formulation_info.variables[:lambda] = JuMP.@variable(
            m,
            [i in 1:num_lambda],
            lower_bound = 0.0,
            upper_bound = 1.0,
            base_name = pre_base_name * "lambda"
        )

    # convex hull constraints
    JuMP.@constraint(m, sum(lambda) == 1)
    JuMP.@constraint(
        m,
        [x..., z] .==
        sum([lambda[i] * extreme_points[i] for i in 1:num_lambda])
    )

    return formulation_info
end

"""
    _get_slice_indices(var_id::Int64, slice_id::Int64, cartesian_indices)::Vector
    
Helper function to get all the linear indices corresponding to a variable and slice 
"""
function _get_slice_indices(
    var_id::Int64,
    slice_id::Int64,
    cartesian_indices,
)::Vector
    shape = size(cartesian_indices)
    dim = length(shape)
    slice = [(i == var_id) ? slice_id : Base.Colon() for i in 1:dim]
    return cartesian_indices[slice...] |>
           vec |>
           x -> LinearIndices(shape)[CartesianIndex.(x)]
end

"""
    _build_multilinear_milp_relaxation!(m, x, z, partitions, pre_base_name)

Build a piecewise polyhedral relaxation for  for ``z = prod(x)`` given partition data. 
"""
function _build_multilinear_sos2_relaxation!(
    m::JuMP.Model,
    x::Tuple,
    z::JuMP.VariableRef,
    partitions::Dict{JuMP.VariableRef,Vector{T}} where {T<:Real},
    pre_base_name::AbstractString,
)::FormulationInfo
    formulation_info = FormulationInfo()
    index_ranges = [1:length(partitions[var]) for var in x]
    indices_iterator = Iterators.product(index_ranges...)

    cartesian_indices = CartesianIndices(indices_iterator |> collect)
    formulation_info.indices[x] = Dict(
        :index_ranges => index_ranges,
        :indices_iterator => indices_iterator,
    )
    num_lambda = length(indices_iterator)
    extreme_points = []
    for i in 1:num_lambda
        cartesian_index = Tuple(cartesian_indices[i])
        pt = [partitions[var][cartesian_index[k]] for (k, var) in enumerate(x)]
        push!(pt, prod(pt))
        push!(extreme_points, pt)
    end

    lambda =
        formulation_info.variables[:lambda] = JuMP.@variable(
            m,
            [i in 1:num_lambda],
            lower_bound = 0.0,
            upper_bound = 1.0,
            base_name = pre_base_name * "lambda"
        )

    bin = formulation_info.variables[:bin] = Dict{JuMP.VariableRef,Any}()
    for var in x
        formulation_info.variables[:bin][var] = JuMP.@variable(
            m,
            [j in 1:(length(partitions[var])-1)],
            binary = true
        )
    end

    # multiplier and binary summation constraints
    JuMP.@constraint(m, sum(lambda) == 1)
    formulation_info.constraints[:bin] = Dict{JuMP.VariableRef,Any}()
    for var in x
        formulation_info.constraints[:bin][var] =
            JuMP.@constraint(m, sum(bin[var]) == 1)
    end

    # convex combination constraints
    JuMP.@constraint(
        m,
        [x..., z] .==
        sum([lambda[i] * extreme_points[i] for i in 1:num_lambda])
    )

    # SOS-2 constraints 
    for (variable_id, var) in enumerate(x)
        for slice_id in index_ranges[variable_id]
            sliced_indices =
                _get_slice_indices(variable_id, slice_id, cartesian_indices)
            if slice_id == 1
                JuMP.@constraint(
                    m,
                    sum(lambda[sliced_indices]) <= first(bin[var])
                )
            elseif slice_id == length(index_ranges[variable_id])
                JuMP.@constraint(
                    m,
                    sum(lambda[sliced_indices]) <= last(bin[var])
                )
            else
                JuMP.@constraint(
                    m,
                    sum(lambda[sliced_indices]) <=
                    sum(bin[var][slice_id-1:slice_id])
                )
            end
        end
    end

    return formulation_info
end
