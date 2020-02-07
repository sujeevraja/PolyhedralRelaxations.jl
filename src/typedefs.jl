using SparseArrays

struct UnivariateFunction
    f::Function
    f_dash::Function
    domain_lb::Real
    domain_ub::Real
    inflection_points::Vector{Real}
end

"Constructors for UnivariateFunction"
function UnivariateFunction(f::Function, f_dash::Function; domain_lb::Real=-Inf, domain_ub::Real=Inf, inflection_points::Vector{Real}=[])::UnivariateFunction
    if isinf(domain_lb) || isinf(domain_ub)
        Memento.error(_LOGGER, "the univariate function's domain has to be a closed interval; please specify the bounds using the domain_lb and domain_ub keyword arguments")
    end
    return UnivariateFunction(f, f_dash, domain_lb, domain_ub, inflection_points)
end
UnivariateFunction(f::Function; domain_lb::Real=-Inf, domain_ub::Real=Inf, inflection_points::Vector{Real}=[]) =
    UnivariateFunction(f, x -> ForwardDiff.derivative(f, x), domain_lb=domain_lb, domain_ub=domain_ub, inflection_points=inflection_points)
UnivariateFunction(f::Function, f_dash::Function, lb, ub) =  UnivariateFunction(f, f_dash, domain_lb=lb, domain_ub=ub)

"Getters for UnivariateFunction"
@inline get_function(uf::UnivariateFunction)::Function = uf.f
@inline get_derivative(uf::UnivariateFunction)::Function = uf.f_dash
@inline get_domain_lb(uf::UnivariateFunction)::Real = uf.domain_lb
@inline get_domain_ub(uf::UnivariateFunction)::Real = uf.domain_ub
@inline get_domain(uf::UnivariateFunction)::Tuple{Real,Real} = get_domain_lb(uf), get_domain_ub(uf)
@inline get_inflection_points(uf::UnivariateFunction)::Vector{<:Real} = uf.inflection_points

struct VariableInfo
    num_variables::Int64
    variable_bounds::Vector{Tuple{Real,Real}}
    independent_variable_id::Int64
    dependent_variable_id::Int64
    partitions::Vector{Tuple{Real,Real}}
    indicator_variable_ids::Union{Int64, Nothing}
    continuous_variable_ids::Vector{Tuple{Int64,Int64}}
end

"Getters for VariableInfo"
@inline get_num_variables(vi::VariableInfo)::Int64 = vi.num_variables
@inline get_variable_bounds(vi::VariableInfo)::Vector{Tuple{<:Real,<:Real}} = vi.variable_bounds
@inline get_variable_bounds(vi::VariableInfo, id::Int64)::Tuple{<:Real,<:Real} = vi.variable_bounds[id]
@inline get_independent_variable_id(vi::VariableInfo)::Int64 = vi.independent_variable_id
@inline get_dependent_variable_id(vi::VariableInfo)::Int64 = vi.dependent_variable_id
@inline get_partitions(vi::VariableInfo)::Vector{Tuple{<:Real,<:Real}} = vi.partitions
@inline get_num_partitions(vi::VariableInfo)::Int64 = length(vi.partitions)
@inline get_indicator_variable_ids(vi::VariableInfo)::Vector{Int64, Nothing} = vi.indicator_variable_ids
@inline get_continuous_variable_ids(vi::VariableInfo)::Vector{Tuple{Int64,Int64}} = vi.continuous_variable_ids

"Vertices in 2D used to build the polyhedral relaxation sequence"
struct Vertex
    x::Real
    y::Real
end

"""
Constraint coefficients and right-hand-side of MIP relaxation.

Variables are ordered as: x, y, delta_1^i, delta_2^i, z_i
"""
struct Model
    A::SparseMatrixCSC{Real,Int64}
    b::SparseVector{Real,Int64}
    x_index::Int64
    y_index::Int64
    delta_1_indices::Array{Int64,1}
    delta_2_indices::Array{Int64,1}
    z_indices::Array{Int64,1}
end
