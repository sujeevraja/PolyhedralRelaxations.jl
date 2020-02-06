using ForwardDiff 

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
        @error "the univariate function's domain has to be a closed interval; please specify the bounds using the domain_lb and domain_ub keyword arguments"
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