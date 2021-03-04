"""
    _build_bilinear_relaxation!(m, x, y, z, partition_x, partition_y)

Build incremental formulation for ``z = xy`` given partition data.
"""
function _build_bilinear_relaxation!(
    m::JuMP.Model,
    x::JuMP.VariableRef,
    y::JuMP.VariableRef,
    z::JuMP.VariableRef,
    partition_x::Vector{<:Real},
    partition_y::Vector{<:Real}
)::FormulationInfo

    return FormulationInfo()
end 