export get_variable_bounds,
    get_variable_names, has_geq_constraints, get_geq_constraint_matrices, get_error_bound

"Functions common to AbstractFormulation"
@inline get_variable_bounds(formulation::AbstractFormulation) =
    formulation.lower_bounds, formulation.upper_bounds
@inline get_variable_names(formulation::AbstractFormulation) = formulation.variable_names
@inline has_geq_constraints(formulation::AbstractFormulation) = false
@inline get_geq_constraint_matrices(formulation::AbstractFormulation) =
    Memento.error(_LOGGER, "both the LP and the MILP relaxation have no >= constraints")
@inline get_error_bound(formulation::AbstractFormulation) = formulation.error_bound
