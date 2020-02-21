using PolyhedralRelaxations, JuMP, Cbc, Test

"""
    example_trig(; verbose = true)

``trig'' is an MINLPLib test instance with trigonometric functions
(see http://www.minlplib.org/trig.html). The best known primal feasible
objective value for this problem is -3.76250036. This example illustrates
how to obtain lower bounds to the optimal objective value using the MILP and
LP relaxations obtained using the PolyhedralRelaxations package.

trig is a non-linear non-convex minimization problem with one decision variable
and one inequality constraint.

Variables: -2 <= x <= 5; start = 1
Constraint: 5 * sin(x) - x <= 0.0
Objective: Minimize sin(11 * x) - cos(13 * x) + sin(17 * x) + cos(19 * x)

We reformulate this problem using additional variables as
Variables:
    -2 <= x <= 5; start = 1
    -1 <= y[1:5] <= 1
Constraints:
    5 * y[1] - x <= 0.0
    y[1] = sin(x)
    y[2] = sin(11 * x)
    y[3] = cos(13 * x)
    y[4] = sin(17 * x)
    y[5] = cos(19 * x)
Objective: Minimize y[2] + y[3] - y[4] - y[5]

"""

function example_trig(; verbose = true)
    best_known_objective = -3.76250036
    cbc_optimizer = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    error_tolerances = [NaN64, 1e-1, 1e-2, 1e-3]
    # base partition holds the inflection points and the end-points for the trigonometric functions
    base_partition = Dict{Int,Vector{Float64}}()
    base_partition[1] = [-2.0, 0.0, π, 5.0]
    base_partition[2] = [-2.0, collect(-7*π/11:π/11:17*π/11)..., 5.0]
    base_partition[3] = [-2.0, collect(-15*π/26:π/13:41*π/26)..., 5.0]
    base_partition[4] = [-2.0, collect(-10*π/17:π/17:27*π/17)..., 5.0]
    base_partition[5] = [-2.0, collect(-23*π/38:π/19:59*π/38)..., 5.0]
    # functions
    functions = Dict{Int,Function}()
    functions[1] = x -> sin(x)
    functions[2] = x -> sin(11 * x)
    functions[3] = x -> cos(13 * x)
    functions[4] = x -> sin(17 * x)
    functions[5] = x -> cos(19 * x)
    for i = 1:length(error_tolerances)
        error_tolerance = error_tolerances[i]
        milp_var = Dict{Int,Any}()
        x_indexes = Dict{Int,Int}()
        y_indexes = Dict{Int,Int}()
        # construct MILP relaxations
        milp = Model(cbc_optimizer)
        for j = 1:5
            var, x_index, y_index = add_milp_relaxation(
                functions[j],
                base_partition[j],
                error_tolerance,
                milp,
                name = string(j),
            )
            milp_var[j] = var
            x_indexes[j] = x_index
            y_indexes[j] = y_index
            (j == 1) && (continue)
            @constraint(milp, var[x_index] == milp_var[1][x_indexes[1]])
        end
        @constraint(milp, 5 * milp_var[1][y_indexes[1]] - milp_var[1][x_indexes[1]] <= 0.0)
        @objective(
            milp,
            Min,
            milp_var[2][y_indexes[2]] + milp_var[3][y_indexes[3]] -
            milp_var[4][y_indexes[4]] - milp_var[5][y_indexes[5]]
        )
        optimize!(milp)
        relaxation_objective = objective_value(milp)
        @test relaxation_objective ≤ best_known_objective
        relative_gap =
            abs(best_known_objective - relaxation_objective) / abs(relaxation_objective)
        if verbose
            println("Optimal solution: $best_known_objective; MILP for $error_tolerance: $relaxation_objective; relative gap: $relative_gap")
        end
    end

    for i = 1:length(error_tolerances)
        error_tolerance = error_tolerances[i]
        lp_var = Dict{Int,Any}()
        x_indexes = Dict{Int,Int}()
        y_indexes = Dict{Int,Int}()
        # construct MILP relaxations
        lp = Model(cbc_optimizer)
        for j = 1:5
            var, x_index, y_index = add_lp_relaxation(
                functions[j],
                base_partition[j],
                error_tolerance,
                lp,
                name = string(j),
            )
            lp_var[j] = var
            x_indexes[j] = x_index
            y_indexes[j] = y_index
            (j == 1) && (continue)
            @constraint(lp, var[x_index] == lp_var[1][x_indexes[1]])
        end
        @constraint(lp, 5 * lp_var[1][y_indexes[1]] - lp_var[1][x_indexes[1]] <= 0.0)
        @objective(
            lp,
            Min,
            lp_var[2][y_indexes[2]] + lp_var[3][y_indexes[3]] - lp_var[4][y_indexes[4]] -
            lp_var[5][y_indexes[5]]
        )
        optimize!(lp)
        relaxation_objective = objective_value(lp)
        @test relaxation_objective ≤ best_known_objective
        relative_gap =
            abs(best_known_objective - relaxation_objective) / abs(relaxation_objective)
        if verbose
            println("Optimal solution: $best_known_objective; LP for $error_tolerance: $relaxation_objective; relative gap: $relative_gap")
        end
    end

end

function add_milp_relaxation(
    f::Function,
    partition::Vector{Float64},
    error_tolerance::Float64,
    m::JuMP.AbstractModel;
    name = "1",
)
    relaxation, _ = construct_milp_relaxation(
        f,
        partition,
        length_tolerance = 1e-5,
        derivative_tolerance = 1e-5,
        error_tolerance = error_tolerance,
    )
    num_variables = get_num_variables(relaxation)
    lb, ub = get_variable_bounds(relaxation)
    names = get_variable_names(relaxation)
    binary = relaxation.binary
    x = @variable(
        m,
        [i = 1:num_variables],
        lower_bound = lb[i],
        upper_bound = ub[i],
        binary = Bool(binary[i]),
        base_name = names[i] * "_" * name
    )
    A, b = get_eq_constraint_matrices(relaxation)
    @constraint(m, A * x .== b)
    A, b = get_leq_constraint_matrices(relaxation)
    @constraint(m, A * x .<= b)
    return x, relaxation.x_index, relaxation.y_index
end

function add_lp_relaxation(
    f::Function,
    partition::Vector{Float64},
    error_tolerance::Float64,
    m::JuMP.AbstractModel;
    name = "1",
)
    relaxation, _ = construct_lp_relaxation(
        f,
        partition,
        length_tolerance = 1e-5,
        derivative_tolerance = 1e-5,
        error_tolerance = error_tolerance,
    )
    num_variables = get_num_variables(relaxation)
    lb, ub = get_variable_bounds(relaxation)
    names = get_variable_names(relaxation)
    x = @variable(
        m,
        [i = 1:num_variables],
        lower_bound = lb[i],
        upper_bound = ub[i],
        base_name = names[i] * "_" * name
    )
    A, b = get_eq_constraint_matrices(relaxation)
    @constraint(m, A * x .== b)
    return x, relaxation.x_index, relaxation.y_index
end


example_trig()
