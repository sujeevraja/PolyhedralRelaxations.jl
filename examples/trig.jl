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
    error_tolerances = [NaN64, 1e-1, 1e-2]
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
        err_tol = error_tolerances[i]
        milp = Model(cbc_optimizer)
        @variable(milp, -2.0 <= x <= 5.0)
        @variable(milp, -1.0 <= y[1:5] <= 1.0)

        # construct MILP relaxations
        for j = 1:5
            f = functions[j]
            p = deepcopy(base_partition[j])
            construct_univariate_relaxation!(
                milp,
                f,
                x,
                y[j],
                p,
                true,
                error_tolerance = err_tol,
            )
        end
        @constraint(milp, 5 * y[1] - x <= 0.0)
        @objective(milp, Min, y[2] + y[3] - y[4] - y[5])
        optimize!(milp)
        relaxation_objective = objective_value(milp)
        @test relaxation_objective ≤ best_known_objective
        relative_gap =
            abs(best_known_objective - relaxation_objective) / abs(relaxation_objective)
        if verbose
            println(
                "Optimal solution: $best_known_objective; MILP for error tolerance $err_tol: $relaxation_objective; relative gap: $relative_gap",
            )
        end
    end

    for i = 1:length(error_tolerances)
        err_tol = error_tolerances[i]
        # construct LP relaxations
        lp = Model(cbc_optimizer)
        @variable(lp, -2.0 <= x <= 5.0)
        @variable(lp, -1.0 <= y[1:5] <= 1.0)
        for j = 1:5
            f = functions[j]
            p = deepcopy(base_partition[j])
            construct_univariate_relaxation!(
                lp,
                f,
                x,
                y[j],
                p,
                false,
                error_tolerance = err_tol,
            )
        end
        @constraint(lp, 5 * y[1] - x <= 0.0)
        @objective(lp, Min, y[2] + y[3] - y[4] - y[5])
        optimize!(lp)
        relaxation_objective = objective_value(lp)
        @test relaxation_objective ≤ best_known_objective
        relative_gap =
            abs(best_known_objective - relaxation_objective) / abs(relaxation_objective)
        if verbose
            println(
                "Optimal solution: $best_known_objective; LP for error tolerance $err_tol: $relaxation_objective; relative gap: $relative_gap",
            )
        end
    end

end

example_trig()
