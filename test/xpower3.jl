@testset "x^3 test" begin
    milp_relaxation, function_data = construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))
    @test milp_relaxation.x_index == 1
    @test milp_relaxation.y_index == 2

    lb, ub = milp_relaxation.lower_bounds, milp_relaxation.upper_bounds
    @test lb[milp_relaxation.x_index] == function_data.partition[1]
    @test ub[milp_relaxation.x_index] == function_data.partition[end]
    @test lb[milp_relaxation.y_index] == -1.0
    @test ub[milp_relaxation.y_index] == 1.0

    num_variables = get_num_variables(milp_relaxation)
    @test num_variables == 8

    # Create variables.
    m = Model(glpk_optimizer)
    @variable(m, lb[i] <= x[i=1:num_variables] <= ub[i],
        binary=Bool(milp_relaxation.binary[i]),
        base_name=milp_relaxation.variable_names[i])

    # Add equality constraints
    A_eq, b_eq = get_eq_constraint_matrices(milp_relaxation)
    @constraint(m, A_eq * x .== b_eq)

    # Add inequality constraints
    A_leq, b_leq = get_leq_constraint_matrices(milp_relaxation)
    @constraint(m, A_leq * x .<= b_leq)

    # Test model solution with different objectives.
    @objective(m, Min, x[milp_relaxation.x_index])
    optimize!(m)
    @test objective_value(m) == -1.0

    @objective(m, Max, x[milp_relaxation.x_index])
    optimize!(m)
    @test objective_value(m) == 1.0
end
