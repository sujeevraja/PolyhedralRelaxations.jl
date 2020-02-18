@testset "x^3 test" begin
    milp_relaxation, function_data = construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))
    @test milp_relaxation.x_index == 1
    @test milp_relaxation.y_index == 2

    lb, ub = get_variable_bounds(milp_relaxation)
    @test lb[milp_relaxation.x_index] == function_data.partition[1]
    @test ub[milp_relaxation.x_index] == function_data.partition[end]
    @test lb[milp_relaxation.y_index] == -1.0
    @test ub[milp_relaxation.y_index] == 1.0

    num_variables = get_num_variables(milp_relaxation)
    @test num_variables == 8
    @test has_eq_constraints(milp_relaxation) == true 
    @test has_leq_constraints(milp_relaxation) == true 
    @test has_geq_constraints(milp_relaxation) == false 

    # Create variables.
    m = Model(glpk_optimizer)
    @variable(m, lb[i] <= x[i=1:num_variables] <= ub[i],
        binary=Bool(get_variable_type(milp_relaxation)[i]),
        base_name=get_variable_names(milp_relaxation)[i])

    # Add equality constraints
    A, b = get_eq_constraint_matrices(milp_relaxation)
    @constraint(m, A * x .== b)

    # Add inequality constraints
    A, b = get_leq_constraint_matrices(milp_relaxation)
    @constraint(m, A * x .<= b)

    # Test model solution with different objectives.
    @objective(m, Min, x[milp_relaxation.x_index])
    optimize!(m)
    @test objective_value(m) == -1.0

    @objective(m, Max, x[milp_relaxation.x_index])
    optimize!(m)
    @test objective_value(m) == 1.0
end
