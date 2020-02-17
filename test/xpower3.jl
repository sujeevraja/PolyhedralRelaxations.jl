@testset "x^3 test" begin
    formulation_data, function_data = PR.construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))
    @test formulation_data.x_index == 1
    @test formulation_data.y_index == 2

    lb, ub = formulation_data.lower_bounds, formulation_data.upper_bounds
    @test lb[formulation_data.x_index] == function_data.partition[1]
    @test ub[formulation_data.x_index] == function_data.partition[end]
    @test lb[formulation_data.y_index] == -1.0
    @test ub[formulation_data.y_index] == 1.0

    num_variables = PR.get_num_variables(formulation_data)
    @test num_variables == 8

    # Create variables.
    m = Model(glpk_optimizer)
    @variable(m, lb[i] <= x[i=1:num_variables] <= ub[i],
        binary=Bool(formulation_data.binary[i]),
        base_name=formulation_data.variable_names[i])

    # Add equality constraints
    A_eq, b_eq = formulation_data.A_eq, formulation_data.b_eq
    @constraint(m, A_eq * x .== b_eq)

    # Add inequality constraints
    A_leq, b_leq = formulation_data.A_leq, formulation_data.b_leq
    @constraint(m, A_leq * x .<= b_leq)

    # Test model solution with different objectives.
    @objective(m, Min, x[formulation_data.x_index])
    optimize!(m)
    @test objective_value(m) == -1.0

    @objective(m, Max, x[formulation_data.x_index])
    optimize!(m)
    @test objective_value(m) == 1.0
end
