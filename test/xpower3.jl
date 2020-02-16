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

    m = JuMP.Model(glpk_optimizer)
    JuMP.@variable(m, lb[i] <= x[i=1:num_variables] <= ub[i], 
        binary=Bool(formulation_data.binary[i]), 
        base_name=formulation_data.variable_names[i])

    # Add equality constraints
    A_eq, b_eq = formulation_data.A_eq, formulation_data.b_eq
    JuMP.@constraint(m, [i=1:formulation_data.num_eq_constraints], dot(A_eq[i, :], x) == b_eq[i])
    
    # Add inequality constraints
    A_leq, b_leq = formulation_data.A_leq, formulation_data.b_leq
    JuMP.@constraint(m, [i=1:formulation_data.num_leq_constraints], dot(A_leq[i, :], x) <= b_leq[i])
    

    JuMP.@objective(m, Min, x[formulation_data.x_index])
    JuMP.optimize!(m)
    println(m)
    @test JuMP.objective_value(m) == -1.0

    JuMP.@objective(m, Max, x[formulation_data.x_index])
    JuMP.optimize!(m)
    @test JuMP.objective_value(m) == 1.0
end
