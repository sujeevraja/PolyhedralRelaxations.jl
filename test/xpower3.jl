@testset "x^3 test" begin
    formulation_data, function_data = PR.construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))
    @test formulation_data.x_index == 1
    @test formulation_data.y_index == 2

    lb, ub = formulation_data.lower_bounds, formulation_data.upper_bounds
    @test lb[formulation_data.x_index] == function_data.partition[1]
    @test ub[formulation_data.x_index] == function_data.partition[end]
    @test lb[formulation_data.y_index] == -Inf
    @test ub[formulation_data.y_index] == Inf

    num_variables = PR.get_num_variables(formulation_data)
    @test num_variables == 8

    m = JuMP.Model(glpk_optimizer)
    JuMP.@variable(m, lb[i] <= x[i=1:num_variables] <= ub[i])

    # Add equality constraints
    A_eq, b_eq = formulation_data.A_eq, formulation_data.b_eq
    for i in 1:formulation_data.num_eq_constraints
        eq_indices, eq_values = SparseArrays.findnz(A_eq[i, :])
        JuMP.@constraint(m, dot(x[eq_indices], eq_values) == b_eq[i])
    end

    # Add inequality constraints
    A_leq, b_leq = formulation_data.A_leq, formulation_data.b_leq
    for i in 1:formulation_data.num_leq_constraints
        leq_indices, leq_values = SparseArrays.findnz(A_leq[i, :])
        JuMP.@constraint(m, dot(x[leq_indices], leq_values) <= b_leq[i])
    end

    JuMP.@objective(m, Min, x[formulation_data.x_index])
    JuMP.optimize!(m)
    @test JuMP.objective_value(m) == -1.0

    JuMP.@objective(m, Max, x[formulation_data.x_index])
    JuMP.optimize!(m)
    @test JuMP.objective_value(m) == 1.0
end
