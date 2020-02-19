get_milp_xcube(; error_tolerance = 1e-2) = construct_milp_relaxation(
    x -> x^3,
    collect(-1.0:1.0:1.0),
    error_tolerance = error_tolerance,
)

@testset "sampling tests for error tolerance" begin
    PR.logger_config!("error")
    milp, milp_function_data = get_milp_xcube()
    lb, ub = get_domain(milp_function_data)
    var_lb, var_ub = get_variable_bounds(milp)
    num_variables = get_num_variables(milp)
    for i = 1:10
        λ = rand()
        x_val = λ * lb + (1 - λ) * ub
        m = Model(glpk_optimizer)
        @variable(
            m,
            var_lb[i] <= x[i = 1:num_variables] <= var_ub[i],
            binary = Bool(get_variable_type(milp)[i]),
            base_name = get_variable_names(milp)[i]
        )
        A, b = get_eq_constraint_matrices(milp)
        @constraint(m, A * x .== b)
        A, b = get_leq_constraint_matrices(milp)
        @constraint(m, A * x .<= b)
        @constraint(m, x[milp.x_index] == x_val)
        # solve for min y
        @objective(m, Min, x[milp.y_index])
        optimize!(m)
        y_min = objective_value(m)
        # solve for max y 
        @objective(m, Max, x[milp.y_index])
        optimize!(m)
        y_max = objective_value(m)
        @test abs(y_max - y_min) <= 1e-2
    end

end
