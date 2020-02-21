f = x -> x^3
partition = collect(-1.0:1.0:1.0)

@testset "sampling tests with x^3 for error tolerance" begin
    PR.silence()
    milp, function_data = construct_milp_relaxation(f, partition, error_tolerance = 1e-2)
    lb, ub = get_domain(function_data)
    var_lb, var_ub = get_variable_bounds(milp)
    num_variables = get_num_variables(milp)
    for i = 1:10
        λ = rand()
        x_val = λ * lb + (1 - λ) * ub
        m = Model(cbc_optimizer)
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
        @test get_error_bound(milp) <= 1e-2
    end
end

@testset "sampling tests with x^3 for binary variable budget" begin
    num_vars = 10
    tol = 1e-5
    for l in [0.1, 0.25, 0.5, 1.0]
        base_partition = collect(-1.0:l:1.0)
        milp, function_data = construct_milp_relaxation(
            f,
            base_partition,
            error_tolerance = tol,
            num_additional_binary_variables = num_vars,
        )
        num_base = length(function_data.base_partition) - 1
        @test PR.get_num_binary_variables(milp) == num_base + num_vars
    end
end

@testset "sampling tests with x^3 for length tolerance" begin
    err_tol = 1e-7
    len_tol = 0.1
    _, function_data = construct_milp_relaxation(
        f,
        partition,
        error_tolerance = err_tol,
        length_tolerance = len_tol,
    )
    for i = 1:length(function_data.partition)-1
        @test function_data.partition[i+1] - function_data.partition[i] >= len_tol
    end
end

@testset "sampling tests with x^3 for derivative tolerance" begin
    err_tol = 1e-7
    der_tol = 0.1
    _, function_data = construct_milp_relaxation(
        f,
        partition,
        error_tolerance = err_tol,
        derivative_tolerance = der_tol,
    )
    for i = 1:length(function_data.partition)-1
        d = function_data.f_dash(function_data.partition[i])
        d_next = function_data.f_dash(function_data.partition[i+1])
        @test abs(d_next - d) >= der_tol
    end
end
