f = x -> x^3
f_dash = x -> 3 * x^2
partition = collect(-1.0:1.0:1.0)

@testset "sampling tests with x^3 for error tolerance" begin
    PR.silence!()
    p = deepcopy(partition)
    for i in 1:10
        lb = -1.0
        ub = 1.0
        λ = rand()
        x_val = λ * lb + (1 - λ) * ub
        m = Model(cbc_optimizer)
        @variable(m, -1.0 <= x <= 1.0)
        @variable(m, y)
        construct_univariate_relaxation!(
            m,
            f,
            x,
            y,
            p,
            true,
            error_tolerance = 1e-2,
        )
        @constraint(m, x == x_val)

        # solve for min y
        @objective(m, Min, y)
        optimize!(m)
        y_min = objective_value(m)

        # solve for max y
        @objective(m, Max, y)
        optimize!(m)
        y_max = objective_value(m)
        @test abs(y_max - y_min) <= 1e-2
    end
end

@testset "sampling tests with x^3 for binary variable budget" begin
    num_vars = 10
    tol = 1e-5
    for l in [0.1, 0.25, 0.5, 1.0]
        base_partition = collect(-1.0:l:1.0)
        num_base = length(base_partition) - 1
        m = Model(cbc_optimizer)
        @variable(m, -1.0 <= x <= 1.0)
        @variable(m, y)
        formulation_info = construct_univariate_relaxation!(
            m,
            f,
            x,
            y,
            base_partition,
            true,
            error_tolerance = tol,
            num_additional_partitions = num_vars,
        )
        @test length(formulation_info.variables[:z]) == num_base + num_vars
    end
end

@testset "sampling tests with x^3 for length tolerance" begin
    err_tol = 1e-7
    len_tol = 0.1
    p = deepcopy(partition)
    m = Model(cbc_optimizer)
    @variable(m, -1.0 <= x <= 1.0)
    @variable(m, y)
    construct_univariate_relaxation!(
        m,
        f,
        x,
        y,
        p,
        true,
        error_tolerance = err_tol,
        length_tolerance = len_tol,
    )
    for i in 1:length(p)-1
        @test p[i+1] - p[i] >= len_tol
    end
end

@testset "sampling tests with x^3 for derivative tolerance" begin
    err_tol = 1e-7
    der_tol = 0.1
    p = deepcopy(partition)
    m = Model(cbc_optimizer)
    @variable(m, -1.0 <= x <= 1.0)
    @variable(m, y)
    construct_univariate_relaxation!(
        m,
        f,
        x,
        y,
        p,
        true,
        error_tolerance = err_tol,
        derivative_tolerance = der_tol,
    )

    for i in 1:length(p)-1
        d = f_dash(p[i])
        d_next = f_dash(p[i+1])
        @test abs(d_next - d) >= der_tol
    end
end
