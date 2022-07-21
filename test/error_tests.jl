@testset "test for errors" begin
    PR.silence!()

    m = Model(milp_optimizer)
    @variable(m, -1.0 <= x <= 1.0)
    @variable(m, y)

    # Test for too few points in base partition.
    @test_throws ErrorException construct_univariate_relaxation!(
        m,
        a -> a^3,
        x,
        y,
        [0.0],
        true,
    )

    # Test for unbounded partition point.
    @test_throws ErrorException construct_univariate_relaxation!(
        m,
        a -> a^3,
        x,
        y,
        [0.0, 1e16],
        true,
    )

    # # Test for unbounded function value.
    @test_throws ErrorException construct_univariate_relaxation!(
        m,
        tan,
        x,
        y,
        [0.0, π / 2.0],
        true,
    )

    # # Test for unbounded derivative value.
    @test_throws ErrorException construct_univariate_relaxation!(
        m,
        a -> sqrt(1 - a^2),
        x,
        y,
        [0.0, 1.0],
        true,
    )

    # # Test for unbounded derivative value.
    @test_throws ErrorException construct_univariate_relaxation!(
        m,
        a -> sqrt(1 - a^2),
        x,
        y,
        [0.0, 1.0],
        true,
    )

    # # Test for equal derivative values in adjacent points of the base partition.
    @test_throws ErrorException construct_univariate_relaxation!(
        m,
        sin,
        x,
        y,
        [0.0, 2 * π],
        true,
    )

    # # Test for invalid ordering of points in base partition.
    @test_throws ErrorException construct_univariate_relaxation!(
        m,
        a -> a^3,
        x,
        y,
        [0.0, -1.0, 1.0],
        true,
    )

    @test_throws ErrorException construct_univariate_relaxation!(
        m,
        a -> a^2,
        x,
        y,
        [0.0, 1e-16],
        true,
    )

    x_lb, x_ub = 10 * rand(2) .* [-1, 1]
    y_lb, y_ub = 10 * rand(2) .* [-1, 1]
    x_mid = (x_lb + x_ub) / 2.0
    y_mid = (y_lb + y_ub) / 2.0
    m = JuMP.Model(milp_optimizer)
    JuMP.@variable(m, x_lb <= x <= x_ub)
    JuMP.@variable(m, y_lb <= y <= y_ub)
    JuMP.@variable(m, z)
    @test_throws ErrorException construct_bilinear_relaxation!(
        m,
        x,
        y,
        z,
        [x_lb, x_mid, x_ub],
        [y_lb, y_mid, y_ub],
    )
end
