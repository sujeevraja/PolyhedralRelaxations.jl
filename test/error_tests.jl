@testset "test for errors" begin
    PR.logger_config!("error")
    logger = Memento.getlogger(PR)

    m = Model(cbc_optimizer)
    @variable(m, -1.0 <= x <= 1.0)
    @variable(m, y)

    # Test for too few points in base partition.
    @test_throws(
        logger,
        ErrorException,
        construct_univariate_relaxation!(m, a -> a^3, x, y, [0.0], true)
    )

    # Test for unbounded partition point.
    @test_throws(
        logger,
        ErrorException,
        construct_univariate_relaxation!(m, a -> a^3, x, y, [0.0, 1e16], true)
    )

    # Test for unbounded function value.
    @test_throws(
        logger,
        ErrorException,
        construct_univariate_relaxation!(m, tan, x, y, [0.0, π / 2.0], true)
    )

    # Test for unbounded derivative value.
    @test_throws(
        logger,
        ErrorException,
        construct_univariate_relaxation!(m, a -> sqrt(1 - a^2), x, y, [0.0, 1.0], true)
    )

    # Test for unbounded derivative value.
    @test_throws(
        logger,
        ErrorException,
        construct_univariate_relaxation!(m, a -> sqrt(1 - a^2), x, y, [0.0, 1.0], true)
    )

    # Test for equal derivative values in adjacent points of the base partition.
    @test_throws(
        logger,
        ErrorException,
        construct_univariate_relaxation!(m, sin, x, y, [0.0, 2 * π], true)
    )

    # Test for invalid ordering of points in base partition.
    @test_throws(
        logger,
        ErrorException,
        construct_univariate_relaxation!(m, a -> a^3, x, y, [0.0, -1.0, 1.0], true)
    )

    @test_throws(
        logger,
        ErrorException,
        construct_univariate_relaxation!(m, a -> a^2, x, y, [0.0, 1e-16], true)
    )
end
