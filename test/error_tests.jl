@testset "test for errors" begin
    PR.logger_config!("error")
    logger = Memento.getlogger(PolyhedralRelaxations)

    # Test for too few points in base partition.
    @test_throws(logger, ErrorException, construct_milp_relaxation(x -> x^2, [0.0]))

    # Test for unbounded partition point.
    @test_throws(logger, ErrorException, construct_milp_relaxation(x -> x^2, [0.0, 1e16]))

    # Test for unbounded function value.
    @test_throws(
        logger,
        ErrorException,
        construct_milp_relaxation(x -> tan(x), [0.0, pi / 2])
    )

    # Test for unbounded derivative value.
    @test_throws(logger, ErrorException, construct_milp_relaxation(x -> sqrt(1-x^2), [0.0, 1e-16]))

    # Test for equal derivative values in adjacent points of the base partition.
    @test_throws(
        logger,
        ErrorException,
        construct_milp_relaxation(x -> sin(x), [0.0, 2 * pi])
    )

    # Test for invalid ordering of points in base partition.
    @test_throws(
        logger,
        ErrorException,
        construct_milp_relaxation(x -> x^3, [0.0, -1.0, 1.0])
    )

    # Test for getting unavailable constraint types.
    lp, lp_function_data = construct_lp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))
    @test_throws(logger, ErrorException, get_geq_constraint_matrices(lp))
    @test_throws(logger, ErrorException, get_leq_constraint_matrices(lp))
    @test_throws(logger, ErrorException, construct_milp_relaxation(x -> x^2, [0.0, 1e-16]))
end
