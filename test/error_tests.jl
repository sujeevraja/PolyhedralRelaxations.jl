@testset "test for errors" begin
    PR.logger_config!("error")
    logger = Memento.getlogger(PR)

    optimizer = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    m = Model(optimizer)
    @JuMP.variable(m, -1.0 <= x <= 1.0)
    @JuMP.variable(m, y)

    # Test for too few points in base partition.
    @test_throws(logger, ErrorException,
        PR.construct_univariate_relaxation!(m, a -> a^3, x, y, [0.0], true))

    # Test for unbounded partition point.
    @test_throws(logger, ErrorException,
        PR.construct_univariate_relaxation!(m, a -> a^3, x, y, [0.0, 1e16], true))

    # Test for unbounded function value.
    # @test_throws(logger, ErrorException,
    #   PR.construct_univariate_relaxation!(m, tan, x, y, [0.0, Base.MathConstants.pi / 2.0], true))
        # PR.construct_univariate_relaxation!(m, a -> tan(a), x, y, [0.0, 1.6], true))
    
    # Test for unbounded derivative value.
    # @test_throws(logger, ErrorException,
    #     PR.construct_univariate_relaxation!(m, a -> sqrt(1 - a^2), x, y, [0.0, 1.0], true))

    #= 
    # Test for unbounded derivative value.
    @test_throws(
        logger,
        ErrorException,
        construct_milp_relaxation(x -> sqrt(1 - x^2), [0.0, 1.0])
    )

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
    @test_throws(logger, ErrorException, construct_milp_relaxation(x -> x^2, [0.0, 1e-16])) =#
end
