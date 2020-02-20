@testset "test for errors" begin
    PR.logger_config!("error")
    logger = Memento.getlogger(PolyhedralRelaxations)
    @test_throws(logger, ErrorException, construct_milp_relaxation(x -> x^2, [0.0]))
    @test_throws(logger, ErrorException, construct_milp_relaxation(x -> x^2, [0.0, 1e16]))
    @test_throws(
        logger,
        ErrorException,
        construct_milp_relaxation(x -> tan(x), [0.0, pi / 2])
    )
    @test_throws(
        logger,
        ErrorException,
        construct_milp_relaxation(x -> sin(x), [0.0, 2 * pi])
    )
    @test_throws(
        logger,
        ErrorException,
        construct_milp_relaxation(x -> x^3, [0.0, -1.0, 1.0])
    )
    lp, lp_function_data = construct_lp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))
    @test_throws(logger, ErrorException, get_geq_constraint_matrices(lp))
    @test_throws(logger, ErrorException, get_leq_constraint_matrices(lp))
    @test_throws(logger, ErrorException, construct_milp_relaxation(x -> x^2, [0.0, 1e-16]))
end
