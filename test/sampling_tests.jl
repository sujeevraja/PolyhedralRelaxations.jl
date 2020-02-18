get_milp_xcube(;error_tolerance=1e-2) = 
    construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0), error_tolerance = error_tolerance)
get_lp_xcube(;error_tolerance=1e-2) = 
    construct_lp_relaxation(x -> x^3, collect(-1.0:1.0:1.0), error_tolerance = error_tolerance)

@testset "sampling tests for error tolerance" begin
    PR.logger_config!("error")
    milp, milp_function_data = get_milp_xcube() 
    lp, lp_function_data = get_lp_xcube()
    @test get_domain(lp_function_data) == get_domain(milp_function_data)
    lb, ub = get_domain(milp_function_data)
    for i=1:10
        # generate x in domain 
        # solve max (min) y 
        # take difference 
        # assert difference is less than error
    end 

end 