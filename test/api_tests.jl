@testset "test lp relaxation" begin
    PR.logger_config!("error")
    optimizer = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    m = Model(optimizer)
    @JuMP.variable(m, -1.0 <= x <= 1.0)
    @JuMP.variable(m, y)
    formulation_info = PR.construct_univariate_relaxation!(
        m, a -> a^3, x, y, collect(-1.0:1.0:1.0), false)

    @test size(formulation_info.variables[:lambda])[1] == 4
    @test haskey(formulation_info.constraints, :sum_lambda)
    @test haskey(formulation_info.constraints, :x)
    @test haskey(formulation_info.constraints, :y)
end

@testset "api tests" begin
    PR.logger_config!("error")
    optimizer = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
    m = Model(optimizer)
    @JuMP.variable(m, -1.0 <= x <= 1.0)
    @JuMP.variable(m, y)
    lp_formulation_info = PR.construct_univariate_relaxation!(
        m, a -> a^3, x, y, collect(-1.0:1.0:1.0), false)

    # milp_relaxation, function_data = construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))

    #= 
    @test get_domain(function_data) == (-1.0, 1.0)
    @test get_derivative(function_data)(0.5) == 0.75
    @test get_function(function_data)(0.5) == 0.125
    @test get_partition(function_data) == [-1.0, 0.0, 1.0]

    @test milp_relaxation.x_index == 1
    @test milp_relaxation.y_index == 2
    @test get_num_variables(milp_relaxation) == 8
    @test has_eq_constraints(milp_relaxation) == true
    @test has_leq_constraints(milp_relaxation) == true
    @test has_geq_constraints(milp_relaxation) == false
    @test milp_relaxation.num_eq_constraints == 2
    @test milp_relaxation.num_leq_constraints == 3
    @test milp_relaxation.delta_1_indices == [3, 4]
    @test milp_relaxation.delta_2_indices == [5, 6]
    @test milp_relaxation.z_indices == [7, 8]

    lp_relaxation, function_data = construct_lp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))

    @test get_num_variables(lp_relaxation) == 6
    @test has_eq_constraints(lp_relaxation) == true
    @test has_leq_constraints(lp_relaxation) == false
    @test has_geq_constraints(lp_relaxation) == false
    @test lp_relaxation.num_constraints == 3
    @test lp_relaxation.λ_indices == [3, 4, 5, 6] =#
end
