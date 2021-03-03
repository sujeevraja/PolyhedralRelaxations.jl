@testset "test MILP relaxation API" begin
    PR.logger_config!("error")
    m = Model(cbc_optimizer)
    @variable(m, -1.0 <= x <= 1.0)
    @variable(m, y)
    formulation_info =
        construct_univariate_relaxation!(m, a -> a^3, x, y, collect(-1.0:1.0:1.0), true)

    @test size(formulation_info.variables[:delta_1])[1] == 2
    @test size(formulation_info.variables[:delta_2])[1] == 2
    @test size(formulation_info.variables[:z])[1] == 2
    @test haskey(formulation_info.constraints, :x)
    @test haskey(formulation_info.constraints, :y)
    @test haskey(formulation_info.constraints, :first_delta)
    @test haskey(formulation_info.constraints, :below_z)
    @test haskey(formulation_info.constraints, :above_z)
end

@testset "test LP relaxation API" begin
    PR.logger_config!("error")
    m = Model(cbc_optimizer)
    @variable(m, -1.0 <= x <= 1.0)
    @variable(m, y)
    formulation_info =
        construct_univariate_relaxation!(m, a -> a^3, x, y, collect(-1.0:1.0:1.0), false)

    @test size(formulation_info.variables[:lambda])[1] == 4
    @test haskey(formulation_info.constraints, :sum_lambda)
    @test haskey(formulation_info.constraints, :x)
    @test haskey(formulation_info.constraints, :y)
end
