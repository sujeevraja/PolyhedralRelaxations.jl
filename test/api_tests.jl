@testset "test MILP relaxation API" begin
    PR.silence!()
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
    PR.silence!()
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


@testset "test bilinear relaxation API (McCormick)" begin
    replicates = 10
    tolerance = 1e-3
    for r = 1:replicates
        x_lb, x_ub = 10 * rand(2) .* [-1, 1]
        y_lb, y_ub = 10 * rand(2) .* [-1, 1]

        m = JuMP.Model(ipopt_optimizer)
        JuMP.@variable(m, x_lb <= x <= x_ub)
        JuMP.@variable(m, y_lb <= y <= y_ub)
        JuMP.@variable(m, z)
        JuMP.@objective(m, Min, z)
        JuMP.@constraint(m, x * y == z)
        status = JuMP.optimize!(m)

        rm = JuMP.Model(cbc_optimizer)
        JuMP.@variable(rm, x_lb <= x <= x_ub)
        JuMP.@variable(rm, y_lb <= y <= y_ub)
        JuMP.@variable(rm, z)
        JuMP.@objective(rm, Min, z)
        construct_bilinear_relaxation!(rm, x, y, z, [x_lb, x_ub], [y_lb, y_ub])
        rstatus = JuMP.optimize!(rm)

        @test(JuMP.objective_value(rm) <= JuMP.objective_value(m) + tolerance)
        @test(rstatus == status)

        JuMP.set_objective_sense(m, MOI.MAX_SENSE)
        JuMP.set_objective_sense(rm, MOI.MAX_SENSE)

        status = JuMP.optimize!(m)
        rstatus = JuMP.optimize!(rm)

        @test(JuMP.objective_value(rm) >= JuMP.objective_value(m) - tolerance)
        @test(rstatus == status)
    end
end

function refine!(partition, val)
    num_partitions = length(partition) - 1
    val_in_partition = 0
    for i = 1:num_partitions
        if val > partition[i] && val < partition[i+1]
            val_in_partition = i
        end
    end
    left = val - (val - partition[val_in_partition]) * 0.5
    right = val + (partition[val_in_partition+1] - val) * 0.5
    insert!(partition, val_in_partition + 1, left)
    insert!(partition, val_in_partition + 2, right)
    return (val - partition[val_in_partition]) * 0.5
end

@testset "test bilinear MIP relaxation API" begin
    gap_x = Float64[]
    gap_y = Float64[]
    point = 0.0

    x_p = [-10.0, 10.0]
    y_p = [-10.0, 10.0]
    for _ = 1:20
        x_lb, x_ub = x_p[1], x_p[end]
        y_lb, y_ub = y_p[1], y_p[end]
        epsilon_x = refine!(x_p, point)
        epsilon_y = 10.0
        m = JuMP.Model(cbc_optimizer)
        JuMP.@variable(m, x_lb <= x <= x_ub)
        JuMP.@variable(m, y_lb <= y <= y_ub)
        JuMP.@variable(m, z)
        construct_bilinear_relaxation!(m, x, y, z, x_p, [y_lb, y_ub])
        JuMP.@constraint(m, x == point)
        JuMP.@constraint(m, y == point)
        JuMP.@objective(m, Min, z)
        JuMP.optimize!(m)
        min_y = JuMP.objective_value(m)
        JuMP.set_objective_sense(m, MOI.MAX_SENSE)
        JuMP.optimize!(m)
        max_y = JuMP.objective_value(m)
        @test abs(max_y) ≈ epsilon_x * epsilon_y atol = 1e-5
        @test abs(min_y) ≈ epsilon_x * epsilon_y atol = 1e-5
    end
    for _ = 1:20
        x_lb, x_ub = x_p[1], x_p[end]
        y_lb, y_ub = y_p[1], y_p[end]
        epsilon_x = 10.0
        epsilon_y = refine!(y_p, point)
        m = JuMP.Model(cbc_optimizer)
        JuMP.@variable(m, x_lb <= x <= x_ub)
        JuMP.@variable(m, y_lb <= y <= y_ub)
        JuMP.@variable(m, z)
        construct_bilinear_relaxation!(m, x, y, z, [x_lb, x_ub], y_p)
        JuMP.@constraint(m, x == point)
        JuMP.@constraint(m, y == point)
        JuMP.@objective(m, Min, z)
        JuMP.optimize!(m)
        min_x = JuMP.objective_value(m)
        JuMP.set_objective_sense(m, MOI.MAX_SENSE)
        JuMP.optimize!(m)
        max_x = JuMP.objective_value(m)
        @test abs(max_x) ≈ epsilon_x * epsilon_y atol = 1e-5
        @test abs(min_x) ≈ epsilon_x * epsilon_y atol = 1e-5
    end
end
