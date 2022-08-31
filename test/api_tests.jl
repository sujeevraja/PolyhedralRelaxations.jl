@testset "test MILP relaxation API" begin
    PR.silence!()
    m = Model(milp_optimizer)
    @variable(m, -1.0 <= x <= 1.0)
    @variable(m, y)
    formulation_info = construct_univariate_relaxation!(
        m,
        a -> a^3,
        x,
        y,
        collect(-1.0:1.0:1.0),
        true,
    )

    @test size(formulation_info.variables[:delta_1])[1] == 2
    @test size(formulation_info.variables[:delta_2])[1] == 2
    @test size(formulation_info.variables[:z])[1] == 2
end

@testset "test LP relaxation API" begin
    PR.silence!()
    m = Model(milp_optimizer)
    @variable(m, -1.0 <= x <= 1.0)
    @variable(m, y)
    formulation_info = construct_univariate_relaxation!(
        m,
        a -> a^3,
        x,
        y,
        collect(-1.0:1.0:1.0),
        false,
    )

    @test size(formulation_info.variables[:lambda])[1] == 4
end

@testset "test bilinear relaxation API (McCormick)" begin
    replicates = 10
    tolerance = 1e-3
    for r in 1:replicates
        x_lb, x_ub = 10 * rand(2) .* [-1, 1]
        y_lb, y_ub = 10 * rand(2) .* [-1, 1]

        m = Model(ipopt_optimizer)
        @variable(m, x_lb <= x <= x_ub)
        @variable(m, y_lb <= y <= y_ub)
        @variable(m, z)
        @objective(m, Min, z)
        @constraint(m, x * y == z)
        status = optimize!(m)

        rm = Model(milp_optimizer)
        @variable(rm, x_lb <= x <= x_ub)
        @variable(rm, y_lb <= y <= y_ub)
        @variable(rm, z)
        @objective(rm, Min, z)
        construct_bilinear_relaxation!(rm, x, y, z, [x_lb, x_ub], [y_lb, y_ub])
        rstatus = optimize!(rm)

        @test(objective_value(rm) <= objective_value(m) + tolerance)
        @test(rstatus == status)

        set_objective_sense(m, MOI.MAX_SENSE)
        set_objective_sense(rm, MOI.MAX_SENSE)

        status = optimize!(m)
        rstatus = optimize!(rm)

        @test(objective_value(rm) >= objective_value(m) - tolerance)
        @test(rstatus == status)
    end
end

function refine!(partition, val)
    num_partitions = length(partition) - 1
    val_in_partition = 0
    for i in 1:num_partitions
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
    for _ in 1:20
        x_lb, x_ub = x_p[1], x_p[end]
        y_lb, y_ub = y_p[1], y_p[end]
        epsilon_x = refine!(x_p, point)
        epsilon_y = 10.0
        m = Model(milp_optimizer)
        @variable(m, x_lb <= x <= x_ub)
        @variable(m, y_lb <= y <= y_ub)
        @variable(m, z)
        construct_bilinear_relaxation!(m, x, y, z, x_p, [y_lb, y_ub])
        @constraint(m, x == point)
        @constraint(m, y == point)
        @objective(m, Min, z)
        optimize!(m)
        min_y = objective_value(m)
        set_objective_sense(m, MOI.MAX_SENSE)
        optimize!(m)
        max_y = objective_value(m)
        @test abs(max_y) ≈ epsilon_x * epsilon_y atol = 1e-5
        @test abs(min_y) ≈ epsilon_x * epsilon_y atol = 1e-5
    end
    for _ in 1:20
        x_lb, x_ub = x_p[1], x_p[end]
        y_lb, y_ub = y_p[1], y_p[end]
        epsilon_x = 10.0
        epsilon_y = refine!(y_p, point)
        m = Model(milp_optimizer)
        @variable(m, x_lb <= x <= x_ub)
        @variable(m, y_lb <= y <= y_ub)
        @variable(m, z)
        construct_bilinear_relaxation!(m, x, y, z, [x_lb, x_ub], y_p)
        @constraint(m, x == point)
        @constraint(m, y == point)
        @objective(m, Min, z)
        optimize!(m)
        min_x = objective_value(m)
        set_objective_sense(m, MOI.MAX_SENSE)
        optimize!(m)
        max_x = objective_value(m)
        @test abs(max_x) ≈ epsilon_x * epsilon_y atol = 1e-5
        @test abs(min_x) ≈ epsilon_x * epsilon_y atol = 1e-5
    end
end
