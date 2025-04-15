@testset "linear objective function tests" begin
    PR.silence!()
    f, partition = a -> a^3, collect(-1.0:0.25:1.0)
    ufd = UnivariateFunctionData(
        f,
        a -> 3 * a^2,
        partition,
        NaN64,
        1e-3,
        1e-3,
        0,
        length(partition),
    )

    # create LP relaxation of the MILP relaxation of the univariate function.
    milp = Model(milp_optimizer)
    @variable(milp, -1.0 <= x <= 1.0)
    @variable(milp, y)
    construct_univariate_relaxation!(milp, f, x, y, partition, true)

    # create LP relaxation of the univariate function using the convex hull formulation.
    lp = Model(milp_optimizer)
    @variable(lp, -1.0 <= x_lp <= 1.0)
    @variable(lp, y_lp)
    construct_univariate_relaxation!(lp, f, x_lp, y_lp, partition, false)

    # univariate function data

    # Collect possible solutions of relaxations.
    sec_verts, tan_verts = PR._collect_vertices(ufd)
    sln_verts = PR.Vertex2d[]
    push!(sln_verts, sec_verts[1])
    append!(sln_verts, tan_verts)
    push!(sln_verts, sec_verts[end])

    tol = 1e-5

    for α in collect((π/20):(π/20):(2*π/5))
        set_objective_function(milp, y - x * tan(α))
        set_objective_function(lp, y_lp - x_lp * tan(α))

        for s in [MOI.MIN_SENSE, MOI.MAX_SENSE]
            # Solve MILP relaxation.
            set_objective_sense(milp, s)
            optimize!(milp)
            @test termination_status(milp) == OPTIMAL

            # Verify that the solution is one of the candidate solutions in `sln_verts`.
            milp_obj = objective_value(milp)
            sln_found = false
            x_sln = value.(x)
            y_sln = value.(y)
            for v in sln_verts
                if isapprox(v[1], x_sln, atol = tol) &&
                   isapprox(v[2], y_sln, atol = tol)
                    sln_found = true
                    break
                end
            end
            @test sln_found == true

            # Verify that the LP relaxation obtained by dropping binary restrictions gives the
            # same objective value as the MILP relaxation.
            relax_integrality(milp)
            optimize!(milp)
            @test termination_status(milp) == OPTIMAL
            lp1_obj = objective_value(milp)
            @test milp_obj ≈ lp1_obj atol = tol

            # Solve LP relaxation formulated in the convex hull form.
            # Verify that the LP relaxation formulated in the convex hull form has the same
            # objective as the above 2 relaxations.
            set_objective_sense(lp, s)
            optimize!(lp)
            @test lp1_obj ≈ objective_value(lp) atol = tol
        end
    end
end
