@testset "linear objective function tests" begin
    PR.logger_config!("error")
    f, partition = x -> x^3, collect(-1.0:0.05:1.0)

    # create LP relaxation of the MILP relaxation of the univariate function.
    milp_relaxation, milp_function_data = construct_milp_relaxation(f, partition)
    milp = Model(cbc_optimizer)
    lb, ub = get_variable_bounds(milp_relaxation)
    num_variables = get_num_variables(milp_relaxation)
    @variable(milp, lb[i] <= x[i = 1:num_variables] <= ub[i])
    A, b = get_eq_constraint_matrices(milp_relaxation)
    @constraint(milp, A * x .== b)
    A, b = get_leq_constraint_matrices(milp_relaxation)
    @constraint(milp, A * x .<= b)
    set_objective_coefficient(milp, x[milp_relaxation.y_index], 1.0)
    binary_indices, _ = findnz(milp_relaxation.binary)

    # create LP relaxation of the univariate function using the convex hull formulation.
    lp_relaxation, lp_function_data = construct_lp_relaxation(f, partition)
    lp = Model(cbc_optimizer)
    lb, ub = get_variable_bounds(lp_relaxation)
    num_variables = get_num_variables(lp_relaxation)
    @variable(lp, lb[i] <= y[i = 1:num_variables] <= ub[i])
    A, b = get_eq_constraint_matrices(lp_relaxation)
    @constraint(lp, A * y .== b)
    set_objective_coefficient(lp, y[lp_relaxation.y_index], 1.0)

    # Collect possible solutions of relaxations.
    sec_verts, tan_verts = PR.collect_vertices(milp_function_data)
    sln_verts = PR.Vertex[]
    push!(sln_verts, sec_verts[1])
    append!(sln_verts, tan_verts)
    push!(sln_verts, sec_verts[end])

    tol = 1e-5

    for i = 1:10
        α = rand() * (π / 2)
        set_objective_coefficient(milp, x[milp_relaxation.x_index], -tan(α))
        set_objective_coefficient(lp, y[lp_relaxation.x_index], -tan(α))

        for s in [MOI.MIN_SENSE, MOI.MAX_SENSE]
            # Solve MILP relaxation.
            set_objective_sense(milp, s)
            for k in binary_indices
                set_binary(x[k])
            end
            optimize!(milp)
            @test termination_status(milp) == MOI.OPTIMAL

            # Verify that the solution is one of the candidate solutions in `sln_verts`.
            milp_obj = objective_value(milp)
            sln_found = false
            x_sln = value.(x[milp_relaxation.x_index])
            y_sln = value.(x[milp_relaxation.y_index])
            for v in sln_verts
                if isapprox(v[1], x_sln, atol = tol) && isapprox(v[2], y_sln, atol = tol)
                    sln_found = true
                    break
                end
            end
            @test sln_found == true

            # Verify that the LP relaxation obtained by dropping binary restrictions gives the
            # same objective value as the MILP relaxation.
            for k in binary_indices
                unset_binary(x[k])
            end
            optimize!(milp)
            @test termination_status(milp) == MOI.OPTIMAL
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
