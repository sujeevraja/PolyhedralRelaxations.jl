@testset "multilinear relaxations" begin
    @testset "test multilinear LP relaxation" begin
        PR.silence!()
        instance = create_multilinear_model()
        m = instance.model
        x = instance.variables[:x]
        y = instance.variables[:y]
        partitions = instance.partitions
        multilinear_terms = instance.multilinear_terms
        gopt_value = instance.objective
        for term in multilinear_terms
            z = first(term)
            vars = last(term)
            construct_multilinear_relaxation!(m, vars, z, partitions)
        end
        set_optimizer(m, milp_optimizer)
        optimize!(m)
        relaxation_obj = round(objective_value(m); digits = 4)
        gap = round(
            abs(relaxation_obj - gopt_value) / abs(gopt_value) * 100.0;
            digits = 2
        )
        time = round(solve_time(m); digits = 2)
        @info "Relative gap (LP): $(gap) %, time: $(time) sec."
        @test objective_value(m) <= gopt_value
    end

    @testset "test multilinear MILP relaxation" begin
        PR.silence!()
        num_discretizations = [3, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        instance =
            create_multilinear_model(num_discretizations = num_discretizations)
        m = instance.model
        x = instance.variables[:x]
        y = instance.variables[:y]
        partitions = instance.partitions
        multilinear_terms = instance.multilinear_terms
        gopt_value = instance.objective
        for term in multilinear_terms
            z = first(term)
            vars = last(term)
            construct_multilinear_relaxation!(m, vars, z, partitions)
        end
        set_optimizer(m, milp_optimizer)
        optimize!(m)
        relaxation_obj = round(objective_value(m); digits = 4)
        gap = round(
            abs(relaxation_obj - gopt_value) / abs(gopt_value) * 100.0;
            digits = 2
        )
        time = round(solve_time(m); digits = 2)
        @info "Relative gap (MILP): $(gap) %, time: $(time) sec."
        @test objective_value(m) <= gopt_value
    end

    @testset "test multilinear LP relaxation with linking constraints" begin
        PR.silence!()
        instance = create_multilinear_model()
        m = instance.model
        x = instance.variables[:x]
        y = instance.variables[:y]
        partitions = instance.partitions
        multilinear_terms = instance.multilinear_terms
        info = Dict(
            last(term) => construct_multilinear_relaxation!(
                m,
                last(term),
                first(term),
                partitions,
            ) for term in multilinear_terms
        )
        gopt_value = instance.objective
        formulation_info =
            add_multilinear_linking_constraints!(m, info, partitions)
        set_optimizer(m, milp_optimizer)
        optimize!(m)
        relaxation_obj = round(objective_value(m); digits = 4)
        gap = round(
            abs(relaxation_obj - gopt_value) / abs(gopt_value) * 100.0;
            digits = 2
        )
        time = round(solve_time(m); digits = 2)
        @info "Relative gap (LP with linking constraints): $(gap) %, time: $(time) sec."
        @test objective_value(m) <= gopt_value
        @test !isempty(formulation_info.extra[:common_subterm_data])
    end

    @testset "test multilinear MILP relaxation with linking constraints" begin
        PR.silence!()
        num_discretizations = [3, 2, 2, 2, 2, 2, 2, 2, 2, 2]
        instance =
            create_multilinear_model(num_discretizations = num_discretizations)
        m = instance.model
        x = instance.variables[:x]
        y = instance.variables[:y]
        partitions = instance.partitions
        multilinear_terms = instance.multilinear_terms
        info = Dict(
            last(term) => construct_multilinear_relaxation!(
                m,
                last(term),
                first(term),
                partitions,
            ) for term in multilinear_terms
        )
        gopt_value = instance.objective
        formulation_info =
            add_multilinear_linking_constraints!(m, info, partitions)
        set_optimizer(m, milp_optimizer)
        optimize!(m)
        relaxation_obj = round(objective_value(m); digits = 4)
        gap = round(
            abs(relaxation_obj - gopt_value) / abs(gopt_value) * 100.0;
            digits = 2
        )
        time = round(solve_time(m); digits = 2)
        @info "Relative gap (MILP with linking constraints): $(gap) %, time: $(time) sec."
        @test objective_value(m) <= gopt_value
        @test !isempty(formulation_info.extra[:common_subterm_data])
    end
end
