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
        @test objective_value(m) <= gopt_value
    end

    @testset "test multilinear MILP relaxation" begin
        PR.silence!()
        num_discretizations = [3, 3, 3, 2, 2, 2, 2, 2, 2, 2]
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
        @test objective_value(m) <= gopt_value
    end
end
