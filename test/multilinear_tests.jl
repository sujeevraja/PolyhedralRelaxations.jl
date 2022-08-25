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
        relaxation_obj = round(objective_value(m); digits=4)
        println(" instance: m_10_3_0_100_1 - Global Optimum: $(gopt_value)")
        println(" instance: m_10_3_0_100_1 - Root LP Relaxation: $(relaxation_obj)")
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
        relaxation_obj = round(objective_value(m); digits=4)
        println(" instance: m_10_3_0_100_1 - MILP Relaxation: $(relaxation_obj)")
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
        linking_info = []
        gopt_value = instance.objective
        for term in multilinear_terms
            z = first(term)
            vars = last(term)
            formulation_info = construct_multilinear_relaxation!(m, vars, z, partitions)
            push!(linking_info, formulation_info)
        end
        set_optimizer(m, milp_optimizer)
        optimize!(m)
        relaxation_obj = round(objective_value(m); digits=4)
        println(" instance: m_10_3_0_100_1 - Global Optimum: $(gopt_value)")
        println(" instance: m_10_3_0_100_1 - Root LP Relaxation: $(relaxation_obj)")
        @test objective_value(m) <= gopt_value
    end
end
