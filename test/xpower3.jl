@testset "x^3 test" begin

    uf = PR.UnivariateFunction(
        x->x^3,  # f
        x->3 * (x^2),  # f'
        domain_lb = -1.0,
        domain_ub = 1.0,
        inflection_points = Vector{Real}(collect(-1.0:1.0:1.0)))
    
    model = PR.build_model(uf)

    @test model.x_index == 1
    @test model.y_index == 2

    m = JuMP.Model(glpk_optimizer)
    num_variables = 2 + length(model.delta_1_indices) + length(model.delta_2_indices) + 
        length(model.z_indices)
    JuMP.@variable(m, x[1:num_variables])
    JuMP.set_lower_bound(x[1], -1.0)
    JuMP.set_upper_bound(x[1], 1.0)
    for i in model.delta_1_indices
        JuMP.set_lower_bound(x[i], 0)
        JuMP.set_upper_bound(x[i], 1)
    end 
    for i in model.delta_2_indices
        JuMP.set_lower_bound(x[i], 0)
        JuMP.set_upper_bound(x[i], 1)
    end 
    for i in model.z_indices
        JuMP.set_lower_bound(x[i], 0)
        JuMP.set_upper_bound(x[i], 1)
    end 
    A = model.A 
    b = model.b
    for i in 1:model.num_constraints 
        indices, values = SparseArrays.findnz(A[i, :])
        if i in model.equality_row_indices
            JuMP.@constraint(m, dot(x[indices], values) == b[i])
        else 
            JuMP.@constraint(m, dot(x[indices], values) <= b[i])
        end 
    end 
    JuMP.@objective(m, Min, x[1])
    JuMP.optimize!(m)
    @test JuMP.objective_value(m) == -1.0

    JuMP.@objective(m, Max, x[1])
    JuMP.optimize!(m)
    @test JuMP.objective_value(m) == 1.0
end