@testset "x^3 test" begin
    formulation_data, function_data = PR.construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))

    @test formulation_data.x_index == 1
    @test formulation_data.y_index == 2

    m = JuMP.Model(glpk_optimizer)
    num_variables = 2 + length(formulation_data.delta_1_indices) +
        length(formulation_data.delta_2_indices) +
        length(formulation_data.z_indices)
    JuMP.@variable(m, x[1:num_variables])
    JuMP.set_lower_bound(x[1], -1.0)
    JuMP.set_upper_bound(x[1], 1.0)
    for i in formulation_data.delta_1_indices
        JuMP.set_lower_bound(x[i], 0)
        JuMP.set_upper_bound(x[i], 1)
    end
    for i in formulation_data.delta_2_indices
        JuMP.set_lower_bound(x[i], 0)
        JuMP.set_upper_bound(x[i], 1)
    end
    for i in formulation_data.z_indices
        JuMP.set_lower_bound(x[i], 0)
        JuMP.set_upper_bound(x[i], 1)
    end
    A = formulation_data.A
    b = formulation_data.b
    for i in 1:formulation_data.num_constraints
        indices, values = SparseArrays.findnz(A[i, :])
        if i in formulation_data.equality_row_indices
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
