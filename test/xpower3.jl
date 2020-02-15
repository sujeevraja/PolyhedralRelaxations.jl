@testset "x^3 test" begin
    formulation_data, function_data = PR.construct_milp_relaxation(x -> x^3, collect(-1.0:1.0:1.0))
    @test formulation_data.x_index == 1
    @test formulation_data.y_index == 2

    lb, ub = formulation_data.lower_bounds, formulation_data.upper_bounds
    @test lb[formulation_data.x_index] == function_data.partition[1]
    @test ub[formulation_data.x_index] == function_data.partition[end]
    @test lb[formulation_data.y_index] == -Inf
    @test ub[formulation_data.y_index] == Inf

    num_variables = PR.get_num_variables(formulation_data)
    @test num_variables == 8

    m = JuMP.Model(glpk_optimizer)
    JuMP.@variable(m, lb[i] <= x[i=1:num_variables] <= ub[i])

    A = formulation_data.A
    b = formulation_data.b
    eq_indices, _ = SparseArrays.findnz(formulation_data.equality_row_indices)
    @test eq_indices == [formulation_data.x_index,formulation_data.y_index]

    for i in 1:formulation_data.num_constraints
        indices, values = SparseArrays.findnz(A[i, :])
        if i in eq_indices
            JuMP.@constraint(m, dot(x[indices], values) == b[i])
        else
            JuMP.@constraint(m, dot(x[indices], values) <= b[i])
        end
    end
    JuMP.@objective(m, Min, x[formulation_data.x_index])
    JuMP.optimize!(m)
    @test JuMP.objective_value(m) == -1.0

    JuMP.@objective(m, Max, x[formulation_data.x_index])
    JuMP.optimize!(m)
    @test JuMP.objective_value(m) == 1.0
end
