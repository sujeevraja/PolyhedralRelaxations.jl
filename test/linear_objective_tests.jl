@testset "linear objective function tests" begin
    PR.logger_config!("error")
    milp_relaxation, function_data = construct_milp_relaxation(x -> x^3, collect(-1.0:0.05:1.0))
    sec_verts, tan_verts = collect_vertices(function_data)
    for α ∈ [0:0.01:π]
        println(α)
    end
end
