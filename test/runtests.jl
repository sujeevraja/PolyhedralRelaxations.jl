using PolyhedralRelaxations
import Memento

const PR = PolyhedralRelaxations

PR.logger_config!("debug")

using Test
using JuMP
using GLPK

GLPK.jl_set_preemptive_check(false)
glpk_optimizer =
    JuMP.optimizer_with_attributes(GLPK.Optimizer, "msg_lev" => GLPK.OFF, "tm_lim" => 100.0)

@testset "PolyhedralRelaxations" begin

    include("api_tests.jl")
    include("sampling_tests.jl")
    include("error_tests.jl")
    include("linear_objective_tests.jl")

end
