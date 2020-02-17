using PolyhedralRelaxations
import Memento

const PR = PolyhedralRelaxations

PR.logger_config!("debug")

using Test
using JuMP
using GLPK
using SparseArrays


GLPK.jl_set_preemptive_check(false)
glpk_optimizer = JuMP.optimizer_with_attributes(
    GLPK.Optimizer, "msg_lev" => GLPK.OFF, "tm_lim" => 100.0)

@testset "PolyhedralRelaxations" begin

    include("xpower3.jl")

end
