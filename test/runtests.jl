using PolyhedralRelaxations
import Memento

# Suppress warnings during testing.
const TESTLOG = Memento.getlogger(PolyhedralRelaxations)
Memento.setlevel!(TESTLOG, "error")

const PR = PolyhedralRelaxations

using Test
import JuMP
import Ipopt
import GLPK
import MathOptInterface 
using SparseArrays
using LinearAlgebra

const MOI = MathOptInterface 
const MOIU = MOI.Utilities

GLPK.jl_set_preemptive_check(false)
glpk_optimizer = JuMP.with_optimizer(GLPK.Optimizer, tm_lim = 100.0, msg_lev = GLPK.OFF)

@testset "PolyhedralRelaxations" begin

    include("xpower3.jl")

end