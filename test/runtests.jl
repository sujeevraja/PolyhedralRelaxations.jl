using PolyhedralRelaxations
using JuMP

const PR = PolyhedralRelaxations

PR.set_logging_level!(:Debug)

using Test
using HiGHS
using Ipopt

milp_optimizer = JuMP.optimizer_with_attributes(
    HiGHS.Optimizer,
    "presolve" => "on",
    "log_to_console" => false,
)
ipopt_optimizer = JuMP.optimizer_with_attributes(
    Ipopt.Optimizer,
    "print_level" => 0,
    "sb" => "yes",
)

@testset "PolyhedralRelaxations" begin
    include("api_tests.jl")
    include("error_tests.jl")
    include("objective_tests.jl")
    include("tolerance_tests.jl")
    include("multilinear_instance.jl")
    include("multilinear_tests.jl")
    include("partition_refinement_tests.jl")
end
