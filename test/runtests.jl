using PolyhedralRelaxations
import Memento
using JuMP

const PR = PolyhedralRelaxations

PR.logger_config!("debug")

using Test
using Cbc
using Ipopt

cbc_optimizer = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)
ipopt_optimizer = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes")

@testset "PolyhedralRelaxations" begin

    include("api_tests.jl")
    # include("error_tests.jl")
    # include("objective_tests.jl")
    # include("tolerance_tests.jl")

end
