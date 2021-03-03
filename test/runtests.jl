using PolyhedralRelaxations
import Memento
using JuMP

const PR = PolyhedralRelaxations

PR.logger_config!("debug")

using Test
using Cbc

cbc_optimizer = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

@testset "PolyhedralRelaxations" begin

    include("api_tests.jl")
    include("error_tests.jl")
    include("objective_tests.jl")
    include("tolerance_tests.jl")

end
