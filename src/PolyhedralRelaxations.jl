module PolyhedralRelaxations

using DataStructures
using JuMP

import ForwardDiff

import Logging 
import LoggingExtras

# Setup Logging
include("logging.jl")
function __init__()
    global _DEFAULT_LOGGER = Logging.current_logger()
    global _LOGGER = Logging.ConsoleLogger(; meta_formatter = PolyhedralRelaxations._pr_metafmt)

    Logging.global_logger(_LOGGER)
end


const EPS = 1e-6
const INF = 1e12

include("common.jl")
include("univariate_lp_relaxation.jl")
include("univariate_milp_relaxation.jl")
include("bilinear_relaxation.jl")
include("api.jl")

end
