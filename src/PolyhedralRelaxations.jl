module PolyhedralRelaxations

using DataStructures
using JuMP

import ForwardDiff
import Memento

# Create our module level logger (this will get precompiled)
const _LOGGER = Memento.getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via
# getlogger(PolyhedralRelaxations)`.
# NOTE: If this line is not included then the precompiled `PolyhedralRelaxations.LOGGER` won't
# be registered at runtime.
__init__() = Memento.register(_LOGGER)

"""
Suppresses information and warning messages output by PolyhedralRelaxations. For fine grained
control use the Memento package.
"""
function silence()
    Memento.info(
        _LOGGER,
        string(
            "Suppressing information and warning messages for the rest ",
            "of this session.  Use the Memento package for more fine-grained control of logging.",
        ),
    )
    Memento.setlevel!(Memento.getlogger(PolyhedralRelaxations), "error")
end

"Allows users to set the logging level without adding Memento."
function logger_config!(level)
    Memento.config!(Memento.getlogger("PolyhedralRelaxations"), level)
end

const EPS = 1e-6
const INF = 1e12

include("common.jl")
include("univariate_lp_relaxation.jl")
include("univariate_milp_relaxation.jl")
include("bilinear_relaxation.jl")
include("api.jl")

end
