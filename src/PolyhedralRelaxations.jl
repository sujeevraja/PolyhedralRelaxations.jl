module PolyhedralRelaxations

    import Memento

    # Create our module level logger (this will get precompiled)
    const _LOGGER = Memento.getlogger(@__MODULE__)

    # Register the module level logger at runtime so that folks can access the logger via `getlogger(PolyhedralRelaxations)`
    # NOTE: If this line is not included then the precompiled `PolyhedralRelaxations.LOGGER` won't be registered at runtime.
    __init__() = Memento.register(_LOGGER)

    "Suppresses information and warning messages output by GasModels, for fine grained control use the Memento package"
    function silence()
        Memento.info(_LOGGER, "Suppressing information and warning messages for the rest of this session.  Use the Memento package for more fine-grained control of logging.")
        Memento.setlevel!(Memento.getlogger(PolyhedralRelaxations), "error")
    end

    include("main.jl")
    include("typedefs.jl")

end 
