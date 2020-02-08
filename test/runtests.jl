using Memento
using PolyhedralRelaxations

# Suppress warnings during testing.
# const TESTLOG = Memento.getlogger(PolyhedralRelaxations)
# Memento.setlevel!(TESTLOG, "error")

const PR = PolyhedralRelaxations

info(PR._LOGGER, "Testing PolyhedralRelaxations.jl ...")
PR.main()