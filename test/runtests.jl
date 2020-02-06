println("Testing PolyhedralRelaxations.jl ...")

using PolyhedralRelaxations

const PR = PolyhedralRelaxations 

import Memento

# Suppress warnings during testing.
const TESTLOG = Memento.getlogger(PolyhedralRelaxations)
Memento.setlevel!(TESTLOG, "error")

PR.main()