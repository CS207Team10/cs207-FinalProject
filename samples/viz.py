from chemkin_g10 import chemkin as ck
import numpy as np

T = 900
R = 8.314
concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
rsystem = ck.ReactionSystem(T, R, "../tests/data/db/nasa.sqlite")
rsystem.buildFromXml("../tests/data/xml/rxns_reversible.xml", concs)

sim = ck.Simulator(rsystem, 0.05)
sim.solveODE()
sim.visualize()

