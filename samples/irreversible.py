from chemkin_g10 import chemkin as ck
import numpy as np

T = 2500.0
R = 8.314
concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
rsystem = ck.ReactionSystem(T, R, "../tests/data/db/nasa.sqlite")
rsystem.buildFromXml("../tests/data/xml/rxns_irreversible.xml", concs)
print("Reaction rate: \n", rsystem.getReactionRate(), "\n")
print("System info: \n", rsystem, "\n")


T = 900.0
rsystem = ck.ReactionSystem(T, R, "../tests/data/db/nasa.sqlite")
rsystem.buildFromXml("../tests/data/xml/rxns_irreversible.xml", concs)
print("Reaction rate: \n", rsystem.getReactionRate(), "\n")
print("System info: \n", rsystem, "\n")

