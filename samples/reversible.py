from chemkin_g10 import chemkin as ck
import numpy as np

T = 1500
R = 8.314
concs = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
rsystem = ck.ReactionSystem(T, R, "../tests/data/db/nasa.sqlite")
rsystem.buildFromXml("../tests/data/xml/rxns_reversible.xml", concs)
print("Progress rate: \n", rsystem.getProgressRate(), "\n")
print("Reaction rate: \n", rsystem.getReactionRate(), "\n")
# print("System info: \b", rsystem, "\n")

    
