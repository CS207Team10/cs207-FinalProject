import package.chemkin_2 as ch
import numpy as np

T_DEFAULT = 1500
R_DEFAULT = 8.314

if __name__ == '__main__':
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = ch.ReactionSystem(T_DEFAULT, R_DEFAULT, concs)
    rsystem.buildFromXml("input files/rxns_reversible.xml")

    print(rsystem.getProgressRate())
    print(rsystem.getReactionRate())
    print(rsystem)
    