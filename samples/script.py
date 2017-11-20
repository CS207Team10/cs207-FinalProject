from chemkin_G10 import chemkin
import numpy as np

T_DEFAULT = 1500
R_DEFAULT = 8.314

if __name__ == '__main__':
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(T_DEFAULT, R_DEFAULT, concs)
    rsystem.buildFromXml("static/xml/rxns_reversible.xml")

    print(rsystem.getProgressRate())
    print(rsystem.getReactionRate())
    print(rsystem)
    
