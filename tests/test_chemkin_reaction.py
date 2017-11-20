import numpy as np
from chemkin_G10 import chemkin

# Test ReactionSystem (irreversible)
def test_system_build_from_xml():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, 'nasa.sqlite')
    rsystem.buildFromXml("../static/xml/rxns_short_2.xml", concs)
    assert(len(rsystem)==3)

def test_system_progress_rates():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, 'nasa.sqlite')
    rsystem.buildFromXml("../static/xml/rxns_short_2.xml", concs)
    assert(np.allclose(rsystem.getProgressRate(), np.array([2.81117621e+08, 5.00000000e+03, 4.48493847e+06])))

def test_system_reaction_rates():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, 'nasa.sqlite')
    rsystem.buildFromXml("../static/xml/rxns_short_2.xml", concs)
    assert(np.allclose(rsystem.getReactionRate(), np.array([-2.81117621e+08, -2.85597559e+08, 5.66715180e+08, 
                                                            4.47993847e+06, -4.47993847e+06])))

def test_system_update_coeff():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, 'nasa.sqlite')
    rsystem.buildFromXml("../static/xml/rxns_short_2.xml", concs)
    rsystem.reactionList[2].updateCoeff(type="modifiedArrhenius", A=100000000.0, b=0.5, E=50000.0) 
    assert(rsystem.reactionList[2].k == rsystem.reactionList[0].k)

def test_system_rebuild():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, 'nasa.sqlite')
    rsystem.buildFromXml("../static/xml/rxns_short_2.xml", concs)
    rsystem.reactionList[2].updateCoeff(type="modifiedArrhenius", A=100000000.0, b=0.5, E=50000.0) 
    rsystem.buildFromList(rsystem.reactionList, rsystem.species, concs)
    assert(len(rsystem)==3)
    assert(np.allclose(rsystem.getProgressRate(), np.array([2.81117621e+08,5.00000000e+03,7.02794052e+07])))

# Test ReactionSystem (reversible)
def test_rsystem_build_from_xml():
    concs = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, 'nasa.sqlite')
    rsystem.buildFromXml("../static/xml/rxns_reversible.xml", concs)
    assert(len(rsystem)==11)

def test_system_reaction_rates():
    concs = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, 'nasa.sqlite')
    rsystem.buildFromXml("../static/xml/rxns_reversible.xml", concs)
    assert(np.allclose(rsystem.getReactionRate(), np.array([6.18148072e+14,  -7.06429290e+14,  -7.24139610e+14,   2.95950938e+13,
                                                            1.35082567e+14,   8.12522426e+14,  -1.06194735e+14,  -5.85845244e+13])))