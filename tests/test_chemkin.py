import numpy as np
import chemkin_g10.chemkin as ck
from chemkin_g10 import simulator as sim
import matplotlib.pyplot as plt
import matplotlib

import os
path = os.path.dirname(os.path.realpath(__file__)) +  "/data/xml/"
path2 = os.path.dirname(os.path.realpath(__file__)) +  "/data/db/"

# Test ReactionSystem (irreversible)
def test_system_build_from_xml():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = ck.ReactionSystem(1500, 8.314, path2+'nasa.sqlite')
    rsystem.buildFromXml(path + "rxns_short_2.xml", concs)
    assert(len(rsystem)==3)

def test_system_progress_rates():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = ck.ReactionSystem(1500, 8.314, path2+'nasa.sqlite')
    rsystem.buildFromXml(path + "rxns_short_2.xml", concs)
    assert(np.allclose(rsystem.getProgressRate(), np.array([2.81117621e+08, 5.00000000e+03, 4.48493847e+06])))

def test_system_reaction_rates():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = ck.ReactionSystem(1500, 8.314, path2+'nasa.sqlite')
    rsystem.buildFromXml(path + "rxns_short_2.xml", concs)
    assert(np.allclose(rsystem.getReactionRate(), np.array([-2.81117621e+08, -2.85597559e+08, 5.66715180e+08, 
                                                            4.47993847e+06, -4.47993847e+06])))

def test_system_update_coeff():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = ck.ReactionSystem(1500, 8.314, path2+'nasa.sqlite')
    rsystem.buildFromXml(path + "rxns_short_2.xml", concs)
    rsystem.reactionList[2].updateCoeff(type="modifiedArrhenius", A=100000000.0, b=0.5, E=50000.0) 
    assert(rsystem.reactionList[2].k == rsystem.reactionList[0].k)

def test_system_update_reaction():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = ck.ReactionSystem(1500, 8.314, path2+'nasa.sqlite')
    rsystem.buildFromXml(path + "rxns_short_2.xml", concs)
    rsystem.reactionList[2].updateReaction(reversible='yes') 
    assert(rsystem.reactionList[2].reactMeta['reversible'] == 'yes')

def test_system_rebuild():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = ck.ReactionSystem(1500, 8.314, path2+'nasa.sqlite')
    rsystem.buildFromXml(path + "rxns_short_2.xml", concs)
    rsystem.reactionList[2].updateCoeff(type="modifiedArrhenius", A=100000000.0, b=0.5, E=50000.0) 
    rsystem.buildFromList(rsystem.reactionList, rsystem.species, concs)
    print(rsystem)
    assert(len(rsystem)==3)
    assert(np.allclose(rsystem.getProgressRate(), np.array([2.81117621e+08,5.00000000e+03,7.02794052e+07])))


# Test ReactionSystem (reversible)
def test_rsystem_build_from_xml():
    concs = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    rsystem = ck.ReactionSystem(1500, 8.314, path2+'nasa.sqlite')
    rsystem.buildFromXml(path + "rxns_reversible.xml", concs)
    assert(len(rsystem)==11)

def test_rsystem_build_from_xml_fail():
    try:
        concs = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
        rsystem = ck.ReactionSystem(1500, 8.314, path2+'nasa.sqlite')
        rsystem.buildFromXml(path + "rxns_reversible.xml", concs)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_rsystem_reaction_rates():
    concs = np.array([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])
    rsystem = ck.ReactionSystem(1500, 8.314, path2+'nasa.sqlite')
    rsystem.buildFromXml(path + "rxns_reversible.xml", concs)
    assert(np.allclose(rsystem.getReactionRate(), np.array([  6.22261584e+14, -7.10493349e+14, -7.28739230e+14, 2.97882825e+13, 
                                                              1.35132846e+14,  8.16829127e+14, -1.06193909e+14, -5.85853515e+13])))


# Test parse 
def test_parse_reactionList():
    reactionList = ck.ReactionSystem.parse(path + "rxns_short_2.xml", 340, 8.314)[0]
    assert(len(reactionList) == 3)
    assert(reactionList[0].reactStr == "2H2 + O2 =] 2OH + H2")
    assert(reactionList[1].reactStr == "OH + HO2 =] H2O + O2")
    assert(reactionList[2].reactStr == "H2O + O2 =] HO2 + OH")

def test_parse_species():
    species = ck.ReactionSystem.parse(path + "rxns_short_2.xml", 340, 8.314)[1]
    assert(species[0] == "H2")
    assert(len(species) == 5)

def test_parse_Not_XML():
    try:
        ck.ReactionSystem.parse(path + "rxns_short_2.pdf", 340, 8.314)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_parse_T_neg():
    try:
        ck.ReactionSystem.parse(path + "rxns_short_2.xml", -340, 8.314)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_parse_R_neg():
    try:
        ck.ReactionSystem.parse(path + "rxns_short_2.xml", 340, -8.314)
    except ValueError as err:
        assert(type(err) == ValueError)

# Test Simulator
def test_solve_ODE():
    T = 900
    R = 8.314
    concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
    rsystem = ck.ReactionSystem(T, R, "tests/data/db/nasa.sqlite")
    rsystem.buildFromXml("tests/data/xml/rxns_reversible.xml", concs)
    s = sim.Simulator(rsystem, 0.05)
    try:
        s.solveODE()
    except ValueError as err:
        assert(type(err) == ValueError)

# Travis doesn't have GUI, so tkinter wont pass tests here
# def test_invalid_plot():
#     T = 900
#     R = 8.314
#     concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
#     rsystem = ck.ReactionSystem(T, R, "tests/data/db/nasa.sqlite")
#     rsystem.buildFromXml("tests/data/xml/rxns_reversible.xml", concs)
#     s = sim.Simulator(rsystem, 0.05)
#     # sim.solveODE()
#     try:
#         s.plot_specie_all()
#         plt.close('all')
#         s.plot_specie(4)
#         plt.close('all')
#         s.plot_reaction_all()
#         plt.close('all')
#     except ValueError as err:
#         assert(type(err) == ValueError)

def test_equilibrium_1():
    T = 900
    R = 8.314
    concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
    rsystem = ck.ReactionSystem(T, R, "tests/data/db/nasa.sqlite")
    rsystem.buildFromXml("tests/data/xml/rxns_reversible.xml", concs)
    s = sim.Simulator(rsystem, 0.1)
    s.solveODE()
    assert(s.check_equilibrium(5, 0.99) == True)




