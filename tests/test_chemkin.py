import numpy as np
import chemkin_g10.chemkin as ck

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
    assert(np.allclose(rsystem.getReactionRate(), np.array([6.18148072e+14,  -7.06429290e+14,  -7.24139610e+14,   2.95950938e+13,
                                                            1.35082567e+14,   8.12522426e+14,  -1.06194735e+14,  -5.85845244e+13])))

def test_rsystem_ode():
    concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
    rsystem = ck.ReactionSystem(900, 8.314, path2+'nasa.sqlite')
    rsystem.buildFromXml(path + "rxns_reversible.xml", concs)
    assert(len(rsystem.ode(0.3)) == 50) 

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
