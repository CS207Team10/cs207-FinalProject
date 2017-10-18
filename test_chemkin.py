import numpy as np
import chemkin

#test k_const
def test_k_const():
    assert chemkin.ChemUtil.k_const(4.0) == 4.0
    
def test_k_const_neg():
    try:
        chemkin.ChemUtil.k_const(-1.0)
    except ValueError as err:
        assert(type(err) == ValueError)

#test k_arr
def test_k_arr():
    assert chemkin.ChemUtil.k_arr(2.0, 3.0, 100.0) == 1.9927962618542914
    
def test_k_arr_A_neg():
    try:
        chemkin.ChemUtil.k_arr(-1.0, 3.0, 100.0)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_k_arr_T_neg():
    try:
        chemkin.ChemUtil.k_arr(1.0, 3.0, -100.0)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_arr_R_neg():
    try:
        chemkin.ChemUtil.k_arr(1.0, 3.0, 100.0, -45.0)
    except ValueError as err:
        assert(type(err) == ValueError)

#test k_mod_arr
def test_k_mod_arr():
    assert chemkin.ChemUtil.k_mod_arr(2.0, -0.5, 3.0, 100.0) == 0.19927962618542916
    
def test_k_mod_arr_A_neg():
    try:
        chemkin.ChemUtil.k_mod_arr(-1.0, 3.0, 4.0, 100.0)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_k_mod_arr_T_neg():
    try:
        chemkin.ChemUtil.k_mod_arr(1.0, 3.0, 4.0, -100.0)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_arr_R_neg():
    try:
        chemkin.ChemUtil.k_mod_arr(1.0, 3.0, 100.0, -4.9)
    except ValueError as err:
        assert(type(err) == ValueError)

#test progress_rate
def test_progress_rate():
    assert(np.all(chemkin.ChemUtil.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0])) == [40., 20.]))
    
def test_progress_rate_rj_neg():
    try:
        chemkin.ChemUtil.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([-10.0, -10.0]), np.array([2.0, 1.0, 1.0]))
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_progress_rate_concs_neg():
    try:
        chemkin.ChemUtil.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([10.0, 10.0]), np.array([-2.0, -1.0, -1.0]))
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_rate_nu_react_neg():
    try:
        chemkin.ChemUtil.progress_rate(np.array([[-2.0, -1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0]))
    except ValueError as err:
        assert(type(err) == ValueError)

# test reaction_rate
def test_reaction_rate():
    assert(np.all(chemkin.ChemUtil.reaction_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([10.0, 10.0]), np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0])) ))

# test parse (To be checked)
def test_parse():
    reactionList = chemkin.ChemUtil.parse( "./test1.xml", 340, 8.314)
    assert(len(reactionList) == 3)
    assert(reactionList[0].reactStr == "2H2 + O2 =] 2OH + H2")
    assert(reactionList[1].reactStr == "OH + HO2 =] H2O + O2")
    assert(reactionList[2].reactStr == "H2O + O2 =] HO2 + OH")

def test_parse_Not_XML():
    try:
        chemkin.ChemUtil.parse( "./test1.pdf", 340, 8.314)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_parse_T_neg():
    try:
        chemkin.ChemUtil.parse( "./test1.xml", -340, 8.314)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_parse_R_neg():
    try:
        chemkin.ChemUtil.parse( "./test1.xml", 340, -8.314)
    except ValueError as err:
        assert(type(err) == ValueError)

# Test ReactionSystem
def test_system_build_from_xml():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, concs)
    rsystem.buildFromXml("test1.xml")
    assert(len(rsystem)==3)

def test_system_progress_rates():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, concs)
    rsystem.buildFromXml("test1.xml")
    assert(np.allclose(rsystem.getProgressRate(), np.array([2.81117621e+08, 5.00000000e+03, 4.48493847e+06])))

def test_system_reaction_rates():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, concs)
    rsystem.buildFromXml("test1.xml")
    assert(np.allclose(rsystem.getReactionRate(), np.array([-2.81117621e+08, -2.85597559e+08, 5.66715180e+08, 
                                                            4.47993847e+06, -4.47993847e+06])))

def test_system_update_coeff():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, concs)
    rsystem.buildFromXml("test1.xml")
    rsystem.reactionList[2].updateCoeff(type="modifiedArrhenius", A=100000000.0, b=0.5, E=50000.0) 
    assert(rsystem.reactionList[2].k == rsystem.reactionList[0].k)

def test_system_rebuild():
    concs = np.array([2.0, 1.0, 0.5, 1.0, 1.0])
    rsystem = chemkin.ReactionSystem(1500, 8.314, concs)
    rsystem.buildFromXml("test1.xml")
    rsystem.reactionList[2].updateCoeff(type="modifiedArrhenius", A=100000000.0, b=0.5, E=50000.0) 
    rsystem.buildFromList(rsystem.reactionList)
    assert(len(rsystem)==3)
    assert(np.allclose(rsystem.getProgressRate(), np.array([2.81117621e+08,5.00000000e+03,7.02794052e+07])))

