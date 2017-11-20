import numpy as np
from chemkin_G10 import chemkin

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
    assert(np.all(chemkin.ChemUtil.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
                                                ,np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0]), 1500, None, [False, False]) == [40., 20.]))
    
def test_progress_rate_rj_neg():
    try:
        chemkin.ChemUtil.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
                                                ,np.array([10.0, -10.0]), np.array([2.0, 1.0, 1.0]), 1500, None, [False, False])
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_progress_rate_concs_neg():
    try:
        chemkin.ChemUtil.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
                                                ,np.array([10.0, 10.0]), np.array([-2.0, 1.0, 1.0]), 1500, None, [False, False])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_rate_nu_react_neg():
    try:
        chemkin.ChemUtil.progress_rate(np.array([[-2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
                                                ,np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0]), 1500, None, [False, False])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_reaction_rate():
    assert(np.all(chemkin.ChemUtil.reaction_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[0.0, 1.0], [1.0, 0.0], [0.0, 2.0]])
                                                ,np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0]), 1500, None, [False, False]) == [-80, 0, 20]))


# test parse 
def test_parse_reactionList():
    reactionList = chemkin.ChemUtil.parse( "../static/xml/rxns_short_2.xml", 340, 8.314)[0]
    assert(len(reactionList) == 3)
    assert(reactionList[0].reactStr == "2H2 + O2 =] 2OH + H2")
    assert(reactionList[1].reactStr == "OH + HO2 =] H2O + O2")
    assert(reactionList[2].reactStr == "H2O + O2 =] HO2 + OH")

def test_parse_species():
    species = chemkin.ChemUtil.parse( "../static/xml/rxns_short_2.xml", 340, 8.314)[1]
    assert(species[0] == "H2")
    assert(len(species) == 5)

def test_parse_Not_XML():
    try:
        chemkin.ChemUtil.parse("../static/xml/rxns_short_2.pdf", 340, 8.314)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_parse_T_neg():
    try:
        chemkin.ChemUtil.parse("../static/xml/rxns_short_2.xml", -340, 8.314)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_parse_R_neg():
    try:
        chemkin.ChemUtil.parse("../static/xml/rxns_short_2.xml", 340, -8.314)
    except ValueError as err:
        assert(type(err) == ValueError)