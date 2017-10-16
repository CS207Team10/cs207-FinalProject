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
   