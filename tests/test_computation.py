import numpy as np
import chemkin_g10.computation as cp

# Test k_const
def test_k_const():
    assert cp.k_const(4.0) == 4.0

def test_k_const_neg():
    try:
        cp.k_const(-1.0)
    except ValueError as err:
        assert(type(err) == ValueError)

# Test k_arr
def test_k_arr():
    assert cp.k_arr(2.0, 3.0, 100.0) == 1.9927962618542914
    
def test_k_arr_A_neg():
    try:
        cp.k_arr(-1.0, 3.0, 100.0)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_k_arr_T_neg():
    try:
        cp.k_arr(1.0, 3.0, -100.0)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_arr_R_neg():
    try:
        cp.k_arr(1.0, 3.0, 100.0, -45.0)
    except ValueError as err:
        assert(type(err) == ValueError)

# Test k_mod_arr
def test_k_mod_arr():
    assert cp.k_mod_arr(2.0, -0.5, 3.0, 100.0) == 0.19927962618542916
    
def test_k_mod_arr_A_neg():
    try:
        cp.k_mod_arr(-1.0, 3.0, 4.0, 100.0)
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_k_mod_arr_T_neg():
    try:
        cp.k_mod_arr(1.0, 3.0, 4.0, -100.0)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_k_mod_arr_R_neg():
    try:
        cp.k_mod_arr(1.0, 3.0, 100.0, -4.9)
    except ValueError as err:
        assert(type(err) == ValueError)

# Test progress_rate
def test_progress_rate():
    assert(np.all(cp.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
                                                ,np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0]), 1500, None, [False, False]) == [40., 20.]))
    
def test_progress_rate_rj_neg():
    try:
        cp.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
                                                ,np.array([10.0, -10.0]), np.array([2.0, 1.0, 1.0]), 1500, None, [False, False])
    except ValueError as err:
        assert(type(err) == ValueError)
        
def test_progress_rate_concs_neg():
    try:
        cp.progress_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
                                                ,np.array([10.0, 10.0]), np.array([-2.0, 1.0, 1.0]), 1500, None, [False, False])
    except ValueError as err:
        assert(type(err) == ValueError)

def test_progress_rate_nu_react_neg():
    try:
        cp.progress_rate(np.array([[-2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]])
                                                ,np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0]), 1500, None, [False, False])
    except ValueError as err:
        assert(type(err) == ValueError)

# Test reaction_rate
def test_reaction_rate():
    assert(np.all(cp.reaction_rate(np.array([[2.0, 1.0], [1.0, 0.0], [0.0, 1.0]]), np.array([[0.0, 1.0], [1.0, 0.0], [0.0, 2.0]])
                                                ,np.array([10.0, 10.0]), np.array([2.0, 1.0, 1.0]), 1500, None, [False, False]) == [-80, 0, 20]))

