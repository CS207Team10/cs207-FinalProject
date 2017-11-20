from chemkin_g10.db import DatabaseOps
import os
path = os.path.dirname(os.path.realpath(__file__)) + "/data/db/"

def test_db_connect_fail():
    filename = path + "foo.sqlite"
    try:
        dbops = DatabaseOps(filename)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_db_get_coeffs():
    filename = path + "nasa.sqlite"
    dbops = DatabaseOps(filename)
    dbops.get_coeffs("H", 700)

def test_db_get_coeffs_none_1():
    filename = path + "nasa.sqlite"
    dbops = DatabaseOps(filename)
    try:
        dbops.get_coeffs("ABC", 1500)
    except ValueError as err:
        assert(type(err) == ValueError)

def test_db_get_coeffs_none_2():
    filename = path + "nasa.sqlite"
    dbops = DatabaseOps(filename)
    try:
        dbops.get_coeffs("O2", 100000)
    except ValueError as err:
        assert(type(err) == ValueError)

