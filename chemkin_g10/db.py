import numpy as np
import sqlite3
from pathlib import Path

class DatabaseOps:

    def __init__(self, fileName):
        my_file = Path(fileName)
        if my_file.is_file():
            self.db = sqlite3.connect(str(my_file))
        else:
            raise ValueError("The db file: {} does not exist!".format(str(my_file)))


    def get_coeffs(self, species, T):
        cursor = self.db.cursor()
        a = []
        for s in species:
            # check table LOW
            query = '''SELECT * FROM LOW WHERE SPECIES_NAME="{}"'''.format(s)
            res = cursor.execute(query).fetchall()
            if len(res) != 1:
                raise ValueError("Specie {} not in the database!".format(s))
            if res[0][2] <= T and  T <= res[0][3]: # in the range
                a.append(list(res[0][-7:]))
                continue

            # check table HIGH if not in the LOW's range
            query = '''SELECT * FROM HIGH WHERE SPECIES_NAME="{}"'''.format(s)
            res = cursor.execute(query).fetchall()
            # if len(res) != 1:
            #     raise ValueError("Specie {} not in the database!".format(s))
            if res[0][2] <= T and  T <= res[0][3]:
                a.append(list(res[0][-7:]))
                continue

            # The temperature T is beyond the range
            raise ValueError("The temperate is beyond the range for species {}!".format(s))
        return np.array(a)
        