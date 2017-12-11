import requests
import numpy as np
import json

concs = np.array([0.5, 0, 0, 2, 0, 1, 0, 0])
url_local = "http://127.0.0.1:8000/visualization/viz/"
files = {'file': open('../tests/data/xml/rxns_reversible.xml', 'rb')}
r = requests.post(url_local, data = {'T':2500, 'concs':json.dumps(concs.tolist())}, files=files)
# print(r.text)
url_viz = "http://127.0.0.1:8000/visualization/demo/" + r.text + "/"