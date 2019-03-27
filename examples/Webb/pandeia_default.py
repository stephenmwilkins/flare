 

try:
    from pandeia.engine.calc_utils import build_default_calc
except:
    print('failed')

from pandeia.engine.calc_utils import build_default_calc
from pandeia.engine.perform_calculation import perform_calculation


import json
with open('default_input.json', 'r') as inf:
    config = json.loads(inf.read())
   
print(config.keys())
   
report = perform_calculation(config)
print(report['scalar'])