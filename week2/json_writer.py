import math
import json
a = {'freq1' : 11, 'freq2' : 12, 'phase1' : 0, 'phase2' : 0, 'count' : 200}
with open('week2/config.json', 'w') as outfile:
    json.dump(a, outfile)
