import pickle
import os, sys

with open(sys.argv[1], 'rb') as f:
    data = pickle.load(f)

iptm = data["ptm"]
print(iptm)
