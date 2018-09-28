#!/usr/bin/env python3
import os
import json
import numpy as np

os.makedirs('.projects/', exist_ok = True)
with open('adj/nets','r') as f: saved = json.load(f)
saved.remove('darth_vader')

for p in saved:
    with open('adj/' + p + '_vars', 'r') as f: proj_vars = json.load(f)
    s = ''
    tmp = proj_vars[1:] + sorted(proj_vars[0])
    for i in tmp[:-1]:
        s += str(i) + '\n'
    s += str(tmp[-1])
    with open('.projects/' + p + '_vars', 'w') as f: print(s,file=f,end='')

    with open('adj/' + p + '_nt', 'r') as f: proj_nt = json.load(f)
    s = ''
    nt = np.array(proj_nt)
    for i in range(nt.shape[0]-1):
        tmp = np.diag(nt,i)
        for j in tmp:
            s += str(j) + '\n'
    s += str(nt[0,-1])
    with open('.projects/' + p + '_nt', 'w') as f: print(s,file=f,end='')

