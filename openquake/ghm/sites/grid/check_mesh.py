import re
import numpy as np
import matplotlib.pyplot as plt

fname = 'out.csv'
dat = []
for line in open(fname, 'r'):
    aa = re.split('\,', line)
    dat.append([float(aa[0]), float(aa[1])])
print('done')
dat = np.array(dat)

fname = 'poly.txt'
pol = []
for line in open(fname, 'r'):
    aa = re.split('\,', line)
    pol.append([float(aa[0]), float(aa[1])])
print('done')
pol = np.array(pol)

plt.plot(dat[:, 0], dat[:, 1], '.')
plt.plot(pol[:, 0], pol[:, 1], '-')
plt.show()
