#! /usr/bin/env python
import numpy as np
import subprocess

data = np.loadtxt('learning-curve.out', skiprows=34)
Ha2eV = 27.211386245988

data[:, 1] *= Ha2eV
data[:, 2] *= Ha2eV
data[:, 9] *= Ha2eV
data[:, 10] *= Ha2eV

E = np.sqrt(np.square(data[:, 1]) + np.square(data[:, 2]))
data[:, 3] = E
data = data[data[:, 3].argsort()[::-1]]
data = np.around(data, 7)

best_energy = f"{data[-1, 1]}  {data[-1, 2]}  {data[-1, 9]}  {data[-1, 10]}"
best_epoch = str(int(data[-1, 0])).zfill(6)

print(f"The picked epoch is: {best_epoch}")
print(best_energy)

subprocess.call(f'for name in *.{best_epoch}.out; do cp ${{name}} ${{name:0:11}}.data; done', shell=True, executable='/bin/bash')
subprocess.call(f'echo epoch {best_epoch} > pickedweight', shell=True, executable='/bin/bash')
subprocess.call(f'echo {best_energy} >> pickedweight', shell=True, executable='/bin/bash')
subprocess.call("echo '    eV/atom             eV/bohr     ' >> pickedweight", shell=True, executable='/bin/bash')
subprocess.call("echo 'Etrain    Etest     Ftrain    Ftest' >> pickedweight", shell=True, executable='/bin/bash')
subprocess.call(f'echo {best_energy} >> pickedweight', shell=True, executable='/bin/bash')
