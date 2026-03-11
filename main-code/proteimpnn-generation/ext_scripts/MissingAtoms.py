import sys

from modeller import *
from modeller.scripts import complete_pdb

PDBToComplete = sys.argv[1]

env = Environ()
env.libs.topology.read('${LIB}/top_heav.lib')
env.libs.parameters.read('${LIB}/par.lib')

m = complete_pdb(env, PDBToComplete, transfer_res_num="true")
print(m)

# Overwrite the input PDB with the completed model
m.write(file=PDBToComplete)


