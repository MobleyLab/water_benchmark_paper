import sys, os
import numpy as np
import mdtraj as md

t = md.load('complex.gro')
ref = md.load('xtal.pdb')
ind1 = ref.topology.select("protein and mass>1.1")
ind2 = t.topology.select("protein and mass>1.1")
new_t = t.superpose(ref,atom_indices=ind2,ref_atom_indices=ind1)
new_t.save('aligned.pdb')

traj = md.load('traj.nc',top='aligned.pdb')
ind3 = traj.topology.select("protein and mass>1.1")
new_traj = traj.superpose(ref,atom_indices=ind3,ref_atom_indices=ind1)
new_traj.save('aligned.xtc')

os.system('python xtraj.py traj=aligned.xtc top=aligned.pdb first=0 last=12000')
# Uncomment if only water density is needed
#os.system('python xtraj.py traj=aligned.xtc top=aligned.pdb selection="water" first=0 last=12000 fcalc=water_fcalc.mtz icalc=water_icalc.mtz diffuse=water_diffuse.hkl')



