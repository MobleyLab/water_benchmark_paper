import sys, os, glob
import numpy as np
import mdtraj
import grand
import mdtraj as md


def analyze(traj,top,resname,atomname,resid,radius,ref,out_traj,out_pdb,cutoff,result):
    traj = md.load(traj,top=top)
    traj_centered = grand.utils.recentre_traj(t=traj,resname=resname,name=atomname,resid=int(resid))
    grand.utils.align_traj(t=traj_centered, output=out_traj)
    grand.utils.write_sphere_traj(radius=radius, ref_atoms=ref, topology=top, trajectory=out_traj, output=out_pdb, initial_frame=True)
    grand.utils.cluster_waters(topology=top,trajectory=out_traj,sphere_radius=radius,ref_atoms=ref,cutoff=cutoff,output=result)


ref = [{'name': 'CA', 'resname': 'ASN', 'resid': '36', 'chain': 0},{'name': 'CA', 'resname': 'THR', 'resid': '169', 'chain': 0}]


trajfile = 'traj.nc'
topfile = 'complex.gro'
resname = 'LIG'
atomname = 'C6'
resid = 8283
radius = 5.5
ref = ref
out_traj = 'center.nc'
out_pdb = 'sphere.pdb'
cutoff = 2.4
result = 'watclusts.pdb'

analyze(trajfile,topfile,resname,atomname,resid,radius,ref,out_traj,out_pdb,cutoff,result)
