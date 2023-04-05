# Add Hs to a ligand for which we have a PDB structure without hydrogens
# and a SMILES

import sys, os
from rdkit import Chem
from rdkit.Chem import AllChem
import copy
from openeye import oechem # OpenEye Python toolkits
from openeye import oeomega # Omega toolkit
from openeye import oequacpac #Charge toolkit


if (1):
    # Parameterize using OpenFF Toolkit and run a  short simulation (in vaccuum)
    from openforcefield.topology import Molecule
    from openforcefield.typing.engines.smirnoff import ForceField
    istream = oechem.oemolistream('LIG.sdf')
    mol = oechem.OEMol()
    oechem.OEReadMolecule(istream, mol)
    oequacpac.OESetNeutralpHModel(mol)
    oechem.OEAddExplicitHydrogens(mol)
    oechem.OEDetermineConnectivity(mol)
    oechem.OEPerceiveBondOrders(mol)
    ligand = Molecule.from_openeye(mol)
    ligand_positions = ligand.conformers[0]
    lig_ff = ForceField('openff_unconstrained-1.2.1.offxml')
    ligand_topology = ligand.to_topology()
    ligand_system = lig_ff.create_openmm_system(ligand_topology)

    from simtk.openmm.app import *
    from simtk.openmm import *
    from simtk import unit
    import parmed as pmd
    from parmed.openmm import NetCDFReporter
    import numpy as np
    import mdtraj as mdt
    from tempfile import NamedTemporaryFile

    omm_forcefield = app.ForceField('amber14/protein.ff14SB.xml', 'amber14/tip3p.xml')
    pdb = app.PDBFile('protein.pdb')
    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.addSolvent(omm_forcefield,
                        model='tip3p',
                        padding=10.*unit.angstrom,
                        ionicStrength=0.2*unit.molar)
    protein_system = omm_forcefield.createSystem(modeller.topology,
                                                 nonbondedMethod=app.PME,
                                                 rigidWater=False)

    pmd_receptor_struct = pmd.openmm.load_topology(modeller.topology,
                                                   protein_system,
                                                   modeller.positions)
    pmd_ligand_struct = pmd.openmm.load_topology(ligand_topology.to_openmm(),
                                                 ligand_system,
                                                 ligand_positions)

    pmd_complex_struct = pmd_receptor_struct + pmd_ligand_struct

    # Assign periodic box vectors from the solvated receptor structure
    pmd_complex_struct.box_vectors = modeller.topology.getPeriodicBoxVectors()


    def find_clashing_water(pmd_struct, lig_resname, distance):
        """
        Find waters that are sterically clashing with a ligand.
        
        Parameters
        ----------
        pmd_struct : parmed.Structure
            The structure to analyze.
        lig_resname : str
            The up-to-three character residue name.
        distance : float
            The distance cutoff (in nanometers) for clash detection.
            
        Returns
        -------
        water_resnums : Iterable[int]
            The residue numbers of waters that are clashing with the ligand.
            
        """
        with NamedTemporaryFile(suffix='.pdb') as tf:
            app.PDBFile.writeFile(pmd_struct.topology, pmd_struct.positions, open(tf.name, 'w'))
            traj = mdt.load(tf.name)
        top = traj.topology
     #   table,bonds = top.to_dataframe()
     #   print(table.head())
        lig_atom_idxs = top.select(f'resname {lig_resname}')
        lig_res_idx = top.atom(lig_atom_idxs[1]).residue.index
        wat_atom_idxs = top.select('resname HOH and name O')
        wat_res_idxs = [top.atom(i).residue.index for i in wat_atom_idxs]
        potential_contacts = [(lig_res_idx, wat_res_idx) for wat_res_idx in wat_res_idxs]
        contacts = mdt.compute_contacts(traj,
                                        contacts=potential_contacts,
                                        scheme='closest',
                                        ignore_nonprotein=False)


        # Note that this is 0-indexed, while the parmed structure is 
        # 1-indexed, therefore we add 1 before returning
        clash_res_idx = [i[1]+1 for i in contacts[1][(contacts[0] < distance)[0,:]]]
        return clash_res_idx


    clashes = find_clashing_water(pmd_complex_struct, 'LIG', 0.15)

    if len(clashes) != 0:
        clash_residues_str = ','.join([str(i) for i in clashes])
        print(f'Removing ligand-clashing water residues {clash_residues_str}')
        pmd_complex_struct.strip(f':{clash_residues_str}')
    else:
        print("No ligand-water clashes to resolve")


    pmd_complex_struct.save('complex.pdb', overwrite=True)

