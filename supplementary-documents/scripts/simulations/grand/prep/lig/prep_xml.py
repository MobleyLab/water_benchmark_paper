import sys
from openforcefield.topology import Molecule
from openforcefield.topology import Topology
from openforcefield.typing.engines.smirnoff import ForceField
import parmed as pmd
from parmed.openmm import NetCDFReporter
from tempfile import NamedTemporaryFile
from openeye import oechem # OpenEye Python toolkits
from openeye import oeomega # Omega toolkit
from openeye import oequacpac #Charge toolkit

istream = oechem.oemolistream('KNI.sdf')
mol = oechem.OEMol()
oechem.OEReadMolecule(istream, mol)
oequacpac.OESetNeutralpHModel(mol)
oechem.OEAddExplicitHydrogens(mol)
oechem.OEDetermineConnectivity(mol)
oechem.OEPerceiveBondOrders(mol)
molecule = Molecule.from_openeye(mol)
topology = Topology.from_molecules(molecule)
force_field = ForceField('openff_unconstrained-1.2.1.offxml')
# Make residue template
from openmmforcefields.generators.template_generators import SMIRNOFFTemplateGenerator
template_generator = SMIRNOFFTemplateGenerator(molecules=molecule,forcefield='openff_unconstrained-1.2.1.offxml')
outffxml=template_generator.generate_residue_template(molecule)

# Write
with open('LIG.xml', 'w') as f:
    f.write(outffxml)

