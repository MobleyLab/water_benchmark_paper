import sys, os
from blues.moves import WaterTranslationMove, MoveEngine
from blues.simulation import *
import json
from blues.settings import Settings

# Parse a YAML configuration, return as Dict
opt = Settings('example.yaml').asDict()
structure = opt['Structure']
# If a restart file is available to continue from a previous run
#restart = parmed.amber.Rst7('water.rst7')
#structure.positions = restart.positions
#structure.box = restart.box
print(json.dumps(opt, sort_keys=True, indent=2, skipkeys=True, default=str))

# Select move type
# atom selected as the center of a spherical region 
water = WaterTranslationMove(structure, water_name='HOH', protein_selection='(resname == "LIG") and (name == "C6")', radius=1*unit.nanometers)

# Initialize object that selects movestep
water_mover = MoveEngine(water)

 #Generate the openmm.Systems outside SimulationFactory to allow modifications
systems = SystemFactory(structure, water.atom_indices, opt['system'])

# Restrain atoms in the MD and alchemical system
systems.md = systems.restrain_positions(structure, systems.md, **opt['restraints'])
systems.alch = systems.restrain_positions(structure, systems.alch, **opt['restraints'])

# Generate the OpenMM Simulations
simulations = SimulationFactory(systems, water_mover, opt['simulation'], opt['md_reporters'], opt['ncmc_reporters'])

# Energy minimize system - Uncomment for  minimization runs
#simulations.md.minimizeEnergy(maxIterations=0) #minimize until convergence is reached
#state = simulations.md.context.getState(getPositions=True, getEnergy=True)
#print('Minimized energy = {}'.format(state.getPotentialEnergy().in_units_of(unit.kilocalorie_per_mole)))

# MD simulation - Uncomment for equilibration and MD production runs only
#simulations.md.step(opt['simulation']['nstepsMD'])

blues = BLUESSimulation(simulations, opt['simulation'])
blues.run()

