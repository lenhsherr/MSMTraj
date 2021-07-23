from sys import stdout
from pathlib import Path
from typing import List, Union

# from openff.toolkit.topology import Molecule
# import parmed
import mdtraj as md
from mdtraj.core.trajectory import load
from simtk import unit
from simtk.openmm import app, openmm, Vec3, openmm
from simtk.openmm.openmm import System, Integrator, State, Context
import numpy as np


def prep_pdb(pdb_path: Path, forcefield: app.ForceField, nb_cutoff: unit.Quantity, water_model: str) -> app.Modeller:
        pdb = app.PDBFile(str(pdb_path))

        zero = 0*unit.nanometer
        padding = 8*unit.angstroms
        salt_concentration=150 * unit.millimolar

        dims = pdb.getPositions(asNumpy=True).max(axis=0) - pdb.getPositions(asNumpy=True).min(axis=0)
        dims = [max(x, 2*nb_cutoff) for x in dims]
        pdb.topology.setPeriodicBoxVectors([Vec3(dims[0],zero, zero), Vec3(zero, dims[1], zero), Vec3(zero, zero, dims[2])])

        modeller = app.Modeller(pdb.topology, pdb.positions)
        modeller.addHydrogens(forcefield)
        modeller.addSolvent(forcefield, model=water_model, padding=padding, ionicStrength=salt_concentration)

        return modeller


def minimize_system(simulation: app.Simulation, start_positions: unit.Quantity, steps: int) -> app.Simulation:
        simulation.context.setPositions(start_positions)
        simulation.minimizeEnergy(steps)
        return simulation


def equilibrate_system(simulation: app.Simulation, steps: int) -> app.Simulation:
        
        simulation.reporters.append(app.StateDataReporter(stdout, int(steps/10), step=True,
                potentialEnergy=True, temperature=True, volume=True))
        simulation.step(steps)
        return simulation


def save_system(out_dir: Path, name: str, simulation: app.Simulation,  system: System, integegrator: Integrator) -> None:   

        state = simulation.context.getState(getPositions=True, getVelocities=True, getEnergy=True, getForces=True) 

        with out_dir.joinpath('state.xml').open('wt') as f:
                xml = openmm.XmlSerializer.serialize(state) 
                f.write(xml)   


        with out_dir.joinpath('system.xml').open('wt') as f:
                xml = openmm.XmlSerializer.serialize(system) 
                f.write(xml)   

        with out_dir.joinpath('integrator.xml').open('wt') as f:
                xml = openmm.XmlSerializer.serialize(system) 
                f.write(xml)   

        with out_dir.joinpath(f'{name}_equil.pdb').open('wt') as f:
                app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)


def load_system(work_dir: Path) -> List[Union[State, System, Integrator]]:

        with work_dir.joinpath('state.xml').open('rt') as f:
                state = openmm.XmlSerializer.deserialize(f.read())  

        with work_dir.joinpath('system.xml').open('rt') as f:
                system = openmm.XmlSerializer.deserialize(f.read())  

        with work_dir.joinpath('integrator.xml').open('rt') as f:
                integrator = openmm.XmlSerializer.deserialize(f.read())  

        return [state, system, integrator]


def prep_system(work_dir: Path):

        pdb_path = list(work_dir.glob('*.pdb'))[0]
        pdb_name = str(pdb_path.stem)

        # General parameters
        temp = 300*unit.kelvin
        pressure = 1*unit.atmospheres
        timestep = 0.004*unit.picoseconds
        friction = 1/unit.picoseconds
        nb_cutoff = 1*unit.nanometers
        forcefield = app.ForceField('amber14-all.xml', f'amber14/tip3p.xml')
        integrator = openmm.LangevinIntegrator(temp, friction, timestep)
        minimize_steps = 1000
        equilibration_length = 1*unit.nanoseconds
        equilibration_steps = int(equilibration_length/timestep)


        # Create model (topology and positions)
        model = prep_pdb(pdb_path, forcefield, nb_cutoff, water_model='tip3p')

        # Create system (the math description)
        system = forcefield.createSystem(model.topology, nonbondedMethod=app.PME,
                nonbondedCutoff=nb_cutoff, constraints=app.HBonds, hydrogenMass=4*unit.amu, rigidWater=True)

        # Add barostat
        barostat = openmm.MonteCarloBarostat(pressure, temp) 
        system.addForce(barostat) 

        # Create simulation (simulations do the sampling)
        simulation = app.Simulation(model.topology, system, integrator)

        # Minimize
        simulation = minimize_system(simulation, model.positions, minimize_steps)

        # Equilibrate
        simulation = equilibrate_system(simulation, equilibration_steps)
        
        # Save 
        save_system(work_dir, pdb_name, simulation,  system, integrator)





# def make_seeds(traj: md.Trajectory, work_dir: Path) -> None:

#         for frame in traj:
#                 make_seed(frame)


if __name__=='__main__':

        prep_system(Path('ala3'))

        #         
        # simulation = app.Simulation()
        # simulation.loadState(f'prepped_system/{name}_minimized.xml')


# simulation.reporters.append(app.PDBReporter(f'trajectories/{name}.pdb', 1000))
# simulation.reporters.append(app.StateDataReporter(stdout, 1000, step=True,
#         potentialEnergy=True, temperature=True))
# simulation.step(10000)



