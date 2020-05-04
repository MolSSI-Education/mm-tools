---
title: "MD Simulation of Alkanes"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

## What are molecular dynamics simulations?

Now that we have our force field defined, we are ready to do a simulation! But, what does that mean and what kind of simulation will we be doing? We will use a simulation method called **molecular dynamics** (MD). In molecular dynamics simulations, we simulate molecules in time by calculating the forces on atoms and updating their positions based on those forces. The output of a molecular dynamics simulation is a trajectory. The trajectory can be visualized to show the movement of a system through time (like watching a movie).

Force is equal to the negative gradient of potential energy. MD simulations use the potential energy function we discussed previously, and Newton's second law to calculate positions of atoms. 

$$ \vec{F} = - \nabla U $$

$$ \vec{F} = m \cdot \vec{a} = m \cdot \frac{d\vec{v}}{dt} = m \cdot \frac{d^{2}\vec{r}}{dt^{2}}$$

Looking at these two equations, we already know the mass for each atom, and with our force field defined, we can calculate the force. This means that we can calculate the acceleration or velocity (and thus a new position) for each atom based on the force field function and mass of the molecule.

$$ \vec{a}_{i}(t_{0}) = \frac{1}{m}\vec{F}_{i}(t_{0}) = -\frac{1}{m}_{i}(\nabla U)_{i} $$

We know that we can calculate positions ($$\vec{x}$$) based on previous positions and acceleration and an amount of time ($$\Delta t$$)

$$\vec{x}_i (t_{0} + \Delta t) = \vec{x}_i (t_{0}) + \vec{v}_{i} (t_{0}) \Delta t + \frac{1}{2} \vec{a}_{i}(t_{0}) \Delta t^{2} $$

This is the general approach to molecular dynamics. You have some atoms which have initial positions and velocities when you start your simulation. You can use these positions to calculate the potential energy and the forces on each atom (because force is the negative gradient of potential energy). With the velocities of the atoms, you can calculate new positions. This process is repeated over and over again, and each iteration is typically called a **timestep**. 

The typical time step size ($$\Delta t$$) represents 1 femtosecond of time and is based on X-H bond vibration frequencies. A timestep that is too small will lead to an inefficient simulation, while a time step that is too large will be inaccurate or cause your simulation to fail. Often, we will constrain hydrogen in a simulation and can use a longer timestep of 2 femtoseconds. 

MD simulations assumes that the [Erogodic hypothesis](https://en.wikipedia.org/wiki/Ergodic_hypothesis) is true, meaning that the time average calculated from the simulation is equal to the ensemble average. For this to be true, the simulation has to sample many conformations. Most simulations have millions to trillions of time steps. 

Usually a simulation protocol follow this general procedure:

1. **Initialization** - Build the system including the topology (connectivity), forcefield parameters, and setting simulation details such as temperature and integrator.
1. **Minimization** - A brief energy minimization to eliminate "bad" interatomic contacts (i.e., ones that would cause high forces) that would result in a numerically unstable simulation.
1. **Equilibration** - A brief MD simulation for the purpose of bringing the system temperature or temperature and volume to the desired equilibrium values.
1. **Production** - A long MD simulation for the purpose of collecting data.

## Running your first simulations

We will now use OpenMM to do a molecular dynamics simulation of the ethane and butane molecules we prepared in the previous lesson. It's important to note at this point that molecular dynamics simulations can be performed using a number of softwares. However, we will be running a simulation with a program called OpenMM. OpenMM has the advantage of being scriptable with Python.

We will first have to make sure we have OpenMM installed. If you are using anaconda, install OpenMM by typing the following command into your terminal or the Anaconda prompt:

~~~
$ conda install -c omnia openmm
~~~
{: .language-bash}


### Simulation Initialization

Once you have OpenMM, we can use it to simulate our molecules. Open a jupyter notebook to run this code. 

Start in your notebook with imports. Here are the python libraries you will need to run simulations with OpenMM.

~~~
from simtk.openmm import app
import simtk.openmm as mm
from simtk import unit
~~~
{: .language-python}

First, we need to read in our structure and our force field. We have to tell the simulation our initial coordinates and the force field we will use. To do this, we use the PDB file we have and the force field file we prepared in the first lesson.

~~~
pdb = app.PDBFile('ethane.pdb')
forcefield = app.ForceField('ethane.gaff2.xml')
~~~
{: .language-python}

We read in the PDB file to get the initial coordinates in the first command. In the second command we give OpenMM the force field file we prepared for ethane before and add it to our `pdb` variable.

Next, we set up the system for our MD simulation. With the following command, we use the `pdb` variable to create a system. The other arguments to the function say that we are not using a cut-off (discussed more later) and that we want to constrain bonds with hydrogens (this allows us to use a larger timestep).

~~~
system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, constraints=app.HBonds)
~~~
{: .language-python}


To run a molecular dynamics simulation we have to integrate Newton's equations of motion. There are different algorithms you might choose to do this, but for this exercise we are going to use a Langevin integrator. Our simulation will be in vacuum at a temperature of 298.15 K. The Langevin integrator is what is called a stochastic integrator. This means that it mimics jostling of air or solvent through random forces. We are using a 5.0 picosecond coupling constant, which is something which controls how often the integrator adds jostling motions.

~~~
integrator = mm.LangevinIntegrator(298.15*unit.kelvin, 5.0/unit.picoseconds, 2.0*unit.femtoseconds)
integrator.setConstraintTolerance(1e-5)
~~~
{: .language-python}


Finally, we initialize the simulation by adding all of the pieces we have prepared:

~~~
platform = mm.Platform.getPlatformByName('Reference')
simulation = app.Simulation(pdb.topology, system, integrator, platform)
simulation.context.setPositions(pdb.positions)
~~~
{: .language-python}

### Energy Minimization

Now, we start calculating energies. First we do an energy minimization. An energy minimization just moves the atoms of the molecule slightly to get to a local minimum in energy. We start in this code block by printing the energy before minimization, doing 100 steps of an energy minimization, then printing the new energy. You should see that the energy decreases:

~~~
print('Minimizing...')

st = simulation.context.getState(getPositions=True,getEnergy=True)
print(F"Potential energy before minimization is {st.getPotentialEnergy()}")

simulation.minimizeEnergy(maxIterations=100)

st = simulation.context.getState(getPositions=True,getEnergy=True)
print(F"Potential energy after minimization is {st.getPotentialEnergy()}")
~~~
{: .language-python}

You can't see it from this code, but the atom positions have changed slightly to cause this change in energy.

### Equilibration

Next, we run an equilibration. The purpose of this equilibration is to get our system to our target temperature and to get the system equilibrated and ready for our production run. 

~~~
print('Equilibrating...')

simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, 
    potentialEnergy=True, temperature=True, separator='\t'))
simulation.context.setVelocitiesToTemperature(150.0*unit.kelvin)
simulation.step(2500)
~~~
{: .language-python}

The first command in this section sets up what information OpenMM will report to us as the simulation runs. We are asking for the potential energy, and temperature every 100 timesteps. By putting `stdout` as the first argument, we are telling the program to just print the information. If we wanted the information in a file instead, you would put the file name.

The second command sets the initial velocities of the system to a temperature equivalent of 150 K. Then, we integrate for 2,500 steps to allow the system to equilibrate.

### Production

~~~
print('Running Production...')

tinit=time.time()
simulation.reporters.clear()
# output basic simulation information below every 250000 steps/500 ps
simulation.reporters.append(app.StateDataReporter(stdout, 250000, 
    step=True, time=True, potentialEnergy=True, temperature=True, 
    speed=True, separator='\t'))
# write out a trajectory (i.e., coordinates vs. time) to a DCD
# file every 100 steps/0.2 ps
simulation.reporters.append(app.DCDReporter('ethane_sim.dcd', 100))

# run the simulation for 1.0x10^7 steps/20 ns
simulation.step(10000000)
tfinal=time.time()
print('Done!')
print('Time required for simulation:', tfinal-tinit, 'seconds')
~~~
{: .language-python}


## Analysis
nglview install from pip. Potential important issue - https://github.com/arose/nglview/issues/718


{% include links.md %}

