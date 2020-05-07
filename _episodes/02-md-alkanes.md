---
title: "MD Simulation of Alkanes"
teaching: 60
exercises: 60
questions:
- "How do I run a simulation?"
- "How do I analyze a simulation?"
objectives:
- "Understand the steps for running a simulation."
- "Be able to perform analysis on  properties such as bond length, angles, and torsions."
- "Understand a PMF calculation."
keypoints:
- "Typical MD simulations consist of simulation initialization, minimization, equilibration, and production."
- "The production simulation gives the data which is analyzed."
- "Results from simulation can be used to calculate the potential of mean force (PMF)."
---

<script type="text/javascript" async
  src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML">
</script>

## What are molecular dynamics simulations?

Now that we have our force field defined, we are ready to do a simulation! But, what does that mean and what kind of simulation will we be doing? We will use a simulation method called **molecular dynamics** (MD). In molecular dynamics simulations, we simulate molecules in time by calculating the forces on atoms and updating their positions based on those forces. The output of a molecular dynamics simulation is a trajectory. The trajectory can be visualized to show the movement of a system through time (like watching a movie).

Force is equal to the negative gradient of potential energy. MD simulations use the potential energy function we discussed previously, and Newton's second law (the force on an object is equal to the object's mass times its acceleration) to calculate positions of atoms. 

$$ \vec{F} = - \nabla U $$

$$ \vec{F} = m \cdot \vec{a} = m \cdot \frac{d\vec{v}}{dt} = m \cdot \frac{d^{2}\vec{r}}{dt^{2}}$$

Looking at these two equations, we se that we we already know the mass for each atom. With our force field defined, we can calculate the force on each atom. The net orce on each atom is the sum of forces from all other atoms in the system:

$$ \vec{F}_i = \sum_{j=1}^{n_{max}}{f_{j}}, (j \neq i) $$

This means that we should be able to calculate positions for each atom based on the force field function, the positions of other atoms in the system, and mass of the molecule. 

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
1. **Analysis** - After you have collected your data from the production run, you must analyze the trajectory to draw conclusions.

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
from sys import stdout

print('Equilibrating...')

simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True, 
    potentialEnergy=True, temperature=True, separator=','))
simulation.context.setVelocitiesToTemperature(150.0*unit.kelvin)
simulation.step(2500)
~~~
{: .language-python}

The first command in this section sets up what information OpenMM will report to us as the simulation runs. We are asking for the potential energy, and temperature every 100 timesteps. By putting `stdout` as the first argument, we are telling the program to just print the information. Note that `stdout` comes from the built in Python module `sys`. If we wanted the information in a file instead, you would put the file name.

The second command sets the initial velocities of the system to a temperature equivalent of 150 K. Then, we integrate for 2,500 steps to allow the system to equilibrate.

### Production

This next block of code is a longer block of simulation called the 'production run'.  We're also added a timer to this code so we can see how long it took our simulation to run.

~~~
import time as time

print('Running Production...')

# Begin timer
tinit=time.time()

# Clear simulation reporters
simulation.reporters.clear()

# Reinitialize simulation reporters. We do this because we want different information printed from the production run than the equilibration run.
# output basic simulation information below every 250000 steps - (which is equal to 2 fs(250,000) = 500,000 fs = 500 ps)
simulation.reporters.append(app.StateDataReporter(stdout, 250000, 
    step=True, time=True, potentialEnergy=True, temperature=True, 
    speed=True, separator=','))

# write out a trajectory (i.e., coordinates vs. time) to a DCD
# file every 100 steps - 0.2 ps
simulation.reporters.append(app.DCDReporter('ethane_sim.dcd', 100))

# run the simulation for 1.0x10^7 steps - 20 ns
simulation.step(10000000)

# End timer
tfinal=time.time()
print('Done!')
print('Time required for simulation:', tfinal-tinit, 'seconds')
~~~
{: .language-python}

> ## Simulation speed - 'ns/day'
> After executing this cell in your notebook, you should see an output which gives you the step number, simulation time, potential energy, temperature, and "speed" for steps in the simulation. The spacing of theses is set in the `simulation.reporters` step where we indicated we wanted information printed every 250,000 timesteps. 
> 
> The "speed" is reported in "ns/day" or "nanoseconds/day". This is a commonly used unit to report how quickly simulations run. It tells you how much simulation time in nanoseconds would pass for 24 hours of computation time. For example, if a simulation is running at 2 ns/day, enough timesteps would be calculated in one day to make 2 nanoseconds of simulation time. If we were using our 2 fs timestep, this would mean that the computer calculated 1,000,000 timesteps over 24 hours.
{: .callout}

> ## Your Turn - Simulation of Butane
>
> Make a copy of the code you wrote to run your ethane simulation. Modify this code to:
> 1. Read in the files `butane.gaff2.xml` and `butane.pdb`
> 1. Carry out a 10 ps MD simulation to bring the butane molecule to an equilibrium temperature of 298 K in which output is printed every 0.5 ps (Leave the minimization portion beforehand unchanged.)
> 1. Carry out a 40 ns MD simulation at 298 K in which output is printed every 1 ns and structures are (still) saved every 0.2 ps into a file called `butane_sim.dcd`.
> 
>> ## Solution
>> ### Simulation Set up
>> ~~~
>> # read in a starting structure for butane and the
>> # corresponding force field file
>> pdb = app.PDBFile('butane.pdb')
>> forcefield = app.ForceField('butane.gaff2.xml')
>>
>> # setup system by taking topology from pdb file;
>> # run gas phase simulation with 2 fs time step (using SHAKE)
>> # at 298.15 K using a Langevin thermostat (integrator) with
>> # coupling constant of 5.0 ps^-1
>> system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.NoCutoff, 
>>                                  constraints=app.HBonds)
>> integrator = mm.LangevinIntegrator(298.15*unit.kelvin, 5.0/unit.picoseconds, 
>>                                    2.0*unit.femtoseconds)
>> integrator.setConstraintTolerance(1e-5)
>> 
>> platform = mm.Platform.getPlatformByName('Reference')
>> simulation = app.Simulation(pdb.topology, system, integrator, platform)
>> simulation.context.setPositions(pdb.positions)
>> ~~~
>> {: .language-python}
>>
>> ### Energy minimization
>> This section should be unchanged from your ethane molecule.
>> ~~~
>> print('Minimizing...')
>> 
>> st = simulation.context.getState(getPositions=True,getEnergy=True)
>> print("Potential energy before minimization is %s" % st.getPotentialEnergy())
>> 
>> simulation.minimizeEnergy(maxIterations=100)
>> 
>> st = simulation.context.getState(getPositions=True,getEnergy=True)
>> print("Potential energy after minimization is %s" % st.getPotentialEnergy())
>> ~~~
>> {: .language-python}
>>
>> ### Equilibration
>> For this section, we double the number of steps in our equilibration since we want a time of 10 ps (5,000 * 2 fs = 10 picoseconds)
>>
>> ~~~
>> print('Equilibrating...')
>> 
>> simulation.reporters.append(app.StateDataReporter(stdout, 250, step=True, potentialEnergy=True, temperature=True, separator=','))
>> simulation.context.setVelocitiesToTemperature(150.0*unit.kelvin)
>> simulation.step(5000)
>> ~~~
>> {: .language-python}
>>
>> #### Production
>> 
>> ~~~
>> print('Running Production...')
>> 
>> tinit=time.time()
>> simulation.reporters.clear()
>> # output basic simulation information below every 500000 steps/1 ns
>> simulation.reporters.append(app.StateDataReporter(stdout, 500000, 
>>     step=True, time=True, potentialEnergy=True, temperature=True, 
>>     speed=True, separator=','))
>> # write out a trajectory (i.e., coordinates vs. time) to a DCD
>> # file every 100 steps/0.2 ps
>> simulation.reporters.append(app.DCDReporter('butane_sim.dcd', 100))
>> 
>> # run the simulation for 2.0x10^7 steps/40 ns
>> simulation.step(20000000)
>> tfinal=time.time()
>> print('Done!')
>> print('Time required for simulation:', tfinal-tinit, 'seconds')
>> ~~~
>> {: .language-python}
> {: .solution}
{: .challenge}

## Analysis

Now that we've performed our computer experiment, it is time to analyze the data we have collected. The main type of data you have collected through this simulation is information on atom positions, or the system trajectory.

As part of our production simulation, we set up a reporter to record atomic positions. The code below shows that code from your previous script, you do not need to execute it.

~~~
simulation.reporters.append(app.DCDReporter('ethane_sim.dcd', 100))
~~~
{: .language-python}

This reporter saved the atomic positions for us every 100 timesteps in a file called `ethane_sim.dcd`. The DCD file format is a binary file (instead of being a text file), so you cannot open it and look at it. However, we will be using certain libraries to analyze and view the file's contents. If you've run your simulation, you should have the file `ethane_sim.dcd` in the same folder as your Jupyter notebook. 

First, we will need to make sure we have a few more Python libraries installed which can help us with analysis. We will use a library called [nglview](http://nglviewer.org/nglview/latest/api.html) to visualize the trajectory, and a library called [MDTraj](http://mdtraj.org/1.9.3/) to analyze the trajectory. Before opening a new notebook for analysis, you may need to install nglview and MDTraj. 

Type the following in your terminal to install nglview and MDTraj:

~~~
$ pip install nglview
$ conda install -c conda-forge mdtraj
~~~
{: .language-bash}

Open a new Jupyter notebook to do the analysis of your ethane trajectory.

First, we will load in our trajectory using MDTraj

~~~
import mdtraj as md

traj = md.load('ethane_sim.dcd', top='ethane.pdb')
~~~
{: .language-python}

The command above reads all of the atomic positions from `ethane_sim.dcd` and keeps track of atom connectivity (topology) which was given in the PDB file. Next, visualize the trajectory using nglview. Nglview has a special function `show_mdtraj` that we can use with our trajectory because it was in a specific format from the MDTraj library.

~~~
import nglview as ngl

visualize = ngl.show_mdtraj(traj)
visualize
~~~
{: .language-python}

This should show you something that looks sort of like a movie of your ethane molecule. These are the atomic positions calculated by OpenMM during the molecular dynamics run. We can now analyze the positions to find out some things about our molecule.

We will use another OpenMM command to pull out our bonds and atoms for analysis

~~~
atoms, bonds = traj.topology.to_dataframe()
atoms
~~~
{: .language-python}

### Analyzing the C-C bond length

Let's look at what C-C bond lengths our ethane molecule had during the simulation. Before we can measure the bond lengths, we have to decide which atoms from our molecule define the bond angle. Below is the output you should have seen above with the `atoms` command (though yours will be styled differently in the Jupyter notebook):

<div>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>serial</th>
      <th>name</th>
      <th>element</th>
      <th>resSeq</th>
      <th>resName</th>
      <th>chainID</th>
      <th>segmentID</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>1</td>
      <td>C1</td>
      <td>C</td>
      <td>1</td>
      <td>ETH</td>
      <td>0</td>
      <td></td>
    </tr>
    <tr>
      <th>1</th>
      <td>2</td>
      <td>H11</td>
      <td>H</td>
      <td>1</td>
      <td>ETH</td>
      <td>0</td>
      <td></td>
    </tr>
    <tr>
      <th>2</th>
      <td>3</td>
      <td>H12</td>
      <td>H</td>
      <td>1</td>
      <td>ETH</td>
      <td>0</td>
      <td></td>
    </tr>
    <tr>
      <th>3</th>
      <td>4</td>
      <td>H13</td>
      <td>H</td>
      <td>1</td>
      <td>ETH</td>
      <td>0</td>
      <td></td>
    </tr>
    <tr>
      <th>4</th>
      <td>5</td>
      <td>C2</td>
      <td>C</td>
      <td>1</td>
      <td>ETH</td>
      <td>0</td>
      <td></td>
    </tr>
    <tr>
      <th>5</th>
      <td>6</td>
      <td>H21</td>
      <td>H</td>
      <td>1</td>
      <td>ETH</td>
      <td>0</td>
      <td></td>
    </tr>
    <tr>
      <th>6</th>
      <td>7</td>
      <td>H22</td>
      <td>H</td>
      <td>1</td>
      <td>ETH</td>
      <td>0</td>
      <td></td>
    </tr>
    <tr>
      <th>7</th>
      <td>8</td>
      <td>H23</td>
      <td>H</td>
      <td>1</td>
      <td>ETH</td>
      <td>0</td>
      <td></td>
    </tr>
  </tbody>
</table>
</div>

We have to pick the atom indices for the C-C bond. An atom's index is the left-most value in the table above. For our torsion, we'll measure `C1-C2` the indices for these are `0` and  `4`. We use the function `compute_distances` in the MDTraj library to measure the distance between these atoms. 

~~~
bond_indices = [0, 4] # atoms to define the bond length
bond_length = md.compute_distances(traj, [bond_indices])
~~~
{: .language-python}

We now have the measurement for this torsion angle in radians for each recorded timestep of the trajectory saved in the array `bond_length`. One way we can examine this data is by plotting it as a histogram using the Python library `matplotlib`.

~~~
import matplotlib.pyplot as plt

bondcounts, binedges, otherstuff = plt.hist(bond_length, bins=120)
plt.title('C-C bond length histogram')
plt.xlabel('Bond length (nm)')
plt.ylabel('Counts')
plt.show()
~~~
{: .language-python}

> ## Exercise - Analyzing the H-C-C-H torsion
>
> A torsion is made up of four atoms which are bonded to each other. Analyze the torsion angle associated with the atoms `H11-C1-C2-H21` for your trajectory. Instead of using the function `compute_distance`, use `compute_dihedrals`. Create a histogram plot of the torsion angles.
>
>> ## Solution
>>
>> First, we need to pick the atom indices of our torsion angle and use the `compute_dihedrals` function to calculate the dihedrals.
>> ~~~
>> phi_indices = [1, 0, 4, 5] # atoms to define the torsion angle
>> phi = md.compute_dihedrals(traj, [phi_indices])
>>
>> print(phi)
>> ~~~
>> {: .language-python}
>> 
>> We now have the measurement for this torsion angle in radians for each recorded timestep of the trajectory. 
>>
>> Next, we can examine this data by plotting it as a histogram using the Python library `matplotlib`.
>> 
>> ~~~
>> import matplotlib.pyplot as plt
>> 
>> phicounts, binedges, otherstuff = plt.hist(phi, bins=90) # create a histogram with 90 bins
>> plt.title('H-C-C-H torsion angle')
>> plt.xlabel(r'$\phi$ (rad)')
>> plt.ylabel('Counts')
>> plt.show()
>> ~~~
>> {: .language-python}
> {: .solution}
{: .challenge}

### Potential of Mean Force Calculation

So far in our analysis, we have looked at the distribution of bond lengths and torsion angles for ethane. However, we can also use our simulations to calculate thermodynamics properties of our system. For example, we can use our calculated distributions along with Boltzmann's constant to calculate the potential of mean force (pmf), or energy change associated with changes in the bond length or torsion angle.

The potential of mean force is defined by the expression

$$ W(x) = -k_{B}T*ln[p(x)] + C $$

where $$p(x)$$ is the probability, or the histogram we calculated previously.

For our torsion angle, we can calculate and plot the PMF:

~~~
kB = 8.31446/1000 # Boltzmann constant in kJ/mol
Temp = 298.15 # simulation temperature
phicounts[phicounts==0] = 0.1 # get rid of any bins with 0 counts/infinite energy
pmf = -kB*Temp*np.log(phicounts) # W(x) = -kT*ln[p(x)] = -kT*ln[n(x)] + C
pmf = pmf - np.min(pmf) # subtract off minimum value so that energies start from 0

bincenters = (binedges[1:] + binedges[:-1])/2 # compute centers of histogram bins

plt.plot(bincenters, pmf)
plt.title('H-C-C-H torsion pmf')
plt.xlabel(r'$\phi$ (rad)')
plt.ylabel('Relative free energy (kJ/mol)')
plt.show()
~~~
{: .language-python}

In the code above, we used the line `pmf = pmf - np.min(pmf)` to subtract the minimum system energy (or $$C$$ from the equation above), giving us the relative free energy.

When we examine the plot, we can see that the PMF is not smooth near the free energy maxima. This is due to finite sampling in our relatively short simulation. This makes sense because the configurations of the molecule with the higher energy would not occur as many times during the simulation. To make this smoother, we could run a longer simulation or use a smoothing function on our data.

#### C-C PMF

Similarly, we can make a plot for the PMF of the C-C bond:

~~~
bondcounts[bondcounts==0] = 0.1
pmf = -kB*Temp*np.log(bondcounts)
pmf = pmf - np.min(pmf)

bincenters = (binedges[1:] + binedges[:-1])/2

pmf_smoothed = sm.nonparametric.lowess(pmf, bincenters, frac=0.05)
pmf_s = pmf_smoothed[:,1] - np.min(pmf_smoothed[:,1])

plt.plot(bincenters, pmf_s)
plt.xlabel('Bond length (nm)')
plt.ylabel('Relative free energy (kJ/mol)')
plt.title('C-C bond length pmf')
plt.show()
~~~
{: .language-python}

Again, we see that our higher energy bond lengths are less smooth. This is because these bond lengths did not occur very much in our simulation (because of the high energy), so our statistics are poor. If we wanted to exclude these areas, we could subset the part of our data which we plot:

~~~
plt.plot(bincenters[pmf_s < 15], pmf_s[pmf_s < 15])
plt.xlabel('Bond length (nm)')
plt.ylabel('Relative free energy (kJ/mol)')
plt.title('C-C bond length pmf')
plt.show()
~~~
{: .language-python}


> ## Your Turn - Analysis of Butane Trajectory
> Create a copy of your ethane analysis code and modify your code to analyze your butane trajectory.
> 1. Read in the files `butane.pdb` and `butane_sim.dcd` and visualize using NGLview.
> 1. Analyze the C-C-C-C torsion angle and compute the PMF.
> 1. Analyze the C-C-C bond angle (use either C1-C2-C3 or C2-C3-C4) and compute the PMF.
> 1. Make only a histogram of one of the C-H bond lengths. Pick any C-H pair. What do you notice about the distribution of this bond length?
>
>> ## Solution
>> ### Read in the MD Trajectory 
>> ~~~
>> traj = md.load('butane_sim.dcd', top='butane.pdb')
>> atoms, bonds = traj.topology.to_dataframe()
>> atoms
>> ~~~
>> {: .language-python}
>>
>> ### Visualize
>> ~~~
>> visualize = ngl.show_mdtraj(traj)
>> visualize
>> ~~~
>> {: .language-python}
>> ### Analyzing the C-C-C-C torsion
>> ~~~
>> phi_indices = [0, 4, 7, 10] # atoms to define the torsion angle
>> phi = md.compute_dihedrals(traj, [phi_indices])
>> 
>> phicounts, binedges, otherstuff = plt.hist(phi, bins=90) # create a histogram with 90 bins
>> plt.title('C-C-C-C torsion angle')
>> plt.xlabel(r'$\phi$ (rad)')
>> plt.ylabel('Counts')
>> plt.show()
>> 
>> print(np.sum(phicounts))
>> ~~~
>> {: .language-python}
>> #### C-C-C-C torsion PMF
>> ~~~
>> B = 8.31446/1000 # Boltzmann constant in kJ/mol
>> Temp = 298.15 # simulation temperature in K
>> phicounts[phicounts==0] = 0.1 # get rid of any bins with 0 counts/infinite energy
>> pmf = -kB*Temp*np.log(phicounts) # W(x) = -kT*ln[p(x)] = -kT*ln[n(x)] + C
>> pmf = pmf - np.min(pmf) # subtract off minimum value so that energies start from 0
>> 
>> bincenters = (binedges[1:] + binedges[:-1])/2 # compute centers of histogram bins
>> 
>> plt.plot(bincenters, pmf)
>> plt.title('C-C-C-C torsion pmf')
>> plt.xlabel(r'$\phi$ (rad)')
>> plt.ylabel('Relative free energy (kJ/mol)')
>> plt.show()
>> ~~~
>> {: .language-python}
>> 
>> ### Analyzing the C-C-C angle
>> ~~~
>> angle_indices = [0, 4, 7] # or could do [4, 7, 10]
>> bondangle = md.compute_angles(traj, [angle_indices])
>> 
>> anglecounts, binedges, otherstuff = plt.hist(bondangle, bins=100)
>> plt.title('C-C-C bond angle')
>> plt.xlabel('Bond angle (rad)')
>> plt.ylabel('Counts')
>> plt.show()
>> ~~~
>> {: .language-python}
>> #### C-C-C angle PMF
>> ~~~
>> anglecounts[anglecounts==0] = 0.1
>> pmf = -kB*Temp*np.log(anglecounts)
>> pmf = pmf - np.min(pmf)
>> 
>> bincenters = (binedges[1:] + binedges[:-1])/2
>> 
>> 
>> plt.plot(bincenters, pmf)
>> plt.xlabel('Bond angle (rad)')
>> plt.ylabel('Relative free energy (kJ/mol)')
>> plt.show()
>> ~~~
>> {: .language-python} 
>> ### Analyzing a C-H bond
>> ~~~
>> bond_indices = [0, 1] # many possibilities!
>> bondlength = md.compute_distances(traj, [bond_indices])
>> 
>> lengthcounts, binedges, otherstuff = plt.hist(bondlength, bins=100)
>> plt.title('C-H bond length')
>> plt.xlabel('Bond length (nm)')
>> plt.ylabel('Counts')
>> plt.show()
>> ~~~
>> {: .language-python}
>> 
>> The C-H bond length does not behave at all like something subject to a harmonic potential.  So what's going on?  Remember that in this simulation we have "frozen" all of the covalent bonds involving H atoms so that we can use a 2 fs time step.  Therefore only the non-H atoms undergo true dynamics; the positions of the H atoms are calculated after each time step using an interative algorithm (SHAKE).  For more information, check out: https://en.wikipedia.org/wiki/Constraint_(computational_chemistry)#The_SHAKE_algorithm
> {: .solution}
{: .challenge}

{% include links.md %}

