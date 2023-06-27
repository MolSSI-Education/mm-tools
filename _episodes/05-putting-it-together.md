---
title: "Putting It All Together"
teaching: 20
exercises: 60
questions:
- "How do I implement a new feature in someone else's code and then share my new feature?"
objectives:
- "Fork an existing code repository into your area on GitHub."
- "Clone a repository from GitHub onto your local computer."
- "Implement a new feature in the code."
- "Share your new code on GitHub."
- "Use a pull request to ask the original code developer to include your new feature."
keypoints:
- "Large coding projects often have hundreds, if not thousands of developers working togehter on code.  Using version control can help many people implement new features at once."
- "Use online documentation from coding projects to learn about existing functions and features when you are trying to change existing code."
- "When you write new code, always document your functions so that other people will be able to understand and build on your work."
---
In this lesson, you will use all of the skills you have used in the previous lessons to make changes to an existing code and share those changes with others.  You will fork an existing repository to your GitHub and then clone that repository to your own computer.  The repository contains an OpenMM script.  You will implement a new feature into the script, commit your changes, and the create a PR to commit your changes to the original repository.  

## Fork the original GitHub repository to your personal account
Log into your GitHub account, visit the original repository located at https://github.com/MolSSI-Education/mm-tools-synthesis and ​click the “Fork” button on the upper right​.

This will make a copy of the original repository under your personal account. It contains all of the code in the original repository; the only difference is that it knows it’s derived from the original one.
Your ​forked repository​ is located at ​https://github.com/[your_username]/mm-tools-synthesis.

[​Note: ​You will be able to make changes to your forked repository, but not the original repository. This is a key aspect of how GitHub manages multiple people contributing features to the same code.]

## Clone the repository to your computer

Make sure you’re in the ​mm-tools folder on your Desktop. Then run the command: git clone https://github.com/[your_username]/mm-tools-synthesis.git

This will create a folder called ​mm-tools-synthesis​ in the folder where you ran the command. Enter the folder​ and you will see several files: ​trpcage.pdb, run_openmm.py, README and
Trp-Cage_AnalysisAndVisualization.ipynb

## Understand the script, run the simulation, and visualize the trajectory

The file trpcage.pdb​ contains a solvated PDB structure of the Trp-cage miniprotein. It’s a bit small for a protein simulation, but designed so that we can get results quickly.

Open the run_openmm.py sript in your text editor.  Open a new jupyter notebook and paste the contents of the script into the cells.  Organize your code like this:
Cell 1: Import statements
Cell 2: The lines pdb = ... through platform = ... (the simulation setup) Cell 3: The rest of the code (running the simulation)

> ## Running from the Command Line
> You could also run the code directly from the command line.  If you choose this option, you will make later edits to the code in your text editor and re-run it from the > command line every time to make a change.
{: .callout}

Run all of the cells in your jupyter notebook.  When the script is running, it will print information to the screen about the simulation time, simulation speed, potential energy, temperature, density, and a progress indicator.

After the script is finished running, you may view the trajectory using the NGLview library in the included jupyter notebook.

## Implement a new feature in the code

The cool thing about OpenMM is that your simulation is inside of the Python script, so you have a lot of flexibility to customize it. You can get creative here and think about different possibilities. These are two  examples of interesting things we can do.

> ## Exercise Option A - Pull apart the protein by applying a force
> Because Trp-cage is a stable protein in solution, running a molecular dynamics simulation will just cause the protein to tumble and vibrate around – perhaps not very interesting. But suppose we want to know how large of a force is needed to pull the protein apart? Imagine applying a spring force between the alpha carbon atoms of the first and last amino acid. Experimentalists do similar things in the lab using​ optical tweezers,​ but they cannot apply precise forces to pairs of atoms; instead they chemically modify the target amino acids, chemically link polystyrene beads to them, and use lasers to apply the forces to the beads. In our simulations we can apply the forces on the atoms themselves.
>
> The alpha carbon (Cα) is the carbon atom on each amino acid that the side chain is bonded to. Open ​trpcage.pdb​ in a text editor.​ Look at the third column containing atom names - the line with CA in the third column is the alpha carbon.  The sixth column is the ​residue ID​ which labels the amino acid or water molecule number. Residue IDs start from 1.
>
> Using NGLview we can visualize the first and last Cα atoms of the protein by adding a new “representation” using the code:
> ~~~
> view.add_representation(repr_type='surface', selection='X.CA or Y.CA', surfaceType='vws')
> ~~~
> {: .language-python}
> where ​X​ and ​Y​ stand for the first and last residue IDs of the protein.  You will need to think about how many amino acids this protein has. In the text editor, scroll down until the residue ID becomes HOH, which is the residue name of water.  The last residue that is not labeled HOH is the last amino acid in the protein.
>
> You can also instantaneously measure the distance between these two atoms using the code:
> ~~~
> view.add_distance(atom_pair=[['X.CA', 'Y.CA']], label_color='black')
> ~~~
> {: .language-python}
>
> There’s a variety of documentation you can use to implement the spring force between these particles. The class you need is called ​CustomBondForce.​ This class allows you to implement any force between selected pairs of atoms, it doesn’t have to represent a covalent bond. To learn about how to use this class, ​visit the Python API (application programing interface) documentation​ located at ​http://docs.openmm.org/latest/api-python/, ​click on Forces, and then CustomBondForce​.
>
> Write the following code after the system object is created, but before the simulation object.
>
> #### 1. Create the CustomBondForce
> Use a command like
> ~~~
> custom_bond = mm.CustomBondForce(“...”)​
> ~~~
> {: .language-python}
> where ... is the energy expression that you want to use for the particles.  For your initial expression use a harmonic spring force where the energy expression is ​
> ~~~
> 0.5*k*(r-r0)^2​
> ~~~
> {: .language-python}
> where k is the spring constant, and r0 is the equilibrium distance. Custom forces are really cool because the user can specify any functional form they want, and OpenMM will actually generate the code for the simulation at runtime.
>
> #### 2. Set your global parameters
> The variables k and r0 are called global parameters. You need to set the force constant using a line such as ​
> ~~~
> custom_bond.addGlobalParameter(“k”, 100)​
> ~~~
> {: .language-python}
> where 100 is value of the parameter in OpenMM’s unit system (kJ mol​$^{-1}$ nm​$^{-1}$ for bond force constants). ​Add a line that sets the r0 parameter to 2 nanometers.​ [OpenMM uses nanometers as the default unit.]
>
> #### 3.  Add a bond between the alpha carbons
> Call the function
> ~~~​
> custom_bond.addBond(a1, a2)​
> ~~~
> {: .language-python}
> to add a bond between the two particles of interest that are indexed a1 and a2. Use the atom indices that you looked up earlier in this step.
>
> #### 4. Add the force to the system
Before the Simulation object is created, add the force to the system, using a line such as
> ~~~
> system.addForce(custom_bond)​
> ~~~
> {: .language-python}
> the same as the syntax that was used for adding the Monte Carlo Barostat.
>
> #### 5. Turn all the previous code into one function
> In your Python script, incorporate the above code into a function ​called AddSpringForce(a1, a2, k, r0)​ so that you may add a spring force between any two atoms with any provided spring constant k and equilibrium distance r0.  Write a docstring to document your function. Add a line to your script that calls the function.
> (Note​: As you become more experienced, you can always define the function first and then start add code to it, rather than writing the code first and then putting it into a function.)
>
> #### 6. Make sure your new script runs
> Once you have finished, run your script to see if your new feature works.  You might want to fist make a copy of your old trajectory so that you can compare the differences.
>
> #### 7. Run your code for several values of k
> Run the simulation using force constants of 10, 100, and 1000​. If your computer is fast enough, you can run for 100,000 or even more steps. If you are using the notebook to run simulations, make sure to set up a new system whenever you modify the spring constant. Otherwise you risk adding multiple custom forces to the same simulation, which will be confusing.
>
> #### 8. Visualize your results with MDTraj
> Use MDTraj to measure the distance vs. time between the pair of Cα atoms in your various trajectories​. How far apart do the atoms get in each simulation? You might also compare results from pulling on different atoms within the protein, such as alpha carbons on residue IDs 5 and 15.
>
>> ## Solution
>> Here is a completed code example for this exercise.
>> ~~~
>> from __future__ import print_function
>> from simtk.openmm import app
>> import simtk.openmm as mm
>> from simtk import unit
>> from sys import stdout
>>
>> # Load the PDB file into an object
>> pdb = app.PDBFile('trpcage.pdb')
>>
>> # Load the force field file into an object
>> forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')
>>
>> # Create system object using information in the force field:
>> # forcefield: contains parameters of interactions
>> # topology: lists of atoms, residues, and bonds
>> system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
>>     nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
>>     ewaldErrorTolerance=0.0005)
>>
>> # Create a Langevin integrator for temperature control
>> integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
>>     2.0*unit.femtoseconds)
>>
>> # Add a Monte Carlo barostat to the system for pressure control
>> barostat = mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin, 25)
>> system.addForce(barostat)
>>
>> # Use the CPU platform
>> platform = mm.Platform.getPlatformByName('CPU')
>>
>> ### If you want to add any forces to your System or modify any
>> ### of the existing forces, you should do it here - after the
>> ### System has been created, but before the Simulation is created.
>>
>> frc = mm.CustomBondForce("0.5*k*(r-r0)^2")
>> frc.addGlobalParameter("k", 1000)
>> frc.addGlobalParameter("r0", 2)
>>
>> frc.addBond(56, 200)
>> system.addForce(frc)
>>
>> # Create a Simulation object by putting together the objects above
>> simulation = app.Simulation(pdb.topology, system, integrator, platform)
>>
>> # Set positions in the Simulation object
>> simulation.context.setPositions(pdb.positions)
>>
>> # Minimize the energy of the system
>> print('Minimizing...')
>> print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
>> simulation.minimizeEnergy(maxIterations=20, tolerance=100)
>> print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
>>
>> # Initialize the random velocities of the system from a Maxwell-Boltzmann distribution
>> simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
>>
>> # Add reporters to the simulation object, which do things at regular intervals
>> # while the simulation is running.
>> # This reporter creates a DCD trajectory file
>> simulation.reporters.append(app.DCDReporter('trajectory-spr.dcd', 100))
>>
>> # This reporter prints information to the terminal
>> simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True,
>>     potentialEnergy=True, temperature=True, density=True, progress=True,
>>     remainingTime=True, speed=True, totalSteps=10000, separator='\t'))
>>
>> # Run the simulation itself
>> print('Running Production...')
>> simulation.step(10000)
>> print('Done!')
>> ~~~
>> {: .language-python}
> {: .solution}
{: .challenge}

> ## Exercise Option B - Zeroing out the charges (this is a more advanced exercise)
> This feature requires a bit more programming knowledge to implement compared to the Option A exercise.
>
> We might be interested in knowing the role that the partial charges play in keeping the protein stable. One might investigate this problem by zeroing out the charges on the protein. There are two ways to do this - you could edit the force field XML file and set all the charge parameters to zero, or you could set the charges to zero in the NonbondedForce (part of the System object). Here we will take the second route.
>
> #### 1. Create a new function
> Create a function called ZeroProteinCharges()​, so that you zero out the protein charges by calling this function. Make sure this function is defined after the system object has been created, but before the Simulation is created. Write the rest of the code into this function.
>
> #### 2. Find the NonbondedForce object and change things
> The System contains a number of Force objects of different kinds (HarmonicBondForcels, HarmonicAngleForce, NonbondedForce etc.) We want to change some parameters in the NonbondedForce. To do this, ​we will iterate over the list of Forces in the System​ and see which object is an instance of the NonbondedForce class.
>
> The method that iterates over the Forces is described in the API documentation. Start from the API documentation: ​http://docs.openmm.org/latest/api-python/index.html​, click on “Core Objects”, go to System, and look through the methods for the right one.
>
> After identifying the correct method, call it in your function​ to obtain the list of forces, and then write a Python loop to iterate over the items in this list. Inside the loop, check if the Force object is an instance of the NonbondedForce class using the syntax:
> ~~~
> if isinstance(force, mm.NonbondedForce):
> ~~~
> {: .language-python}
>
> #### 3. Set up your loop to count over the particles
Inside the if statement, you are now ready to iterate over the particles and set the charge parameters to zero. First determine the number of atoms in the protein (by reading the PDB file), and ​write a for loop from 0 (zero) up to the number of atoms in the protein​.
>
> #### 4. Use the documenation for the function to find the right parameters
Using the API documentation, look up the method in NonbondedForce for getting and setting the nonbonded parameters for particles.
>
> #### 5. Set the charges to zero
> e. We wish to set the charge parameters for the atoms in the protein to zero. Inside your for loop over protein atoms, first ​get the charge, sigma, and epsilon parameters for the atom (make sure to set three return values), and then set the parameters for that atom to 0 (zero), sigma, and epsilon respectively​.
>
> #### 6. Visualize your results
> Visualize the trajectory and analyze some interatomic distances to see if/how zeroing out the atomic charges affects the protein conformation. ​[One handy reaction coordinate for protein folding – at least in globular proteins – is radius of gyration. See if you can find the relevant MDTraj function for computing this quantity.]
>
>> ## Solution
>> Here is a completed code example for this exercise.
>> ~~~
>> from __future__ import print_function
>> from simtk.openmm import app
>> import simtk.openmm as mm
>> from simtk import unit
>> from sys import stdout
>>
>> # Load the PDB file into an object
>> pdb = app.PDBFile('trpcage.pdb')
>>
>> # Load the force field file into an object
>> forcefield = app.ForceField('amber99sbildn.xml', 'tip3p.xml')
>>
>> # Create system object using information in the force field:
>> # forcefield: contains parameters of interactions
>> # topology: lists of atoms, residues, and bonds
>> system = forcefield.createSystem(pdb.topology, nonbondedMethod=app.PME,
>>     nonbondedCutoff=1.0*unit.nanometers, constraints=app.HBonds, rigidWater=True,
>>     ewaldErrorTolerance=0.0005)
>>
>> # Create a Langevin integrator for temperature control
>> integrator = mm.LangevinIntegrator(300*unit.kelvin, 1.0/unit.picoseconds,
>>     2.0*unit.femtoseconds)
>>
>> # Add a Monte Carlo barostat to the system for pressure control
>> system.addForce(mm.MonteCarloBarostat(1*unit.atmospheres, 300*unit.kelvin, 25))
>>
>> # Use the CPU platform
>> platform = mm.Platform.getPlatformByName('CPU')
>>
>> ### If you want to add any forces to your System or modify any
>> ### of the existing forces, you should do it here - after the
>> ### System has been created, but before the Simulation is created.
>>
>> atoms = list(pdb.topology.atoms())
>>
>> for f in system.getForces():
>>     if isinstance(f, mm.NonbondedForce):
>>         print("Found nonbonded force")
>>         for i in range(f.getNumParticles()):
>>             if atoms[i].residue.name != 'HOH':
>>                 chg, sig, eps = f.getParticleParameters(i)
>>                 f.setParticleParameters(i, 0.0, sig, eps)
>>
>> # Create a Simulation object by putting together the objects above
>> simulation = app.Simulation(pdb.topology, system, integrator, platform)
>>
>> # Set positions in the Simulation object
>> simulation.context.setPositions(pdb.positions)
>>
>> # Minimize the energy of the system
>> print('Minimizing...')
>> simulation.minimizeEnergy()
>>
>> # Initialize the random velocities of the system from a Maxwell-Boltzmann distribution
>> simulation.context.setVelocitiesToTemperature(300*unit.kelvin)
>>
>> # Add reporters to the simulation object, which do things at regular intervals
>> # while the simulation is running.
>> # This reporter creates a DCD trajectory file
>> simulation.reporters.append(app.DCDReporter('trajectory-noChg.dcd', 100))
>>
>> # This reporter prints information to the terminal
>> simulation.reporters.append(app.StateDataReporter(stdout, 100, step=True,
>>     potentialEnergy=True, temperature=True, density=True, progress=True,
>>     remainingTime=True, speed=True, totalSteps=10000, separator='\t'))
>>
>> # Run the simulation itself
>> print('Running Production...')
>> simulation.step(10000)
>> print('Done!')
>> ~~~
>> {: .language-python}
> {: .solution}
{: .challenge}

## Commit your changes on the local machine, and push them to your repository fork on GitHub.
If you have been working in the jupyter notebook to write your code modifications, copy and paste you completed code back into a text file.  Save your code.  Make sure the file extension is .py.

Now you want to share your modified code by pushing it to your fork on GitHub.  

The commands are:
~~~
$ git commit -m “Type any commit message you want”
$ git push origin master
~~~
{: .language-bash}

Here, ​origin ​stands for the remote repository, i.e. your repository fork on GitHub. You can request to display the remote repository names and their URLs using the command:
~~~
$ git remote -v
~~~
{: .language-bash}

Your changes are now uploaded to your repository fork on GitHub, but they are not incorporated into the main code on the MolSSI GitHub.

## Create a pull request on GitHub for your feature to be added to the original repository.
At this point, your code has a new feature that isn’t part of the original repository; you may request the owner to incorporate your code changes using the pull request button. Click the pull request button to initiate a conversation thread; the owner may accept the pull request right away (updating the original repository with your codes), or ask you to make more changes.

{% include links.md %}
