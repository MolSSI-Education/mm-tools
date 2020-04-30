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

Now that we have our force field defined, we are ready to do a simulation! But, what does that mean and what kind of simulation will we be doing? You may remember from the previous unit that molecular mechanics simulations can be used to predict the time dependent properties of a system. We will use a simulation method called **molecular dynamics** (MD). In molecular dynamics simulations, we simulate molecules in time by calculating the forces on a molecule and updating atom positions based on those forces. The output of a molecular dynamics simulation is a trajectory. The trajectory can be visualized to show the movement of a system through time (like watching a movie).

Force is equal to the negative gradient of potential energy. MD simulations use the potential energy function we discussed previously, and Newton's second law to calculate positions of atoms. 

$$ \vec{F} = - \nabla U $$

$$ \vec{F} = m \cdot \vec{a} = m \cdot \frac{d\vec{v}}{dt} = m \cdot \frac{d^{2}\vec{r}}{dt^{2}}$$

Looking at these two equations, we already know the mass for each atom, and with our force field defined, we can calculate the force. This means that we can calculate the acceleration or velocity (and thus a new position) for each atom based on the force field function and mass of the molecule.

$$ \vec{a}_{i}(t_{0}) = \frac{1}{m}\vec{F}_{i}(t_{0}) = -\frac{1}{m}_{i}(\nabla U)_{i} $$

We know that we can calculate positions ($$\vec{x}$$) based on previous positions and acceleration and an amount of time ($$\Delta t$$)

$$\vec{x}_i (t_{0} + \Delta t) = \vec{x}_i (t_{0}) + \vec{v}_{i} (t_{0}) \Delta t + \frac{1}{2} \vec{a}_{i}(t_{0}) \Delta t^{2} $$

This is the general approach to molecular dynamics. You have some atoms which have initial positions **and** velocities when you start your simulation. You can use these positions to calculate the potential energy and the forces on each atom (because force is the negative gradient of potential energy). With the velocities of the atoms, you can calculate new positions. This process is repeated over and over again, and each iteration is typically called a **timestep**. 

The typical time step size ($$\Delta t) represents 1 femtosecond of time and is based on X-H bond vibration frequencies. A timestep that is too small will lead to an inefficient simulation, while a time step that is too large will be inaccurate or crash your software. 

MD simulations assumes that the [Erogodic hypothesis](https://en.wikipedia.org/wiki/Ergodic_hypothesis) is true, meaning that the time average calculated from the simulation is equal to the ensemble average. For this to be true, the simulation has to sample many conformations. Most simulations have millions to trillions of time steps. 

### Challenges and limitations

### Running your first simulations

{% include links.md %}

