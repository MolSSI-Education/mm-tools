---
title: "Introduction"
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

## What does "molecular mechanics" mean?

Computational chemistry is a field which uses calculations performed on computers to predict real world properties of molecules. Because we can't have "real" molecules in a computer, we have to approximate their behavior using mathematical models in simulations. Within the field of computational chemistry, the types of simulation are grouped into two broad categories based on the type of physics the models are based on. Quantum chemistry simulations are based on quantum mechanics (Schrodinger's equation if you have taken quantum chemistry). The other broad category of computational chemistry is based on **molecular mechanics**, or classical(Newtonian) physics. This is essentially the same kind of physics you learn in introductory physics course where you might predict the path of a ball thrown through the air. 

Quantum mechanics calculations are more accurate and detailed than molecular mechanics calculations. However, quantum calculations take longer (are more "expensive") and are limited in what can be studied. Molecular mechanics methods are most often used to study the time dependent properties of molecules. The figure below illustrates this trade off between accuracy and cost. The lower left corner shows simulations which use classical physics. They are less expensive, but also less detailed. The upper right represents quantum mechanics methods - they are very accurate, but also very expensive. This workshop focuses on molecular mechanics simulations.

<img src="/fig/simulation-scale.png" height=500>

## The "force field"

As stated above, all computational simulations depend on a mathematical model to describe the molecule. The full expression of this model is often called a **force field**  or potential energy function in molecular mechanics simulations. This force field describes the energy associated with molecular movements such as bond stretching, angle bending, or dihedral angle rotation. If you are a chemistry student, you have likely discussed many of these molecular motions in your classes (who can forget talking about the cis-trans isomerism of butane in organic chemistry class?). The force field describes the energies associated with these movements mathematically. The potential energy is commonly represented by the letter $$U$$, and the total potential energy or "force field" is the sum of terms related to bond stretching, angle bending, torsional rotation (that cis-trans isomerism), electrostatic interaction, and nonbonded Van der Waals interactions.

$$ U = U_{bond} + U_{angle} + U_{torsion} + U_{elec} + U_{vdw} $$

Below we show commonly used forms from these terms. You will recognize many of these from other classes you have taken:

$$ U_{bond} = \frac{1}{2}k_{l}(l-l_{eq})^2 $$

$$ U_{angle} = \frac{1}{2}k_{\theta}(\theta - \theta_{eq})^2 $$

Note that both bonds and angles are commonly described using a harmonic potential (also used to describe spring-mass systems.) The variable $$l$$ represents the distance between two bonded atoms. The parameters $$k_{l}$$ and $$l_{eq}$$ represent the bond stiffness and equilibrium bond length respectively. For example, two double bonded carbons will have a shorter bond length ($$l_{eq}$$) and stiffer bond (stiffer spring, or higher $$k_{l}$$) than single bonded carbons.

$$ U_{torsion} = \sum_{n=1}^{n_{max}}{U_{n}[1+cos(n\phi-\gamma_{n})]} $$

$$ U_{elec} = k_{e}\frac{q_{i}q_{j}}{r_{ij}} $$

$$ U_{vdw} = 4\epsilon_{ij}[(\frac{\alpha_{ij}}{r_{ij}})^{12} - (\frac{\alpha_{ij}}{r_{ij}})^{6}] $$

## Getting ready to simulate molecules

We are going to perform simulations of molecules using the software [OpenMM](http://openmm.org/). As discussed, we must first set our model, or force field. Consider the following parameters which we can use in the equations above to describe ethane.

| Interaction       | Parameters                                                                            |
|-------------------|---------------------------------------------------------------------------------------|
|C-C bond           | $$l_{eq} = 0.1538 nm$$, $$k_{l} = 1.946 x 10^5 kJ/mol \cdot nm^2$$                    |
|C-H bond           | $$l_{eq} = 0.1097 nm$$, $$k_{l} = 3.146 x 10^6 kJ/mol \cdot nm^2$$                    |
|H-C-H angle        | $$\theta_{eq} = 1.878$$ $$rad$$, $$k_{\theta} = 326.0 kJ/mol \cdot rad^2$$            |
|H-C-C angle        | $$\theta_{eq} = 1.916$$ $$rad$$, $$k_{\theta} = 391.8 kJ/mol \cdot rad^2$$            |
|H-C-C-H torsion    | $$\gamma_{3} = 0$$ $$rad$$, $$U_{3} = 0.5021 kJ/mol$$                                 |
|C atom (nonbonded) | $$q_{C} = -0.0951e$$, $$\alpha_{C} = 0.3398$$ $$nm$$, $$\epsilon_{C} = 0.4510 kJ/mol$$|
|H atom (nonbonded) | $$q_{H} = 0.0317e$$, $$\alpha_{H} = 0.2600$$ $$nm$$, $$\epsilon_{H} = 0.0870 kJ/mol$$ |


You can imagine how we could fill in these parameters to the equations given in the section above, but we also have to get them into a format that our force field understands. OpenMM understands forcefields in a file format called XML. The XML file gives the force field parameters, and also tells the program the connectivity (or **topology**) of the atoms. Before we can simulate a molecule with OpenMM, we must first make an XML file to tell the software about the molecule we want to simulate. We will build the XML files by hand for this exercise, however, you will usually not have to do this by hand (so don't let this part scare you away).

Consider the XML file for ethane:

~~~
<ForceField>
 <AtomTypes>
  <Type name="0" class="c3" element="C" mass="12.01078"/>
  <Type name="1" class="hc" element="H" mass="1.007947"/>
 </AtomTypes>
 <Residues>
  <Residue name="ETH">
   <Atom name="C1" type="0"/>
   <Atom name="H11" type="1"/>
   <Atom name="H12" type="1"/>
   <Atom name="H13" type="1"/>
   <Atom name="C2" type="0"/>
   <Atom name="H21" type="1"/>
   <Atom name="H22" type="1"/>
   <Atom name="H23" type="1"/>
   <Bond atomName1="C1" atomName2="H11"/>
   <Bond atomName1="C1" atomName2="H12"/>
   <Bond atomName1="C1" atomName2="H13"/>
   <Bond atomName1="C1" atomName2="C2"/>
   <Bond atomName1="C2" atomName2="H21"/>
   <Bond atomName1="C2" atomName2="H22"/>
   <Bond atomName1="C2" atomName2="H23"/>
  </Residue>
 </Residues>
 <HarmonicBondForce>
  <Bond class1="c3" class2="c3" length="0.15380" k="1945727.36"/>
  <Bond class1="c3" class2="hc" length="0.10970" k="3145687.56"/>
 </HarmonicBondForce>
 <HarmonicAngleForce>
  <Angle class1="c3" class2="c3" class3="hc" angle="1.91637152" k="391.756288"/>
  <Angle class1="hc" class2="c3" class3="hc" angle="1.87762521" k="326.01728"/>
 </HarmonicAngleForce>
 <PeriodicTorsionForce>
  <Proper class1="hc" class2="c3" class3="c3" class4="hc" periodicity1="3" phase1="0.0" k1="0.50208"/>
 </PeriodicTorsionForce>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="0" charge="-0.094100" sigma="0.3397710" epsilon="0.4510352"/>
  <Atom type="1" charge="0.031700" sigma="0.2600177" epsilon="0.0870272"/>
 </NonbondedForce>
</ForceField>
~~~
{: .language-xml}

Let's look at each section of the file.

First, we tell OpenMM that this is a force field file by using this syntax:

~~~
<ForceField>
~~~
{: .language-xml}

Next, we identify atoms in our molecules by setting atom types. In our ethane molecule, we have two carbon atoms and six hydrogen atoms. The two carbon atoms are equivalent, and the six hydrogen atoms are all also the same. Therefore, we have two types of atoms. We tell OpenMM this by using starting an AtomTypes section `<AtomTypes>`, and listing our atom types.

~~~
 <AtomTypes>
  <Type name="0" class="c3" element="C" mass="12.01078"/>
  <Type name="1" class="hc" element="H" mass="1.007947"/>
 </AtomTypes>
~~~
{: .language-xml}

Next, we define the connectivity or topology of our system. We might simulate many molecules in a simulation, and each molecule will be called a 'residue'. We will tell OpenMM what atoms belong to each molecule, and which other atoms each atom is bonded to. Because we are only describing one molecule, we have only one residue.

~~~
<Residues>
  <Residue name="ETH">
   <Atom name="C1" type="0"/>
   <Atom name="H11" type="1"/>
   <Atom name="H12" type="1"/>
   <Atom name="H13" type="1"/>
   <Atom name="C2" type="0"/>
   <Atom name="H21" type="1"/>
   <Atom name="H22" type="1"/>
   <Atom name="H23" type="1"/>
   <Bond atomName1="C1" atomName2="H11"/>
   <Bond atomName1="C1" atomName2="H12"/>
   <Bond atomName1="C1" atomName2="H13"/>
   <Bond atomName1="C1" atomName2="C2"/>
   <Bond atomName1="C2" atomName2="H21"/>
   <Bond atomName1="C2" atomName2="H22"/>
   <Bond atomName1="C2" atomName2="H23"/>
  </Residue>
 </Residues>
~~~
{: .language-xml}

Finally, we tell OpenMM the parameters for our force field equation.

~~~
<HarmonicBondForce>
  <Bond class1="c3" class2="c3" length="0.15380" k="1945727.36"/>
  <Bond class1="c3" class2="hc" length="0.10970" k="3145687.56"/>
 </HarmonicBondForce>
 <HarmonicAngleForce>
  <Angle class1="c3" class2="c3" class3="hc" angle="1.91637152" k="391.756288"/>
  <Angle class1="hc" class2="c3" class3="hc" angle="1.87762521" k="326.01728"/>
 </HarmonicAngleForce>
 <PeriodicTorsionForce>
  <Proper class1="hc" class2="c3" class3="c3" class4="hc" periodicity1="3" phase1="0.0" k1="0.50208"/>
 </PeriodicTorsionForce>
 <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
  <Atom type="0" charge="-0.094100" sigma="0.3397710" epsilon="0.4510352"/>
  <Atom type="1" charge="0.031700" sigma="0.2600177" epsilon="0.0870272"/>
 </NonbondedForce>
</ForceField>
~~~
{: .language-xml}

The last line tells OpenMM we are done telling it about the ForceField. We have prepared our file for OpenMM, and are now ready to simulate our molecule.

> ## Exercise
> Prepare an XML file for the butane molecule. The butane molecule is a little more complicated. Consider the following:
> 
> How many atom types are there in butane? Are there more than ethane? (think about whether carbons and hydrogens are equivalent)
> Are there more bond types? Ethane had C-H bonds and C-C bonds. Are the bonds in butane different?
>
> Once you have thought about these questions, consider the parameters below. Use these parameters and the ones for ethane to build an XML file for butane.
>
> | Interaction       | Parameters                                                                            |
> |-------------------|---------------------------------------------------------------------------------------|
> |C-C-C angle        | $$\theta_{eq} = 1.946$$ $$rad$$, $$k_{\theta} = 543.0 kJ/mol \cdot rad^2$$            |
> |H-C-C-C torsion    | $$\gamma_{3} = 0$$ $$rad$$, $$U_{3} = 0.3347 kJ/mol$$                                 |
> |C-C-C-C torsion    | $$\gamma_{1} = 0$$ $$rad$$, $$U_{1} = 0.4602 kJ/mol$$, $$\gamma_{2} = 3.146$$ $$rad$$, $$U_{2} = 1.2134 kJ/mol$$, $$\gamma_{3} = 0$$ $$rad$$, $$U_{3} = 0.5439 kJ/mol    |
> |"outer" C atom (nonbonded) | $$q_{C} = -0.0932e$$, $$\alpha_{C}$$ and $$\epsilon_{C}$$ same as ethane|
> |"inner" C atom (nonbonded) | $$q_{C} = -0.0814e$$, $$\alpha_{C}$$ and $$\epsilon_{C}$$ same as ethane|
> |"outer" H atom (nonbonded) | $$q_{H} = 0.0324e$$, $$\alpha_{H}$$ and $$\epsilon_{H}$$ same as ethane|
> |"inner" H atom (nonbonded) | $$q_{H} = 0.0387e$$, $$\alpha_{H}$$ and $$\epsilon_{H}$$ same as ethane|
>
> Copy this XML File to get started
> ~~~
> <ForceField>
>  <AtomTypes>
>   <Type name="0" class="c3" element="C" mass="12.01078"/>
>   <Type name="1" class="c3" element="C" mass="12.01078"/>
>   <Type name="2" class="hc" element="H" mass="1.007947"/>
>   <Type name="3" class="hc" element="H" mass="1.007947"/>
>  </AtomTypes>
>  <Residues>
>   <Residue name="NBU">
>    <Atom name="C1" type="0"/>
>    <Atom name="H11" type="2"/>
>    <Atom name="H12" type="2"/>
>    <Atom name="H13" type="2"/>
>    <Atom name="C2" type="1"/>
>    <Atom name="H21" type="3"/>
>    <Atom name="H22" type="3"/>
>    <Atom name="C3" type="1"/>
>    <Atom name="H31" type="3"/>
>    <Atom name="H32" type="3"/>
>    <Atom name="C4" type="0"/>
>    <Atom name="H41" type="2"/>
>    <Atom name="H42" type="2"/>
>    <Atom name="H43" type="2"/>
>    ***FILL IN THE TOPOLOGY***
>   </Residue>
>  </Residues>
>  <HarmonicBondForce>
>   ***FILL IN THE BOND PARAMETERS***
>  </HarmonicBondForce>
>  <HarmonicAngleForce>
>   ***FILL IN THE ANGLE PARAMETERS***
>  </HarmonicAngleForce>
>  <PeriodicTorsionForce>
>   ***FILL IN THE TORSION PARAMETERS***
>  </PeriodicTorsionForce>
>  <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
>   ***FILL IN THE NONBONDED PARAMETERS***
>  </NonbondedForce>
> </ForceField>
> ~~~ 
> {: .language-xml}
>
>> ## Solution
>>
>> ~~~
>> <ForceField>
>>  <AtomTypes>
>>   <Type name="0" class="c3" element="C" mass="12.01078"/>
>>   <Type name="1" class="c3" element="C" mass="12.01078"/>
>>   <Type name="2" class="hc" element="H" mass="1.007947"/>
>>   <Type name="3" class="hc" element="H" mass="1.007947"/>
>>  </AtomTypes>
>>  <Residues>
>>   <Residue name="NBU">
>>    <Atom name="C1" type="0"/>
>>    <Atom name="H11" type="2"/>
>>    <Atom name="H12" type="2"/>
>>    <Atom name="H13" type="2"/>
>>    <Atom name="C2" type="1"/>
>>    <Atom name="H21" type="3"/>
>>    <Atom name="H22" type="3"/>
>>    <Atom name="C3" type="1"/>
>>    <Atom name="H31" type="3"/>
>>    <Atom name="H32" type="3"/>
>>    <Atom name="C4" type="0"/>
>>    <Atom name="H41" type="2"/>
>>    <Atom name="H42" type="2"/>
>>    <Atom name="H43" type="2"/>
>>    <Bond atomName1="C1" atomName2="H11"/>
>>    <Bond atomName1="C1" atomName2="H12"/>
>>    <Bond atomName1="C1" atomName2="H13"/>
>>    <Bond atomName1="C1" atomName2="C2"/>
>>    <Bond atomName1="C2" atomName2="H21"/>
>>    <Bond atomName1="C2" atomName2="H22"/>
>>    <Bond atomName1="C2" atomName2="C3"/>
>>    <Bond atomName1="C3" atomName2="H31"/>
>>    <Bond atomName1="C3" atomName2="H32"/>
>>    <Bond atomName1="C3" atomName2="C4"/>
>>    <Bond atomName1="C4" atomName2="H41"/>
>>    <Bond atomName1="C4" atomName2="H42"/>
>>    <Bond atomName1="C4" atomName2="H43"/>
>>   </Residue>
>>  </Residues>
>>  <HarmonicBondForce>
>>   <Bond class1="c3" class2="c3" length="0.15380" k="1945727.36"/>
>>   <Bond class1="c3" class2="hc" length="0.10970" k="3145687.56"/>
>>  </HarmonicBondForce>
>>  <HarmonicAngleForce>
>>   <Angle class1="c3" class2="c3" class3="c3" angle="1.94621665" k="542.982784"/>
>>   <Angle class1="c3" class2="c3" class3="hc" angle="1.91637152" k="391.756288"/>
>>   <Angle class1="hc" class2="c3" class3="hc" angle="1.87762521" k="326.01728"/>
>>  </HarmonicAngleForce>
>>  <PeriodicTorsionForce>
>>   <Proper class1="c3" class2="c3" class3="c3" class4="c3" periodicity1="3" phase1="0.0" k1="0.5439" 
>>           periodicity2="2" phase2="3.1416" k2="1.2134" periodicity3="1" phase3="0.0" k3="0.4602"/>
>>   <Proper class1="c3" class2="c3" class3="c3" class4="hc" periodicity1="3" phase1="0.0" k1="0.3347"/>
>>   <Proper class1="hc" class2="c3" class3="c3" class4="hc" periodicity1="3" phase1="0.0" k1="0.50208"/>
>>  </PeriodicTorsionForce>
>>  <NonbondedForce coulomb14scale="0.833333" lj14scale="0.5">
>>   <Atom type="0" charge="-0.0932" sigma="0.3397710" epsilon="0.4510352"/>
>>   <Atom type="1" charge="-0.0814" sigma="0.3397710" epsilon="0.4510352"/>
>>   <Atom type="2" charge="0.0324" sigma="0.2600177" epsilon="0.0870272"/>
>>   <Atom type="3" charge="0.0387" sigma="0.2600177" epsilon="0.0870272"/>
>>  </NonbondedForce>
>> </ForceField>
>> ~~~
>> {: .language-xml}
>>
> {: .solution}
{: .challenge}




{% include links.md %}

