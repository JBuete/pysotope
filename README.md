# Hello, and welcome to Pysotope!

Pysotope is a basic implementation of a Python class for the known isotopes.

## Introduction

There are two main objects that are worth knowing if you are going to use this 
module for anything nuclear physics related. Your primary interface to the information
contained in the module will be through the Isotope class.

### Isotopes

Isotopes can be created in two main ways, by their typical written definition,

``` python
>>> pb = Isotope("208pb")
>>> pb = Isotope("pb208") 
```

or by supplying the mass and charge numbers,

```python
​```python
>>> pb = Isotope(82, 208)
>>> pb = Isotope(208, 82)
​```
```

Once assigned to a variable you can begin to perform basic linear operations with the isotopes,

```python
>>> alpha = Isotope("4he")
>>> pb + alpha
Isotope: 212Po
>>> pb - alpha
Isotope: 204Hg
```

Each isotope contains:

* its charge and mass number,
* its mass in micro-amu,
* its binding energy,
* its mass excess,

and all associated uncertainties. This information is retrieved from the 2016 Atomic
Mass Evaluation table, the NUBASE2016 evaluation, and calculations from FRLDM 2012,  when the module is first imported and accessed whenever an isotope is created. You can access this information yourself through the following functions

```python
>>> pb.get_mass()  # returns the mass in micro-amu
207976651.918  
```

```python
>>> pb.get_mass_uncertainty()  # returns the uncertainty in the mass in micro-amu
1.231
```

```python
>>> pb.get_binding_energy()  # returns the binding energy in keV
7867.453
```

```python
>>> pb.get_be_uncertainty()  # returns the uncertainty in the binding energy in keV
0.006
```

```python
>>> pb.get_mass_excess()  # returns the mass excess in keV
-21748.598
```

```python
>>> pb.get_me_uncertainty()  # returns the uncertainty in the mass excess in keV
1.148
```

```python
>>> pb.get_z()  # returns the charge number
82
```

```python
>>> pb.get_a()  # returns the mass number
208
```

You can also access a simply print out of the information by using the provided get_info()
function,

```python
>>> pb.get_info()  # prints some basic information about the isotope
Isotope: 208Pb
Mass = 207.9766519180 +/- 0.0000012310 u
M_ex = -21.748598 +/- 0.001148 MeV
B/A  = 7.867453 +/- 0.000006 MeV
S_p  = 8.005489 MeV
S_2p = 15.385295 MeV
S_n  = 7.370052 MeV
S_2n = 14.109828 MeV
```

If you want to print out the isotope's name or representation you can do that in two different
ways. The first is to convert it to a string implicitly using the standard str() function

```python
>>> print(str(pb))
208Pb
```

You can also use python's format operator

```python
>>> print("{}".format(pb))
208Pb
```

This may not seem like the most efficient way to print out the isotope's name, but being able to
use the format operator allows very clean interactions with the next important object in
Pysotope, Reactions.

### Reactions

Reactions allow you to create a reaction and get key information, such as the Q-value, without
going through the hassle of actually doing the maths. You can initialise a reaction in the
following way

```python
>>> my_cool_reaction = Reaction("6li + 207pb -> 208pb + p + a")
```

The reaction string only requires the names of the reactants and products, merely
separated by an arrow "->". Notice that we can also use standard shortcuts like a for an alpha
particle (4He) or p for a proton (1H). You can also use d (2H), t (3H), and n (1n). (This will
also work when defining isotopes) 

You can also specify the beam energy in MeV as a second (optional) argument to the initialisation

```python
>>> my_cool_reaction = Reaction("6li + 207pb -> 208pb + p + a", 29)
```

or as an argument to any function that returns a reaction value that depends on the beam energy,
we will discuss those later. 

In the case where you are checking multiple reactions automatically or want to write a
generalised piece of code that checks reaction values I recommend using the format operator as this enables easy generation of reactions, 

```python
>>> target = Isotope("207pb")
>>> projectile = Isotope("6li")
>>> n = Isotope("n")
>>> Reaction("{} + {} -> {} + {}".format(projectile, target, target + n, projectile - n))
Reaction: 6li + 207pb -> 208pb + 5li
```

Reactions also allow you to access the $q_{value}$, and other information about the reaction via the
following functions

```python
>>> my_cool_reaction.get_qvalue()
3.670631556440236
```

```python
>>> my_cool_reaction.get_q_opt()
0.0 
```

```python
>>> my_cool_reaction.get_cm_energy()
28.181005076906924
```

Notice that both $q_{opt}$ and the centre-of-mass energy depend on the energy of the beam, which we
supplied when we first initialised the reaction, but if we wanted to check the centre of mass
energy for a different beam energy we can simply add it as an argument

```python
>>> my_cool_reaction.get_cm_energy(31)
30.124522668417747
```

but this doesn't change the internal ebeam.

```python
>>> my_cool_reaction.get_ebeam()
29
```

If instead we wanted to change the beam energy for the reaction (or we forgot to do so in the
first place we can simply use set_ebeam

```python
>>> my_cool_reaction.set_ebeam(31)
>>> my_cool_reaction.get_ebeam()
31
```



## Sources

The data used in Pysotope comes from the following sources:

* Wang, M., Audi,  G., Kondev, F. G., Huang, W. J., Naimi, S., & Xu, X. (2017). The  AME2016 atomic mass evaluation (II). Tables, graphs and references. *Chinese Physics C*, *41*(3), 030003.
* Audi, G., Kondev, F. G., Wang, M., Huang, W. J., & Naimi, S. (2017). The NUBASE2016 evaluation of nuclear properties. *Chinese Physics C*, *41*(3), 030001.
* Möller, P., Sierk, A. J., Ichikawa, T., & Sagawa, H. (2016). Nuclear ground-state masses and deformations: FRDM (2012). *Atomic Data and Nuclear Data Tables*, *109*, 1-204.