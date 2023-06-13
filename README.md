
# Integration of Tinker and LFMM method

 This documentation is specific for the modified version of <a href="http://dasher.wustl.edu/tinker/">Tinker v6.3.3</a> that can run calculations according to the Ligand Field Molecular Mechanics (LFMM) method as developed by Prof. Rob J. Deeth – University of Warwick, Inorganic Computational Chemistry Group (ICCG).

## Integration Strategy
TinkerLFMM is the result of integrating the capabilities of Dommino (the original software developed by Prof. Rob J. Deeth for Ligand field molecular mechanics - LFMM) into Tinker and, therefore, it  consists of two main building blocks:

* *Tinker v6.3.3*: takes care of all regular molecular mechanics energy term, and behaves as sole interface with the user.
* *ICCGLFSE*, ex-*Dommino*: calculates Ligand Field Stabilization Energy.

Presently (June 2023), this repository contains only the modifier Tinker code.

According to the LFMM method, the potential energy is calculated from an ensemble of regular MM-style energy terms and additional contributions describing each metal center. These extra terms consist of metal-ligands bond stretches, ligand-ligand interaction, electron pairing energy and ligand field stabilization energy. All such contributions to the total potential energy of the system, as well as  their contributions to the gradient and Hessian, are calculated by TinkerLFMM whenever the modelling task requires so, calling the `iccglfse` routine (see #Compilation).

## Limitations
The use of LFMM potential is compatible only with some of the methods that are available in the Tinker package. The following Tinker routines that have been tested with LFMM calculations:
* `analyze`
* `minimize`
* `optimize`
* `newton`
* `monte`
* `dynamic`
* `testgrad`
* `testhess`

Since the definition of the LFMM potential requires that for an LFMM center the metal, the ligand atoms, and, depending on the ligand type, also subsidiary ligand atoms, are treated as a single unit, the use of `GROUP`, `ACTIVE`, and `INACTIVE` keywords may lead to unexpected results if such options are used to alter the consistency of the LFMM center.

Potential energy terms for LFMM are not compatible with the `SMOOTHING` option.

Anything that is not explicitly reported here should be considered as not compatible.

## Compilation
Both Tinker and ICCGLFSE are written in Fortran (mixing Fortran 77 and 95) and are compiled independently. The present version of Tinker contains a placeholder routine `iccglfse` that is overwritten upon integration. The integrated TinkerLFMM is obtained by linking the compiled Tinker object files (*.o) with the corresponding object files of ICCGLFSE, which contain the actual `iccglfse` routine, and which should have been previously collected into a library during the compilation step.

Makefiles are provided to allow a make-based building and installation, but you may have to edit the `make/Makefile_TinkerLFMM` and `make/Makefile` to use the compile or your choice:

```
make -f make/Makefile_TinkerLFMM
make -f make/Makefile_TinkerLFMM install
```

## Test and Examples
Test that compare the results of TinkerLFMM against the reference results from [DommiMOE](https://doi.org/10.1002/jcc.20137) are available in [a dedicated repository](https://github.com/marco-foscato/TinkerLFMM_tests).

## Usage
The following lists contain all the documentation needed to use the LFMM method in Tinker. For worked out examples, please see [the repository with tests](https://github.com/marco-foscato/TinkerLFMM_tests).

### Keywords for LFMM calculations

* `LFMMTERM [ONLY/NONE]` This keyword controls use of all energy terms related to transition metal ligand field potentials developed by Rob Deeth's group. In the absence of a modifying option, this keyword makes Tinker use all the potentials associated with this approach. The NONE option turns OFF use of all these potential energy terms. The ONLY option turns OFF all potential energy terms except for these ones. This keyword corresponds to the combination of HARM-LFMMTERM, MORSE-LFMMTERM, LLR-LFMMTERM, VWLL-LFMMTERM, PAIR-LFMMTERM, and LFSE-LFMMTERM.

* `HARM-LFMMTERM [ONLY/NONE]` This keyword controls use of the harmonic metal-ligand bond stretching term of LFMM method. In the absence of a modifying option, this keyword turns ON use of the potential energy term. The NONE option turns OFF use of this potential energy term. The ONLY option turns OFF all potential energy terms except for this one.

* `MORSE-LFMMTERM [ONLY/NONE]` This keyword controls use of the Morse metal-ligand bond stretching of term LFMM method. In the absence of a modifying option, this keyword turns ON use of the potential energy term. The NONE option turns OFF use of this potential energy term. The ONLY option turns OFF all potential energy terms except for this one.

* `LLR-LFMMTERM [ONLY/NONE]` This keyword controls use of the Ligand-Ligand pure repulsion term of LFMM method. In the absence of a modifying option, this keyword turns ON use of the potential energy term. The NONE option turns OFF use of this potential energy term. The ONLY option turns OFF all potential energy terms except for this one.

* `VWLL-LFMMTERM [ONLY/NONE]` This keyword controls use of the Ligand-Ligand Lennard-Jones-style interaction term of LFMM method. In the absence of a modifying option, this keyword turns ON use of the potential energy term. The NONE option turns OFF use of this potential energy term. The ONLY option turns OFF all potential energy terms except for this one.

* `PAIR-LFMMTERM [ONLY/NONE]` This keyword controls use of the electron pairing energy term of LFMM method. In the absence of a modifying option, this keyword turns ON use of the potential energy term. The NONE option turns OFF use of this potential energy term. The ONLY option turns OFF all potential energy terms except for this one.

* `LFSE-LFMMTERM [ONLY/NONE]` This keyword controls use of the ligand field stabilization energy term of LFMM method. In the absence of a modifying option, this keyword turns ON use of the potential energy term. The NONE option turns OFF use of this potential energy term. The ONLY option turns OFF all potential energy terms except for this one.

* `LFMMCENTER [3 integer]` This keyword allows to define a single LFMM center which is identified by a single metal atom. The three integers specify i) the atom number of the metal in the atom list, ii) the number of d electrons, and iii) the spin state (0=low; 1=high; 2=intermediate). Coordinating atoms are identified by Tinker on the basis of the connectivity defined in the input XYZ file.

* `BARY-SJ-LFSE [integer]` Sets the barycenter of Shaffer and Jorgensen d/s-mixing. Default value is 2.

* `NO-SJ-LFSE` Disable Shaffer and Jorgensen d/s-mixing. By default d/s-mixing is active.

* `DERIV-FD-LFSE` Causes LFSE routine to calculate the first derivatives of the ligand field stabilization energy numerically by finite differences. By default analytical derivatives are used instead.

* `DISPL-FD-LFSE [real]` Sets the number of displacement used by LFSE routine for numerical differentiation. Default value is 2.

* `STEP-FD-LFSE [real]` Sets the step size used for finite differences method by LFSE routine. Default value is 0.02d0.

* `TNK-FDSD-LFSE` Causes Tinker to calculate second derivatives numerically by finite differences approximation of the LFSE gradient. This approach is meant for testing and debugging.

* `STB-ENG-FORM-LFSE` Causes LFSE routine to use Burdett's formulation of Molecular Orbital Stabilization Energy (MOSE), see Coord. Chem. Rev. 212 (2001) 11-34. By default such condition is not applied.

* `1210VDWTERM` Causes Tinker to use 12-10 van der Waals functional form for all the pairs of atom types for which such parameters are available. By default 12-10 van der Waals interaction are not calculated.

* `DIST-DEP-DIELECTRIC` Causes Tinker to use distance dependent dielectric for charge-charge interactions. When this keyword is used, charge-charge energies are calculated as k*(q_1*q_2)/(r_12+buffer)**n where n = 2. By default n = 1.

* `CHG-TAPER-TYPE [integer]` Selection of energy switch function for charge-charge interactions if near the cutoff distance. Specify '1' to require standard tapering type involving both multiplicative and additive switch functions, or '2' to require only multiplicative switch function. Default tapering type is 1.

* `LOW-ML-RESTRAINT [real]` Sets the reduction to the ideal M-L bond length of the Morse-type, metal-ligand bond stretching term of the LFMM method. This value controls the interatomic distance at which the restraint function is applied to destabilize short M-L bond length. Default value is 0.3d0.

* `UP-ML-RESTRAINT [real]` Sets the increment to the ideal M-L bond length of the Morse-type, metal-ligand bond stretching term of the LFMM method. This value controls the interatomic distance at which the restraint function is applied to destabilize long M-L bond length. Default value is 1.2d0.

* `W-ML-RESTRAINT [real]` Sets the weight of the restraint applied to the Morse metal-ligand bond stretching of term LFMM method. Default value is 1.0d8.

* `REPRODUCE-PRE-2014` Makes Tinker reproduce the behaviour of pre-2014 LFMM calculation, including those from [DommiMOE](https://doi.org/10.1002/jcc.20137). In particular, causes LFSE routines to mix low and high precision cm**-1 to kcal/mol conversion factors (2.859d-3 and 2.8591435d-3).


### LFMM Parameters
The parameters needed for LFMM calculations are defined in an additional section of a normal Tinker parameters file. The following keyword-controlled lines are recognized and interpreted:

* `LFMMBOND [3*integer,4*real]` Polynomial Bond Stretching Parameters. The first two integers correspond to the pair of atom classes (NOT atom types) involved in the bond, while the third is used for the MMFF bond type. The first of the real numbers is the equilibrium bond length, while the other three real numbers are the force constants of the squared, cubic, and quartic term.

* `LFMMANGLE [4*integer,4*real]` Polynomial Angle Bending Parameters. The first three integers correspond to the atom classes (NOT atom types) defining the angle. The last integer defines the MMFF type of the angle. The real numbers define the equilibrium angle, and the force constants of the squared, cubic, and quartic term.

* `LFMMTORS [5*integer,12*real]` Polynomial Torsional Parameters. The first four integers correspond to the atom classes (NOT atom types). The last integer defines the MMFF type of the torsion. Next, each pair of real numbers defines force constant and phase factor for the six term of the polynomial.

* `NONBND1210 [2*integer,2*real]` Parameters for 12-10 Van der Waals Term. The first two integers correspond to the atom classes (NOT atom types) of the non-bonded atoms. The first real is the radius parameter and the second the epsilon.

* `LFMMLLREP [2*integer,real,integer]` Parameters for Pure Repulsive Ligand-Ligand Interaction. The first two integers define the atom classes (NOT atom types) of the metal-ligand pair. Note that the parameters for the coordinating atoms are dependent on the connected metal. The real number is the constant factor and the last integer the exponent.

* `LFMMVWREP [2*integer,3*real]` Parameters for Lennard-Jones Ligand-Ligand Interaction. The first two integers define the atom classes (NOT atom types) of the metal-ligand pair. Note that the parameters for the coordinating atoms are dependent on the connected metal. The real numbers represent the exponent of the repulsive term, the factor for the repulsive term, and the factor for the attractive term.

* `LFMMHARM [2*integer,2*real]` Parameters for Metal-Ligand bond (harmonic). The first two integers correspond to the atom classes (NOT atom types) of metal and ligand atom. Next, the equilibrium bond length and the force constant.

* `LFMMMORSE [2*integer,3*real]` Parameters for Metal-Ligand bond (Morse). The first two integers correspond to the atom classes (NOT atom types) of metal and ligand atom. Next the reference bond length, dissociation energy and curvature.

## Citation
Please, in addition to giving proper credit to Tinker and the LFMM method, remember to cite also <a href="https://doi.org/10.1021/acs.jcim.5b00098">J. Chem. Inf. Model. 2015, 55, 6, 1282–1290</a> when using this modified version of Tinker.

## License
Tinker's licence terms apply. See [doc/license.pdf](doc/license.pdf).
