# Spatially-Resolved-Stochastic-Cluster-Dynamics-SRSCD-simulator-with-Qianran-Yu-
The SRSCD simulation package contains a computer code written in C++. It is an alternative mean field rate theory model that dynamically simulates microstructure (as clusters) evolutions of multi-species systems in solid materials in 1-D space. We have used the code to study Zr-H and W-H system.

****Developed by****

Qianran Yu

Research projects on:
1. Zr-Hydride nucleation and growth under LWR fission condition
2. Hydrogen retention in heavy ion irradiated W materials

Please send an email to the following address for more information:

Qianran Yu (yuqianran0709@gmail.com)

Jaime Marian(jmarian@ucla.edu)

****Introduction****

The SRSCD is a stochastic variant of the mean field rate theory method that is typically developed to simulate microstructure evolutions during permeation of long-term migrating elements, such as hydrogen, into metallic materials. Please see the following journal papers for detail:

[1] Qianran Yu, Michael Reyes, Nachiket Shah and Jaime Marian, "Kinetic Model of Incipient Hydride Formation in Zr Clad under Dynamic Oxide Growth Conditions", Materials 13(5), 1088 (2020). (link: https://www.mdpi.com/1996-1944/13/5/1088)

[2] Qianran Yu, Michael J. Simmonds, Russ. Doerner, George R. Tynan, Li Yang, Brian D. Wirth and Jaime Marian, “Understanding hydrogen retention in damaged tungsten using experimentally-guided models of complex multispecies evolution”, Nuclear Fusion 60, 096003 (2020). (link: https://iopscience.iop.org/article/10.1088/1741-4326/ab9b3c)

****How to use****

Before use: please install gcc/g++ for the newest version (c++11 or newer)

Set up parmeters, check equations in relavant function, and copy/create input files in "src" folder, and then simply type "make". An
executable file named "scdexe" will be generated. Run the simulations by using command "./scdexe".

****Physical mechanisms****

1. (0th) External particle insertion :
   - Once a particle insertion event is selected, a set of random pka energies are selected based on cpdf (obtained from SRIM) until the sum of the energies reaches total incident energy that is expended on lattice damage. Then the code performs collision cascade damage by generating certain number and/or size of defect clusters with statistical variations. Such statistical variation as well as the fraction of SIA/V clusters being generated comes from MD studies.  
2. (1st) Monomer dissociation:
   - single SIA/V atom departs from a cluster.
3. (1st) Defect absorption by dislocations (or sinks)
   - Defects come and get trapped into sinks (such as dislocations, grain boundaries, precipitates etc.)
4. (2nd) Binary combination:
   - Two clusters combine and become a larger cluster, or get annhilated. At least one of the reactants need to be mobile.
5. Long term migration (diffusion):
   - Clusters' long-term migration by Fick's law. 

****Program structure****
A. In the "src" folder:
- Object.cpp/Object.h: store information of species.
- constant.h: store parameters. 
- Damage.cpp/Damage.h: store external particle insertion rates.
- CascadeDamage.cpp/CascadeDamage.h: functions that process cascade damage.
- cpdf.cpp/cpdf.h: sample pka energies using cpdf function.
- rvgs.cpp/rvgs.h: store statistical functions.
- gnuplot_i.h: a library for plotting figures by gnuplot.
- OneLine.cpp/OneLine.h: compuate 1st, 2nd order reaction and diffusion reaction rates. This class handles the information of one species in one spatial element 
- Bundle.cpp/Bundle.h: a class that links information of one species in all spatial elements together.
- SCDWrapper.cpp/SCDWrapper.h: the class that handles the whole rate matrix. It includes ways to select and process events and functions to update rates. There are also output functions. 
- main.cpp: main function
- makefile

B. example_input.zip:
This folder includes the input files to run the example of 3.4 MeV Cu ion irradiation on tungsten materials

****Input files****

The users need to manually create input files to run specific cases. These input files are:
- damage.txt: see "damage.cpp/damage()"
- cpdf#.txt : the number of cpdf files should be equal to the declared number of points (constant.h), see "cpdf.cpp/cpdf()"
- restart.txt : this file is needed when you want to continue the simulation from a previous result, see "SCDWrapper.cpp/restart()"
- sink.txt: store sink absorption information from previous result, see "SCDWrapper.cpp/restart()"

****Additional remarks****

The current code is used to simulate irradiation damage. We also used it to simulate Zr-hydride formation, Hydrogen doposition into W and thermal desorption processes. Each case includes different/additional functions. The user is welcomed to contact Qianran for specific version of the SRSCD code.
