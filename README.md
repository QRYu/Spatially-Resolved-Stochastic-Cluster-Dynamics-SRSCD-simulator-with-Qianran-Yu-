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
4. (2st) Binary combination:
   - Two clusters combine and become a larger cluster, or get annhilated. At least one of the reactants need to be mobile.
5. Long term migration (diffusion):
   - Clusters' long-term migration by Fick's law. 

****Program structure****
- Object.cpp/Object.h: 
