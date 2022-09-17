Molecular dynamics

A molecular dynamics simulation code has been created for relativistic and 
electromagnetic fields in three dimensions, and is used such as 
nanoscale particle phenomena. 
Maxwell equations are solved, and relativistic particles are then advanced in time.
But, it is noted that the Gauss's equation must be corrected as finite errors accumulate 
in the divergence term. 
This is true if a discrete coordinated scheme is utilized in any method.

Four physical CGS units are used in this code: a_unit= 1.00d-08 cm, 
t_unit= 1.00d-15 sec, electron mass m_unit= 0.9110d-27 g and its charge 
e_unit= 4.8032d-10 esu. 
The mass of hydrogen, for example, is 1.6726d-24 g. 
One needs files in the simulation: 1) @cnt3-3p8Ca.f03: Molecular dynamics 
simulation code, 2) param_em3p8_Ca.h: Common parameters of this simulation, 
3) Cntemp_config.STARTC: figure parameters, 
4) p_config_ss.xyz_D150 and P135: Pellet electrons, H, C and Au ions.

All explicit simulation code must satisfy the Courant condition, that is, 
Dx(length)/Dt(time step) > c, the speed of light. 
Simulations of ejecting nanotube accelerator is discussed in the latter 
half of CPC simulation paper (2019).

References:

1. M. Tanaka and M. Murakami, Relativistic and electromagnetic molecular dynamics simulations 
for a carbon-gold nanotube accelerator, Computer Physics Communications, 241, 56-63 (2019).

2. M. Tanaka, A simulation of low-frequency electromagnetic phenomena in kinetic plasmas of 
three dimensions, J.Comput. Phys., 107, 124-145 (1993).

3. M. Murakami and M. Tanaka, Generation of high-quality mega-electron volt proton beams with 
intense-laser-driven nanotube accelerator, Applied Phys. Letters, 102, 163101 (2013).

