
molecular_dynamics

A molecular dynamics simulation code has been created 
for applications of relativistic and electromagnetic fields 
in three dimensions, as for nanoscale phenomena. 
Maxwell equations are solved with Gauss's law for necessary 
corrections. Relativistic particles are the advanced in time.

The CGS units are used: a_unit= 1.0000d+00 cm, m_unit= 0.9110d-27 g, 
electron mass e_unit= 4.8032d-10 esu, and t_unit= 1.0000d+00 sec.
One expects four files: 1) @cnt3-3p8Ca.f03: Molecular dynamics 
simulation code, 2) param_em3p8_Ca.h: Common parameters of 
this simulation, 3) Cntemp_config.STARTC: Configure parameters,  
4) p_config_ss.xyz_D150, P135: Pellets electrons and ions.

All explicit simulation code must satisfy the Courant condition, 
that is, Dx(length)/Dt(time step) > c, the speed of light. 
Simulations of ejecting nanotube accelerator is discussed 
in the later half of CPC 2019 paper.

Ref: Molecular dynamics, M. Tanaka, Computer Physics Communicatios, 
vol.241, 56 (2019); Implicit particle simulation code, M. Tanaka, 
Comput.Phys.Comm., vol.87, 117 (1995); M.Tanaka, J.Comput.Phys., 
vol.107, 124 (1993).

