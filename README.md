## Relativistic and Electromagnetic Molecular Dynamics Simulation for Nanoscale Phenomena ##

As "Open Internet Access by Molecular Dynamics Simulations", a couple of various codes are shown in https://github.com/Mtanaka77/, which are "Relativistic and Electromagnetic Molecular Dynamics Simulation for Nanoscale Phenomena", "Large-scale Electromagnetic Particle-in-Cell Simulation", "SIESTA on Vector-Parallel Clusters", and 
"Molecular Dynamics of Water and Ice by TIP5P Code".

This page is discussed on the relativistic and nanoscale molecular dynamics simulations, Computer Physics Communications (2019, Ref. 1). Updated files of @cnt3em_03Aa.f03, parameters and configuration files are uploaded which are dated of Sep. 2024.


### Molecular Dynamics Simulation: CGS Units and Necessary Files ###

A molecular dynamics simulation code is implemented for relativistic and electromagnetic fields 
in three dimensions. It is applied to nanoscale particle phenomena such as nanotube accelerators. 
In the code, Maxwell equations are solved and momentum equations of relativistic particles are advanced in time. 
Four physical CGS units are used in this code: a_unit= 1.00d-08 cm, t_unit= 1.00d-15 sec, 
electron mass m_unit= 0.9110d-27 g and its charge e_unit= 4.8032d-10 esu. 
The mass of hydrogen, for example, is 1.6726d-24 g.

One needs files in the simulation: 1) @cnt3em_03Aa.f03: Molecular dynamics simulation code, 
2) param_em3p7_Aa.h: Common parameters of this simulation, 
3) Cntemp_config.STARTA: configuring parameters, 
4) p_config_ss.xyz_D150 and P135 of pellet electrons: H, C and Au ions and electrons. 
The program is written in Fortran 2003/Fortran 2008 (write format in the same line) and MPI of Ver.3 for parallelization.

The description of each subroutine and important lines of @cnt3em_03Aa.f03 and @3ddisppC.f03 
(to be shown later), is written as comments of the simulation code and post-processing programs. 
Initial 70 lines of the file @cnt3em_03Aa.f03 are devoted to give the title, references, 
summary of subroutines and remarks of the simulation code. 
In the major subroutine /moldyn/, (i) the magnetic field is advanced, (ii) current density is calculated 
and the transverse electric field is advanced, 
(iii) the correction of the longitudinal electric field is made, (iv) the longitudinal electric field is added, 
(v) the forces are calculated, and (vi) positions and momenta of particles are advanced toward the next time step.

### Gauss's Law, Courant Condition and Realistic Mass Simulation ###

It is noted, however, that the Ampere's law becomes inaccurate in
the Cartesian coordinate space, (1/c)\partial{\textbf E}/\partial t=
rot{\textbf B} -(4\pi/c){\textbf J}, residing from the longitudinal
electric field.
Thus, the Gauss's law, div{\textbf E}=4\pi q, must be used for correction
because finite errors accumulate in the divergence term.
That is why the longitudinal electric field must be solved in the discrete coordinate space
(Refs. 1, 2 and 3).
But, the relativistic formulae of velocity and momeutum
/vec{v}= \vec{p}/(sqrt(m^2 +(px^2 +py^2 +pz^2)/c^2)) is valid
in the nanoscale cases (Ref. 1 and Ref. 2).
Also, all the explicit simulation code must satisfy the Courant condition,
that is, Dx(length)/Dt(time step) > c, the speed of light.
Otherwise, a simulation is overflown quite shortly.

A simulation of the nanotube accelerator is set up by putting pellets of H, C and Au atoms 
and associated electrons at null velocity. 
Electromagnetic monochromatic waves at the wavelength 800 nm are travelling from 
the negative direction toward the origin and then go out to the positive direction. 
The pellets at the origin are irradiated by these waves and are ejected to ion perpenducular 
and electron parallel directions toward an open space. 
The final energies for laser intensity 10^22 W/cm^2 are around 30-40 MeV in 20-40 fs, 
which is shown by animation movies at my homepage, http://photon.isc.chubu.ac.jp/.

### Execution Scripts ###

(1) Linux (PGI Fortran): MPI and FFTW; configure, make, make install for installation..

  %  mpich-4.0.2: ./configure --prefix=/opt/pgi/mpich-4.0.2 2>&1 | tee conf.txt

  % fftw3-3.3.10: ./configure --disable-shared --enable-maintainer-mode --enable-threads --prefix=/opt/pgi/fftw3

Compilation by mpif90: 

The script 'mpif90 @a_cnt3-3p8Ca.f03' needs param_em3p8_Ca.h, Cntemp_config.STARTC, p_config_ss.xyz_P135 and p_config_ss.xyz_D150.

  % mpif90 -byteswapio -mcmodel=medium -fast @a_cnt3-3p8Ca.f03 -I/opt/pgi/fftw3/include -L/opt/pgi/fftw3/lib -lfftw3

(2) Linux (gcc gfortran): MPI and FFTW; configure, make, and make install for installation.

  %  mpich-4.0.2: ./configure --prefix=/opt/mpich-4.0.2 2>&1 | tee conf.txt

  % fftw3-3.3.10: ./configure --disable-shared --enable-maintainer-mode --enable-threads --prefix=/opt/fftw3

  % mpif90  -mcmodel=medium @a_cnt3-3p8Ca.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3

Two different fortrans, PGI and gfotran, are incompatible with processors of the same name but different styles.

  % Execution by mpiexec (with lots of co-processors): Only for a test case. % mpiexec -n 6 a.out &


### Post-processing Simulation Analysis ###

To analyze simulation results, this program provides the post-processing tool. 
They are named @3dfdisp.f03 and @3ddisp.f03, for examples. 
The velocity distributions in parallel and perpendicular directions, @3dfdispC.f03, are plotted 
in sequential times of ions and electrons. The @3ddisppC.f03 program is time sequential plots 
of H, C, Au ions and electrons from side and top views with energy histories as well at the end of the run. 
These graphic outputs by PDF files are shown on the PC screen, either cntemp.77Cfb.pdf or 
cntemp.77Csa.pdf. They are discussed in the latter half of the CPC paper in 2019 (Ref. 1 below).

### Parallelization of Fields ###

The heavy load of many particles is generally divided on parallel processors which is easily coded. 
On the other hand, the electromagnetic fields are parallelized for now the one-dimensional case 
(the z-direction) where the long axis of pellets is open to eject heavy ions in that direction. 
Finally, the machine run time depends on the physical time and cpu's architecture. 
For the elapsed time of parallel simulation, it executed in 3.2 sec/step for 52 ranks and 
4.0 10^5 particles (equal to 1 fs for elapsed 1.8 hours) by Fujitsu FX100 Supercomputer.

### References: ###

1. M. Tanaka and M. Murakami, Relativistic and electromagnetic molecular dynamics simulations for a carbon-gold nanotube accelerator, Computer Physics Communications, 241, 56-63 (2019).

2. M. Tanaka, A simulation of low-frequency electromagnetic phenomena in kinetic plasmas of three dimensions, J.Comput. Phys., 107, 124-145 (1993).

3. M. Tanaka, Macro-EM particle simulation method and a study of collisionless magnetic reconnection, Comput.Phys.Commun., 87, 117-138 (1995).

4. M. Murakami and M. Tanaka, Generation of high-quality mega-electron volt proton beams with intense-laser-driven nanotube accelerator, Applied Phys. Letters, 102, 163101 (2013).

