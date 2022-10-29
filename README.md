## Relativistic and Electromagnetic Molecular Dynamics Simulation for Nanoscale Phenomena ##

### Molecular Dynamics Simulation and Necessary Files ###

A molecular dynamics simulation code is implemented for relativistic and electromagnetic fields in three dimensions. It is applied to nanoscale particle phenomena such as nanotube accelerators. Maxwell equations are solved, and momentum equations of relativistic particles are then advanced in time. 
Four physical CGS units are used in this code: a_unit= 1.00d-08 cm, t_unit= 1.00d-15 sec, electron mass m_unit= 0.9110d-27 g and its charge e_unit= 4.8032d-10 esu. The mass of hydrogen, for example, is 1.6726d-24 g. One needs files in the simulation: 1) @cnt3-3p8Ca.f03: Molecular dynamics simulation code, 2) param_em3p8_Ca.h: Common parameters of this simulation, 3) Cntemp_config.STARTC: figure parameters, 4) p_config_ss.xyz_D150 and P135 of pellet electrons, H, C and Au ions. The program is written in Fortran 2003 and MPI of ver.3 for parallelization.

### Courant Condition and Simulation ###

It is noted, however, that the Gauss's equation must be corrected for t>0 time steps as finite errors accumulate in the divergence term. This is true if a discrete coordinate space is used in any method.
Also, all explicit simulation code must satisfy the Courant condition, that is, Dx(length)/Dt(time step) > c, the speed of light. Otherwise, a simulation is overflown shortly. 

A simulation of the nanotube accelerator is set up by putting pellets of H, C and Au atoms and associated electrons at null velocity. Electromagnetic monochromatic waves at the wavelength 800 nm are travelling from the negative direction toward the origin and then go out to the positive direction. The pellets at the origin are irradiated by these waves and are ejected to ion perpenducular and electron parallel directions toward an open space. The final energies for laser intensity 10^22 W/cm^2 are around 30-40 MeV in 20-40 fs, which is shown by animation movies at my homepage.

### Execution Scripts ###

### (1) Linux (PGI): MPI and FFTW by PGI fortan; configure, make, make install.

mpich-4.0.2: ./configure --prefix=/opt/pgi/mpich-4.0.2 2>&1 | tee conf.txt

fftw3-3.3.10: ./configure --disable-shared --enable-maintainer-mode --enable-threads --prefix=/opt/pgi/fftw3

(Old) mpich-3.2: env CC=pgcc FC=pgfortran F77=pgfortran CXX=pgcpp CFLAGS=-fast FCFLAGS=-fast FFLAGS=-fast CXXFLAGS=-fast ./configure --prefix=/opt/pgi/mpich-3.2 --disable-cxx & configure.log

(Old) fftw-3.3.5: env CC=pgcc CFLAGS="-fast -Minfo -fPIC" F77=pgfortran FFLAGS="-fast -Minfo" MPICC=mpicc ./configure --enable-threads --enable-sse2 --enable-openmp --enable-shared --enable-mpi --prefix=/opt/pgi/fftw3

>mpif90 needs param_em3p8_Ca.h, Cntemp_config.STARTC, p_config_ss.xyz_P135 and p_config_ss.xyz_D150.   

>% mpif90 -byteswapio -mcmodel=medium -fast @a_cnt3-3p8Ca.f03 -I/opt/pgi/fftw3/include -L/opt/pgi/fftw3/lib -lfftw3

### (2) Linux (gfortran); configure, make, and make install.

mpich-4.0.2: ./configure --prefix=/opt/mpich-4.0.2 2>&1 | tee conf.txt

fftw3-3.3.10: ./configure --disable-shared --enable-maintainer-mode --enable-threads --prefix=/opt/fftw3

>% mpif90 @a_cnt3-3p8Ca.f03 -I/opt/fftw3/include -L/opt/fftw3/lib -lfftw3

>Execution: Only for a test. % mpiexec -n 6 a.out &


### Simulation Analysis ###

To analyze simulation results, some programs are provided here as the post-processing tool. They are named @3dfdisp.f03 and @3ddisp.f03, for examples. The velocity distributions in parallel and perpendicular directions, @3dfdispC.f03, are plotted in sequential times of ions and electrons. The @3ddisppC.f03 program is time sequential plots of H, C, Au and electrons in side and top views with energy histories as well at the end. These graphic outputs by PDF files are shown on the PC screen, either cntemp.77Cfb.pdf or cntemp.77Csa.pdf. They are discussed in the latter half of the CPC paper in 2019 (Ref. 1 below).

### Parallelization of Fields ###

The heavy load of particles is generally divided into parallel processors which is easily coded.
On the other hand, the electromagnetic fields are parallelized for now the one-dimensional case (the z-direction) where the long axis of pellets is open to eject heavy ions in that direction. Finally, the machine run time depends on the physical time and cpu's architecture. For the elapsed time of parallel simulation it executed in 3.2 sec/step for 52 ranks and 4.0 10^5 particles (equal to 1 fs for elapsed 1.8 hours) by Fujitsu FX100 Supercomputer.

### References: ###

1. M. Tanaka and M. Murakami, Relativistic and electromagnetic molecular dynamics simulations for a carbon-gold nanotube accelerator, Computer Physics Communications, 241, 56-63 (2019).

2. M. Tanaka, A simulation of low-frequency electromagnetic phenomena in kinetic plasmas of three dimensions, J.Comput. Phys., 107, 124-145 (1993).

3. M. Murakami and M. Tanaka, Generation of high-quality mega-electron volt proton beams with intense-laser-driven nanotube accelerator, Applied Phys. Letters, 102, 163101 (2013).

