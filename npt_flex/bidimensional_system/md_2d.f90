! simulations of the solid-solid phase transition in 2d model in the NPT ensemble by the Langevin dynamics
! author: Oscar Samuel Cajahuaringa Macollunco and Alex Antonelli
! copyright 2018 Samuel Cajahuaringa <samuelif@ifi.unicamp.br
! article:  Stochastic sampling of the isothermal-isobaric ensemble: phase diagram of crystalline solids 
! from molecular dynamics simulation
! this is free software: you can redistribute it and/or modify it under the terms of 
! the GNU General Public License as published by the Free Software Foundation.
! The code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
! the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.

program md
use routines
implicit none

integer              :: natoms          ! number of atoms
real,    allocatable :: positions(:,:)  ! atomic positions
real,    allocatable :: velocities(:,:) ! velocities
real,    allocatable :: masses(:)       ! masses
real,    allocatable :: forces(:,:)     ! forces   
real                 :: h(2,2)          ! box array in 2d 
real                 :: O(2,2)
real                 :: OT(2,2)
real                 :: vg(2,2)
real                 :: ve1
real                 :: ve2
real                 :: Obox(4,4)
real                 :: OTbox(4,4)
real                 :: h1,h2,h3,h4 
real                 :: tmp(3)
real                 :: tmp0(3)
real                 :: ve(2,2)           ! velocity of barostat in 2d ve[1] = vexx, ve[2] = veyy, ve[3] = vexy = veyx
! neighbour list variables
! see Allen and Tildesey book for details
integer              :: listsize        ! size of the list array
integer, allocatable :: list(:)         ! neighbour list
integer, allocatable :: point(:)        ! pointer to neighbour list
real,    allocatable :: positions0(:,:) ! reference atomic positions, i.e. positions when the neighbour list

! input parameters
! all of them have a reasonable default value, set in read_input()
real           :: tstep          ! simulation timestep
real           :: Tstart         ! temperature start
real           :: Tstop          ! temperature stop
real           :: Tp             ! thermostat period
real           :: temperature    ! temperature
real           :: friction       ! friction for Langevin dynamics (for NVE, use 0)
real           :: Pstart         ! pressure starp
real           :: Pstop          ! pressure stop 
real           :: pressure       ! external pressure (for NPT)
real           :: Tb             ! barostat period
real           :: frictionL      ! friction for batostat 
real           :: wmass          ! mass of barostat
real           :: listcutoff     ! cutoff for neighbour list
real           :: forcecutoff    ! cutoff for forces
integer        :: nstep          ! number of steps
integer        :: nequil
integer        :: nconfig        ! stride for output of configurations
integer        :: nstat          ! stride for output of statistics
integer        :: maxneighbour   ! maximum average number of neighbours per atom
integer        :: idum           ! seed
character(256) :: ensemble
character(256) :: npt_style
character(256) :: inputfile      ! name of file with starting configuration (xyz)
character(256) :: outputfile     ! name of file with final configuration (xyz)
character(256) :: trajfile       ! name of the trajectory file (xyz)
character(256) :: statfile       ! name of the file with statistics
integer :: n1
integer :: n2

integer :: ensemble_flag
integer :: pressure_flag
real    :: engkin         ! kinetic energy
real    :: engkin_tensor(3)
real    :: engconf        ! configurational energy
real    :: engint         ! integral for conserved energy in Langevin dynamics
real    :: area           ! volume of system
real    :: virial(2,2)    ! pressure from virial theorem
real    :: rcm(2)         ! position of center mass of sytem
real    :: pinst(2,2)     ! instantaneous pressure
real    :: pxx,pyy,pxy    
logical :: recompute_list ! control if the neighbour list have to be recomputed
integer :: istep          ! step counter
integer :: iatom
real    :: start,finish
real    :: C0
real    :: C2

call cpu_time(start)

call read_input(Tstart,Tstop,tstep,Tp,Pstart,Pstop,Tb, &
                forcecutoff,listcutoff,nstep,nequil,nconfig,nstat, &
                inputfile,outputfile,trajfile,statfile, &
                maxneighbour,ensemble,npt_style,idum)

if (ensemble == "nve") then
   ensemble_flag = 1
else if (ensemble == "nvt") then
   ensemble_flag = 2
else if (ensemble == "npt") then
   ensemble_flag =3
end if

if (npt_style == "iso") then
   pressure_flag = 1
else if (npt_style == "flex") then
   pressure_flag = 2
end if


! number of atoms is read from file inputfile
call read_natoms(inputfile,natoms)

! particle friction coefficient
friction=0.5/Tp

! barostat friction coefficient
frictionL=0.1/Tb

temperature=Tstart
pressure=Pstart

! barostat mass
if (ensemble_flag == 3) then ! only for npt
   call mass_barostat(wmass,natoms,temperature,Tb,2,pressure_flag)
end if

! write the parameters in output so they can be checked
write(*,*) "Starting configuration           : ",trim(inputfile)
write(*,*) "Final configuration              : ",trim(outputfile)
write(*,*) "Number of atoms                  : ",natoms
write(*,*) "Temperature                      : ",Tstart," ",Tstop
write(*,*) "Time step                        : ",tstep
write(*,*) "Friction                         : ",friction
write(*,*) "Pressure                         : ",Pstart," ",Pstop
write(*,*) "FrictionL                        : ",frictionL
write(*,*) "Mass of barostat                 : ",wmass
write(*,*) "Cutoff for forces                : ",forcecutoff
write(*,*) "Cutoff for neighbour list        : ",listcutoff
write(*,*) "Number of steps                  : ",nstep
write(*,*) "Stride for trajectory            : ",nconfig
write(*,*) "Trajectory file                  : ",trim(trajfile)
write(*,*) "Stride for statistics            : ",nstat
write(*,*) "Statistics file                  : ",trim(statfile)
write(*,*) "Max average number of neighbours : ",maxneighbour
write(*,*) "Seed                             : ",idum
!write(*,*) "Are atoms wrapped on output?     : ",wrapatoms

! Since each atom pair is counted once, the total number of pairs
! will be half of the number of neighbours times the number of atoms
listsize=maxneighbour*natoms/2

! allocation of dynamical arrays
allocate(positions(2,natoms))
allocate(positions0(2,natoms))
allocate(velocities(2,natoms))
allocate(forces(2,natoms))
allocate(masses(natoms))
allocate(point(natoms))
allocate(list(listsize))

! masses are hard-coded to 1
masses=1.0

! energy integral initialized to 0
engint=0.0

! positions are read from file inputfile
call read_positions(inputfile,natoms,positions,h)

!call square(1.08,n1,n2,natoms,positions,h)
!call triangular(1.464,n1,n2,natoms,positions,h)

! velocities are randomized according to temperature
call randomize_velocities(natoms,temperature,masses,velocities,idum)
call compute_remove_vcm(natoms,masses,velocities)

! neighbour list are computed, and reference positions are saved
call compute_list(natoms,listsize,positions,h,listcutoff,point,list)
!write(*,*) "List size: ",point(natoms)-1
!positions0=positions

call coeff(1.0,forcecutoff,c0,C2)

! forces are computed before starting md
call compute_forces(natoms,listsize,positions,h,c0,c2,forcecutoff,point,list,forces,engconf,virial)

area=h(1,1)*h(2,2)-h(2,1)*h(1,2)   

! velocity of barostat
ve=0.0

do istep=1,nequil

  if (ensemble_flag == 3) then ! only for npt
     if (pressure_flag == 1) then
        call thermostat_barostat_iso(wmass,0.5*tstep,frictionL,temperature,ve(1,1),engint,idum)
     else if (pressure_flag == 2) then
        call thermostat_barostat_flex(wmass,0.5*tstep,frictionL,temperature,ve,engint,idum)  
     end if
  end if

  if (ensemble_flag /= 1) then ! for npt and nvt
     call thermostat_particles(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,idum)
     call compute_remove_vcm(natoms,masses,velocities)
  end if
  
  if (ensemble_flag == 3) then ! only for npt 
     call compute_engkin(natoms,masses,velocities,engkin,engkin_tensor)
     call compute_pressure(pinst,area,virial,engkin_tensor)
     if (pressure_flag == 1) then
        call velocity_barostat_iso(ve(1,1),area,pinst,pressure,natoms,engkin,tstep,wmass)
        call velocities_particles_iso(natoms,velocities,masses,forces,ve(1,1),tstep)
        call positions_particles_iso(natoms,positions,velocities,ve(1,1),tstep)
        call box_propagated_iso(h,ve(1,1),tstep)  
     else if (pressure_flag == 2) then
        call velocity_barostat_flex(ve,area,pinst,pressure,natoms,engkin,tstep,wmass)
        call orthogonal_array(ve,O,OT,ve1,ve2)
        call velocities_particles_flex(natoms,velocities,masses,forces,O,OT,ve1,ve2,tstep)   
        call positions_particles_flex(natoms,positions,velocities,O,OT,ve1,ve2,tstep)
        call orthogonal_box_array(ve,Obox,OTbox,h1,h2,h3,h4)
        call box_propagated_flex(h,Obox,OTbox,h1,h2,h3,h4,tstep)
     end if
  else  ! for nvt and nve 
     call compute_verlet_velocities(natoms,velocities,masses,forces,tstep)
     call compute_verlet_positions(natoms,positions,velocities,tstep) 
  end if

  area=h(1,1)*h(2,2)-h(2,1)*h(1,2)
 
  call compute_list(natoms,listsize,positions,h,listcutoff,point,list)

  call compute_forces(natoms,listsize,positions,h,c0,c2,forcecutoff,point,list,forces,engconf,virial)

  if (ensemble_flag == 3) then
     if (pressure_flag == 1) then
        call velocities_particles_iso(natoms,velocities,masses,forces,ve(1,1),tstep)
     else if (pressure_flag == 2) then
        call velocities_particles_flex(natoms,velocities,masses,forces,O,OT,ve1,ve2,tstep)
     end if
     call compute_engkin(natoms,masses,velocities,engkin,engkin_tensor)
     call compute_pressure(pinst,area,virial,engkin_tensor)
     if (pressure_flag == 1) then
        call velocity_barostat_iso(ve(1,1),area,pinst,pressure,natoms,engkin,tstep,wmass)
     else if (pressure_flag == 2) then
        call velocity_barostat_flex(ve,area,pinst,pressure,natoms,engkin,tstep,wmass)
     end if
  else                              
     call compute_verlet_velocities(natoms,velocities,masses,forces,tstep)
  end if

  if (ensemble_flag /= 1) then ! for npt and nvt
     call thermostat_particles(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,idum)
     call compute_remove_vcm(natoms,masses,velocities)
  end if

  if (ensemble_flag == 3) then ! only for npt
     if (pressure_flag == 1) then
        call thermostat_barostat_iso(wmass,0.5*tstep,frictionL,temperature,ve(1,1),engint,idum)
     else if (pressure_flag == 2) then
        call thermostat_barostat_flex(wmass,0.5*tstep,frictionL,temperature,ve,engint,idum)
     end if
  end if
 
  ! kinetic energy is calculated
  !call compute_engkin(natoms,masses,velocities,engkin,engkin_tensor)

  ! instantaneus pressure
  !call compute_pressure(pinst,area,virial,engkin_tensor)

end do

! energy integral initialized to 0
engint=0.0

! here is the main md loop
! Langevin thermostat is applied before and after a velocity-Verlet integrator
! the overall structure is:
!   thermostat
!   update velocities
!   update positions
!   (eventually recompute neighbour list)
!   compute forces
!   update velocities
!   thermostat
!   (eventually dump output informations)
do istep=0,nstep

  if (ensemble_flag /= 1) then
     call temperature_system(temperature,Tstart,Tstop,istep,nstep)
  end if

  if (ensemble_flag == 3) then ! only for npt
     call pressure_system(pressure,Pstart,Pstop,istep,nstep) 
     call mass_barostat(wmass,natoms,temperature,Tb,2,pressure_flag) 
     if (pressure_flag == 1) then
        call thermostat_barostat_iso(wmass,0.5*tstep,frictionL,temperature,ve(1,1),engint,idum)
     else if (pressure_flag == 2) then
        call thermostat_barostat_flex(wmass,0.5*tstep,frictionL,temperature,ve,engint,idum)
     end if
  end if

  if (ensemble_flag /= 1) then ! for npt and nvt
     call thermostat_particles(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,idum)
     call compute_remove_vcm(natoms,masses,velocities)
  end if

  if (ensemble_flag == 3) then ! only for npt 
     call compute_engkin(natoms,masses,velocities,engkin,engkin_tensor)
     call compute_pressure(pinst,area,virial,engkin_tensor)
     if (pressure_flag == 1) then
        call velocity_barostat_iso(ve(1,1),area,pinst,pressure,natoms,engkin,tstep,wmass)
        call velocities_particles_iso(natoms,velocities,masses,forces,ve(1,1),tstep)
        call positions_particles_iso(natoms,positions,velocities,ve(1,1),tstep)
        call box_propagated_iso(h,ve(1,1),tstep)
     else if (pressure_flag == 2) then
        call velocity_barostat_flex(ve,area,pinst,pressure,natoms,engkin,tstep,wmass)
        call orthogonal_array(ve,O,OT,ve1,ve2)
        call velocities_particles_flex(natoms,velocities,masses,forces,O,OT,ve1,ve2,tstep)
        call positions_particles_flex(natoms,positions,velocities,O,OT,ve1,ve2,tstep)
        call orthogonal_box_array(ve,Obox,OTbox,h1,h2,h3,h4)
        call box_propagated_flex(h,Obox,OTbox,h1,h2,h3,h4,tstep)
     end if
  else  ! for nvt and nve 
     call compute_verlet_velocities(natoms,velocities,masses,forces,tstep)
     call compute_verlet_positions(natoms,positions,velocities,tstep)
  end if

  area=h(1,1)*h(2,2)-h(2,1)*h(1,2)

  call compute_list(natoms,listsize,positions,h,listcutoff,point,list)

  call compute_forces(natoms,listsize,positions,h,c0,c2,forcecutoff,point,list,forces,engconf,virial)

  if (ensemble_flag == 3) then
     if (pressure_flag == 1) then
        call velocities_particles_iso(natoms,velocities,masses,forces,ve(1,1),tstep)
     else if (pressure_flag == 2) then
        call velocities_particles_flex(natoms,velocities,masses,forces,O,OT,ve1,ve2,tstep)
     end if
     call compute_engkin(natoms,masses,velocities,engkin,engkin_tensor)
     call compute_pressure(pinst,area,virial,engkin_tensor)
     if (pressure_flag == 1) then
        call velocity_barostat_iso(ve(1,1),area,pinst,pressure,natoms,engkin,tstep,wmass)
     else if (pressure_flag == 2) then
        call velocity_barostat_flex(ve,area,pinst,pressure,natoms,engkin,tstep,wmass)
     end if
  else
     call compute_verlet_velocities(natoms,velocities,masses,forces,tstep)
  end if

  if (ensemble_flag /= 1) then ! for npt and nvt
     call thermostat_particles(natoms,masses,0.5*tstep,friction,temperature,velocities,engint,idum)
     call compute_remove_vcm(natoms,masses,velocities)
  end if

  if (ensemble_flag == 3) then ! only for npt
     if (pressure_flag == 1) then
        call thermostat_barostat_iso(wmass,0.5*tstep,frictionL,temperature,ve(1,1),engint,idum)
     else if (pressure_flag == 2) then
        call thermostat_barostat_flex(wmass,0.5*tstep,frictionL,temperature,ve,engint,idum)
     end if
  end if

  ! kinetic energy is calculated
  call compute_remove_vcm(natoms,masses,velocities)
  call compute_engkin(natoms,masses,velocities,engkin,engkin_tensor)

  ! instantaneus pressure
  call compute_pressure(pinst,area,virial,engkin_tensor)

  !call compute_center_mass(natoms,masses,positions,rcm)

  if(modulo(istep,nconfig)==0) call write_positions(trajfile,natoms,positions,h,istep)
  if(modulo(istep,nstat)==0)   call write_statistics(statfile,istep,tstep,natoms,engkin,engconf,pressure,pinst,h,wmass,ve,engint)
!  if(modulo(istep,nstat)==0) write(6,*) istep,rcm(1),rcm(2)
end do

call write_final_positions(outputfile,natoms,positions,h)
write(*,*) "Execution completed"
call cpu_time(finish)
write(*,*) "Time =",finish-start,"seconds"


! deallocation of all allocatable array
deallocate(positions)
deallocate(velocities)
deallocate(forces)
deallocate(masses)
deallocate(positions0)
deallocate(point)
deallocate(list)

end program md
