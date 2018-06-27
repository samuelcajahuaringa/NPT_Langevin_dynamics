! simulations of the jensen model in the NPT ensemble by the Langevin dynamics
! author: Oscar Samuel Cajahuaringa Macollunco and Alex Antonelli
! copyright 2018 Samuel Cajahuaringa <samuelif@ifi.unicamp.br
! article: Stochastic sampling of the isothermal-isobaric ensemble: phase diagram of crystalline solids 
! this is free software: you can redistribute it and/or modify it under the terms of 
! the GNU General Public License as published by the Free Software Foundation.
! The code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
! the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.

program npt
use routines
implicit none

integer              :: natoms        ! number of atoms
real,    allocatable :: positions(:)  ! atomic positions
real,    allocatable :: velocities(:) ! velocities
                                      ! was calculated last time
real,    allocatable :: masses(:)     ! masses
real,    allocatable :: forces(:)     ! forces   
real                 :: cell          ! cell size
real                 :: cellnew
real                 :: cellold
real                 :: ve            ! velocity barostat
real                 :: pos_scale
real                 :: pos_poly
real                 :: vel_scale
real                 :: vel_poly
    
! input parameters
! all of them have a reasonable default value, set in read_input()
real           :: tstep          ! simulation timestep
real           :: temperature    ! temperature in kB units
real           :: friction       ! friction for Langevin dynamics (for NVE, use 0)
real           :: pressure       ! external pressure for NPT dynamics
real           :: frictionL      ! friction for Langevin dynamics of barostat
real           :: wmass          ! barostat mass
integer        :: nequil         ! number of steps of equilibration
integer        :: nstep          ! number of steps of thermalization
integer        :: nconfig        ! stride for output of configurations
integer        :: nstat          ! stride for output of statistics
integer        :: idum           ! seed
real           :: deep           ! potential deep
logical        :: wrapatoms      ! if true, atomic coordinates are written wrapped in minimal cell
character(256) :: inputfile      ! name of file with starting configuration (xyz)
character(256) :: outputfile     ! name of file with final configuration (xyz)
character(256) :: trajfile       ! name of the trajectory file (xyz)
character(256) :: statfile       ! name of the file with statistics
character(256) :: cellfile       ! name of the file with the length

real    :: engkin         ! kinetic energy
real    :: engconf        ! configurational energy
real    :: engkinbar      ! kinetic energy of barostat
real    :: engint         ! integral for conserved energy in Langevin dynamics
real    :: lbox           ! length of box simulation
real    :: virial         ! virial part from force 
real    :: pint           ! instantaneous pressure
real    :: alpha           
integer :: istep          ! step counter
integer :: iatom
real    :: Hi             ! effective energy    
real    :: H0             ! initial effective energy
real    :: dH
real    :: drift
real    :: taub
real    :: Tstart,Tstop
real    :: Pstart,Pstop


call read_input(Tstart,Tstop,tstep,friction,Pstart,Pstop,frictionL,nstep,nequil,nconfig,nstat, &
                inputfile,cellfile,outputfile,trajfile,statfile,idum,deep)

! number of atoms is read from file inputfile
call read_natoms(inputfile,natoms)

! mass of barostat
taub=0.1/frictionL
temperature=Tstart
pressure=Pstart
call mass_barostat(wmass,natoms,temperature,Taub)

! write the parameters in output so they can be checked
write(*,*) "Starting configuration           : ",trim(inputfile)
write(*,*) "Final configuration              : ",trim(outputfile)
write(*,*) "Number of atoms                  : ",natoms
write(*,*) "Temperature                      : ",Tstart,Tstop
write(*,*) "Time step                        : ",tstep
write(*,*) "Friction                         : ",friction
write(*,*) "Pressure                         : ",Pstart,Pstop
write(*,*) "FrictionL                        : ",frictionL
write(*,*) "Mass of barostat                 : ",wmass
write(*,*) "Number of steps                  : ",nstep
write(*,*) "Number of steps of equilibration : ",nequil
write(*,*) "Stride for trajectory            : ",nconfig
write(*,*) "Trajectory file                  : ",trim(trajfile)
write(*,*) "Stride for statistics            : ",nstat
write(*,*) "Statistics file                  : ",trim(statfile)
write(*,*) "Seed                             : ",idum
write(*,*) "Deep of potential energy         : ",deep
write(*,*) "Cell file                        : ",cellfile

! allocation of dynamical arrays
allocate(positions(natoms))
allocate(velocities(natoms))
allocate(forces(natoms))
allocate(masses(natoms))
! masses are hard-coded to 1
masses=1.0

alpha=1.0+1.0/natoms

ve=0.0

! positions are read from file inputfile
call read_positions(inputfile,natoms,positions,cell)

call boundary_conditions(natoms,cell,positions)

! velocities are randomized according to temperature
call randomize_velocities(natoms,temperature,masses,velocities,idum)

call compute_remove_vcm(natoms,masses,velocities)

! forces are computed before starting md
call compute_forces(natoms,positions,cell,forces,engconf,deep,virial)

call compute_engkin(natoms,masses,velocities,engkin)

pint=(2.0*engkin+virial)/cell

do istep=1,nstep

  call thermostat(natoms,masses,wmass,0.5*tstep,friction,frictionL,temperature,velocities,ve,engint,idum)

  call compute_remove_vcm(natoms,masses,velocities)

  call compute_engkin(natoms,masses,velocities,engkin)

  pint=(2.0*engkin+virial)/cell

  ve=ve+(cell*(pint-pressure)+2.0*engkin/natoms)*0.5*tstep/wmass

! for positions
  pos_scale=exp(ve*0.5*tstep)
  call sinhx((ve*0.5*tstep)**2,pos_poly)

! for velocities
  vel_scale=exp(-alpha*ve*0.25*tstep)
  call sinhx((-alpha*ve*0.25*tstep)**2,vel_poly)

! velocities at dt/2
  do iatom=1,natoms
    velocities(iatom)=vel_scale*(vel_scale*velocities(iatom)+vel_poly*forces(iatom)*0.5*tstep/masses(iatom))
  end do  

! position at dt
  do iatom=1,natoms
    positions(iatom)=pos_scale*(pos_scale*positions(iatom)+pos_poly*velocities(iatom)*tstep)
  end do
  
  cell=exp(ve*tstep)*cell

  call compute_forces(natoms,positions,cell,forces,engconf,deep,virial)

! velocities at dt/2
  do iatom=1,natoms
    velocities(iatom)=vel_scale*(vel_scale*velocities(iatom)+vel_poly*forces(iatom)*0.5*tstep/masses(iatom))
  end do

  call compute_engkin(natoms,masses,velocities,engkin)

  pint=(2.0*engkin+virial)/cell

  ve=ve+(cell*(pint-pressure)+2.0*engkin/natoms)*0.5*tstep/wmass
  
  call thermostat(natoms,masses,wmass,0.5*tstep,friction,frictionL,temperature,velocities,ve,engint,idum)

  call compute_remove_vcm(natoms,masses,velocities)
! kinetic energy is calculated
 call compute_engkin(natoms,masses,velocities,engkin)

! instantaneus pressure
  pint=(2.0*engkin+virial)/cell

end do

engint=0.0
engkinbar=0.5*wmass*ve**2
H0=engkin+engconf+cell*pressure+engkinbar+engint

do istep=1,nequil

  if (Tstart /= Tstop) then 
     call temperature_system(temperature,Tstart,Tstop,istep,nequil)
     call mass_barostat(wmass,natoms,temperature,Taub)
  end if

  if (Pstart /= Pstop) call pressure_system(pressure,Pstart,Pstop,istep,nequil)

  call thermostat(natoms,masses,wmass,0.5*tstep,friction,frictionL,temperature,velocities,ve,engint,idum)

  call compute_remove_vcm(natoms,masses,velocities)

  call compute_engkin(natoms,masses,velocities,engkin)

  pint=(2.0*engkin+virial)/cell

  ve=ve+(cell*(pint-pressure)+2.0*engkin/natoms)*0.5*tstep/wmass

! for positions
  pos_scale=exp(ve*0.5*tstep)
  call sinhx((ve*0.5*tstep)**2,pos_poly)

! for velocities
  vel_scale=exp(-alpha*ve*0.25*tstep)
  call sinhx((-alpha*ve*0.25*tstep)**2,vel_poly)

! velocities at dt/2
  do iatom=1,natoms
    velocities(iatom)=vel_scale*(vel_scale*velocities(iatom)+vel_poly*forces(iatom)*0.5*tstep/masses(iatom))
  end do

! position at dt
  do iatom=1,natoms
    positions(iatom)=pos_scale*(pos_scale*positions(iatom)+pos_poly*velocities(iatom)*tstep)
  end do

  cell=exp(ve*tstep)*cell

  call compute_forces(natoms,positions,cell,forces,engconf,deep,virial)

! velocities at dt/2
  do iatom=1,natoms
    velocities(iatom)=vel_scale*(vel_scale*velocities(iatom)+vel_poly*forces(iatom)*0.5*tstep/masses(iatom))
  end do

  call compute_engkin(natoms,masses,velocities,engkin)

  pint=(2.0*engkin+virial)/cell

  ve=ve+(cell*(pint-pressure)+2.0*engkin/natoms)*0.5*tstep/wmass

  call thermostat(natoms,masses,wmass,0.5*tstep,friction,frictionL,temperature,velocities,ve,engint,idum)

  call compute_remove_vcm(natoms,masses,velocities)

! kinetic energy is calculated
  call compute_engkin(natoms,masses,velocities,engkin)

! instantaneus pressure
  pint=(2.0*engkin+virial)/cell
 
  engkinbar=0.5*wmass*ve**2

! eventually, write positions and statistics
  if(modulo(istep,nconfig)==0) call write_positions(trajfile,natoms,positions,cell,istep)
  if(modulo(istep,nstat)==0) call write_statistics(statfile,istep,tstep,natoms,engkin,engconf,engint,pint,cell,&
  engkinbar,ve,pressure)

  if (Pstart == Pstop) then
     call write_cell_only(cellfile,cell,istep)
  else 
     call write_cell(cellfile,cell,pressure,istep)
  end if
  
end do
  
call write_final_positions(outputfile,natoms,positions,cell)
write(*,*) "Execution completed"

! deallocation of all allocatable array
deallocate(positions)
deallocate(velocities)
deallocate(forces)
deallocate(masses)

end program npt
