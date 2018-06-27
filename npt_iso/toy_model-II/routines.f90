module routines
implicit none
contains
subroutine read_input(Tstart,Tstop,tstep,friction,Pstart,Pstop,frictionL,nstep,nequil,nconfig,nstat, &
                      inputfile,cellfile,outputfile,trajfile,statfile,idum,deep)
  implicit none
  real,           intent(out) :: Tstart,Tstop
  real,           intent(out) :: tstep
  real,           intent(out) :: friction
  real,           intent(out) :: Pstart,Pstop
  real,           intent(out) :: frictionL
  integer,        intent(out) :: nstep
  integer,        intent(out) :: nequil
  integer,        intent(out) :: nconfig
  integer,        intent(out) :: nstat
  character(256), intent(out) :: inputfile
  character(256), intent(out) :: cellfile
  character(256), intent(out) :: outputfile
  character(256), intent(out) :: trajfile
  character(256), intent(out) :: statfile
  integer,        intent(out) :: idum
  real,           intent(out) :: deep

  integer :: iostat
  character(256) :: line,keyword,keyword1
  integer :: i
  logical :: foundsharp
! default values
  Tstart=1.0
  Tstop=1.0
  tstep=0.005
  friction=0.0
  Pstart=0.0
  Pstop=0.0
  frictionL=0.0
  nstep=1
  nconfig=10
  nstat=1
  idum=0
  deep=1.0
  statfile=""
  trajfile=""
  outputfile=""
  inputfile=""
  cellfile=""
  nequil=1

  do
    read(*,"(a)",iostat=iostat) line
! when the file finishes, exit the loop
    if(iostat/=0) exit
! delete everything past an eventual "#" comment
    foundsharp=.false.
    do i=1,len(line)
      if(line(i:i)=="#") foundsharp=.true.
      if(foundsharp) line(i:i)=" "
    end do
! if the remaining line is empty, skip it
    if(len_trim(line)==0) cycle
! read the first word from line
    read(line,*) keyword
! the second word is then read to the proper variable
    select case(keyword)
    case("temperature")
      read(line,*) keyword1,Tstart,Tstop
    case("tstep")
      read(line,*) keyword1,tstep
    case("friction")
      read(line,*) keyword1,friction
    case("pressure")
      read(line,*) keyword1,Pstart,Pstop
    case("frictionL")
      read(line,*) keyword1,frictionL
    case("nstep")
      read(line,*) keyword1,nstep
    case("nequil")
      read(line,*) keyword1,nequil
    case("nconfig")
      read(line,*) keyword1,nconfig,trajfile
    case("nstat")
      read(line,*) keyword1,nstat,statfile
    case("inputfile")
      read(line,*) keyword1,inputfile
    case("outputfile")
      read(line,*) keyword1,outputfile
    case("cellfile")
      read(line,*) keyword1,cellfile
    case("deep")
      read(line,*) keyword1,deep
    case("seed")
      read(line,*) keyword1,idum
      idum=-idum ! idum for ran1() needs to be negative
    case default
! an unknown word will stop the execution
      write(0,*) "Unknown keyword :",trim(keyword)
      stop
    end select
  end do
  if(inputfile=="") then
    write(0,*) "Specify input file"
    stop
  end if
  if(outputfile=="") then
    write(0,*) "Specify output file"
    stop
  end if
  if(trajfile=="") then
    write(0,*) "Specify traj file"
    stop
  end if
  if(statfile=="") then
    write(0,*) "Specify stat file"
    stop
  end if
end subroutine read_input

subroutine read_positions(inputfile,natoms,positions,cell)
! read positions and cell from a file called inputfile
! natoms (input variable) and number of atoms in the file should be consistent
  implicit none
  character(*), intent(in)  :: inputfile
  integer,      intent(in)  :: natoms
  real,         intent(out) :: positions(natoms)
  real,         intent(out) :: cell
  integer :: iatom
  open(10,file=inputfile)
  read(10,*)
  read(10,*) cell
  do iatom=1,natoms
    read(10,*) positions(iatom)
  end do
! note: atomname is read but not used
  close(10)
end subroutine read_positions

subroutine read_natoms(inputfile,natoms)
! read the number of atoms in file "input.xyz"
  implicit none
  character(*),     intent(in) :: inputfile
  integer,          intent(out):: natoms
  open(10,file=inputfile)
  read(10,*) natoms
  close(10)
end subroutine read_natoms

subroutine write_positions(trajfile,natoms,positions,cell,istep)
! write positions on file trajfile
! positions are appended at the end of the file
  implicit none
  character(*), intent(in) :: trajfile
  integer,      intent(in) :: natoms
  real,         intent(in) :: positions(natoms)
  real,         intent(in) :: cell
  integer,      intent(in) :: istep
  integer :: iatom
  real :: pos(natoms)
  logical, save :: first=.true.
  if(first) then
    open(10,file=trajfile)
    first=.false.
  else
    open(10,file=trajfile,position="append")
  end if
! lammps format of trajectory file
1000 FORMAT (A14)
1001 FORMAT (I6)
1002 FORMAT (A21)
1003 FORMAT (I4)
1004 FORMAT (A25)
1005 FORMAT (F10.4,2X,F10.4)
1006 FORMAT (A20)
1007 FORMAT (I4,1X,F10.4,1X,F10.4,1X,F10.4)
  write(10,1000)"ITEM: TIMESTEP"
  write(10,1001)istep
  write(10,1002)"ITEM: NUMBER OF ATOMS"
  write(10,1003)natoms 
  write(10,1004)"ITEM: BOX BOUNDS pp ff ff"
  write(10,1005)0.0,cell
  write(10,1005)-0.25,0.25
  write(10,1005)-0.25,0.25
  write(10,1006)"ITEM: ATOMS id x y z"   
  
  pos=positions


call boundary_conditions(natoms,cell,pos)

  do iatom=1,natoms
! usually, it is better not to apply pbc here, so that diffusion
! is more easily calculated from a trajectory file:
!    if(wrapatoms) then
!      call pbc(cell,positions(:,iatom),pos)
!    else
    
!    end if
  write(10,1007)iatom,pos(iatom),0.0,0.0
  end do
  close(10)
end subroutine write_positions

subroutine write_final_positions(outputfile,natoms,positions,cell)
! write positions on file outputfile
  character(*), intent(in) :: outputfile
  integer,      intent(in) :: natoms
  real,         intent(in) :: positions(natoms)
  real,         intent(in) :: cell
  real :: pos(natoms)
  integer :: iatom

  pos=positions

call boundary_conditions(natoms,cell,pos)

  open(10,file=outputfile)
  write(10,*) natoms
  write(10,*) cell
  do iatom=1,natoms
! usually, it is better not to apply pbc here, so that diffusion
! is more easily calculated from a trajectory file:
!    if(wrapatoms) then
!      call pbc(cell,positions(:,iatom),pos)
!    else
!      pos(iatom)
!    end if
    write(10,*) pos(iatom)
  end do
  close(10)
end subroutine write_final_positions

subroutine write_statistics(statfile,istep,tstep,natoms,engkin,engconf,engint,pint,cell,engkinbar,ve,pressure)
! write statistics on file statfile
  character(*), intent(in) :: statfile
  integer,      intent(in) :: istep
  real,         intent(in) :: tstep
  integer,      intent(in) :: natoms
  real,         intent(in) :: engkin
  real,         intent(in) :: engconf
  real,         intent(in) :: engint
  real,         intent(in) :: pint
  real,         intent(in) :: cell
  real,         intent(in) :: ve
  real,         intent(in) :: pressure
  real,         intent(in) :: engkinbar
  logical, save :: first=.true.
  integer, save :: last_time_reopened=0
  real :: energy
  real :: enthalpy
  real :: nonhamiltonian
  real :: drift

  energy = engkin+engconf
  enthalpy = energy+cell*pressure
  nonhamiltonian = enthalpy+engkinbar
 
  if(first) then
! first time this routine is called, open the file
    open(666,file=statfile)
    first=.false.
  end if
  if(istep-last_time_reopened>100) then
! every 100 steps, reopen the file to flush the buffer
    close(666)
    open(666,file=statfile,position="append")
    last_time_reopened=istep
  end if
  write(666,"(i10,99g15.7)") istep,istep*tstep,2.0*engkin/natoms,engconf,energy,pint,cell,enthalpy,ve,engkinbar,nonhamiltonian
end subroutine write_statistics

subroutine write_cell_only(cellfile,cell,istep)
  character(*), intent(in) :: cellfile
  real,         intent(in) :: cell
  integer,      intent(in) :: istep
  real :: H
  logical, save :: first=.true.
  integer, save :: last_time_reopened=0
 
  if(first) then
! first time this routine is called, open the file
    open(777,file=cellfile)
    first=.false.
  end if
  if(istep-last_time_reopened>100) then
! every 100 steps, reopen the file to flush the buffer
    close(777)
    open(777,file=cellfile,position="append")
    last_time_reopened=istep
  end if
  write(777,"(F15.7)") cell  
end subroutine write_cell_only

subroutine write_cell(cellfile,cell,pressure,istep)
  character(*), intent(in) :: cellfile
  real,         intent(in) :: cell
  integer,      intent(in) :: istep
  real,         intent(in) :: pressure
  real :: H
  logical, save :: first=.true.
  integer, save :: last_time_reopened=0

  if(first) then
! first time this routine is called, open the file
    open(777,file=cellfile)
    first=.false.
  end if
  if(istep-last_time_reopened>100) then
! every 100 steps, reopen the file to flush the buffer
    close(777)
    open(777,file=cellfile,position="append")
    last_time_reopened=istep
  end if
1000 format (F9.6,2X,F8.3)

  write(777,1000)pressure,cell
end subroutine write_cell

subroutine randomize_velocities(natoms,temperature,masses,velocities,idum)
! randomize the velocities according to the temperature
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: temperature
  real,    intent(in)    :: masses(natoms)
  real,    intent(out)   :: velocities(natoms)
  integer, intent(inout) :: idum
  real, external :: gasdev
  integer :: iatom
  do iatom=1,natoms
      velocities(iatom)=sqrt(temperature/masses(iatom))*gasdev(idum)
  end do
end subroutine randomize_velocities

subroutine compute_engkin(natoms,masses,velocities,engkin)
! calculate the kinetic energy from the velocities
  implicit none
  integer, intent(in) :: natoms
  real,    intent(in) :: masses(natoms)
  real,    intent(in) :: velocities(natoms)
  real,    intent(out):: engkin
  integer :: iatom
  engkin=0.0
  do iatom=1,natoms
    engkin=engkin+0.5*masses(iatom)*velocities(iatom)**2
  end do
end subroutine compute_engkin

subroutine boundary_conditions(natoms,cell,positions)
! apply periodice boundary condition for the positions inside the box
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: cell
  real,    intent(inout) :: positions(natoms)
  real :: xi
  integer :: i
  do i=1,natoms
     xi = positions(i)
     if (xi > cell) then
     xi = xi - cell
     else if (xi < 0.0) then
     xi = xi + cell
     end if
     positions(i) = xi
  end do  
end subroutine boundary_conditions

subroutine pbc(cell,vin,vout)
! apply periodic boundary condition to a vector
  implicit none
  real, intent(in) :: cell
  real, intent(in) :: vin
  real, intent(out) :: vout
    vout=vin-nint(vin/cell)*cell
end subroutine pbc

subroutine compute_list(natoms,list)
  implicit none
  integer, intent(in)  :: natoms
  integer, intent(out) :: list(2,natoms)
  integer :: i,j
  list(1,1)=natoms
  list(2,1)=3
  do i=2,natoms-1
     list(1,i)=i-1
     list(2,i)=i+1
  end do
  list(1,natoms)=natoms-1
  list(2,natoms)=1   
end subroutine compute_list

subroutine compute_forces(natoms,positions,cell,forces,engconf,deep,virial)
  implicit none
  integer, intent(in)  :: natoms
  real,    intent(in)  :: positions(natoms)
  real,    intent(in)  :: cell
!  integer, intent(in)  :: list(2,natoms)
  real,    intent(out) :: forces(natoms)
  real,    intent(out) :: engconf
  real,    intent(in)  :: deep
  real,    intent(out) :: virial ! virial part from the force  of instantaneous pressure 
  integer :: iatom,jatom  ! indexes of the two involved atoms
  integer :: jlist        ! counter for the neighbours of iatom
  real :: dist            ! distance of the two atoms
  real :: rij             ! minimum-image distance of the two atoms
  real :: fij             ! force between two atoms
  real :: engcorrection   ! energy necessary shift the potential avoiding discontinuities
  engconf=0.0
  forces=0.0
  virial=0.0
 
  do iatom=1,natoms-1
     jatom=iatom+1
     dist=positions(iatom)-positions(jatom)
     call pbc(cell,dist,rij)
     engconf = engconf + deep/abs(rij) + 0.5*log(abs(rij))
     fij = (deep/(rij**2) - 0.5/abs(rij))*rij/abs(rij) 
     forces(iatom) = forces(iatom) + fij
     forces(jatom) = forces(jatom) - fij    
     virial = virial + rij*fij
  end do
     iatom=natoms
     jatom=1
     dist=positions(iatom)-positions(jatom)
     call pbc(cell,dist,rij)
     engconf = engconf + deep/abs(rij) + 0.5*log(abs(rij))
     fij = (deep/(rij**2) - 0.5/abs(rij))*rij/abs(rij)
     forces(iatom) = forces(iatom) + fij
     forces(jatom) = forces(jatom) - fij
     virial = virial + rij*fij 
end subroutine compute_forces

subroutine thermostat(natoms,masses,wmass,dt,friction,frictionL,temperature,velocities,ve,engint,idum)
! Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
! it is a linear combination of old velocities and new, randomly chosen, velocity,
! with proper coefficients
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: masses(natoms)
  real,    intent(in)    :: wmass
  real,    intent(in)    :: dt
  real,    intent(in)    :: friction
  real,    intent(in)    :: frictionL
  real,    intent(in)    :: temperature
  real,    intent(inout) :: velocities(natoms)
  real,    intent(inout) :: ve
  real,    intent(inout) :: engint ! contribution of phase-space compression
                                   ! its increment is equal to minus the kinetic-energy increment due to the
                                   ! thermostat
  integer, intent(inout) :: idum
  real :: c1 ! coefficient for the old velocity
  real :: c2 ! coefficient for the new velocity
  real :: b1
  real :: b2 
  real, external :: gasdev
  integer :: i,iatom
  c1=exp(-friction*dt)
  b1=exp(-frictionL*dt)
  b2=sqrt((1.0-b1**2)*temperature/wmass)
  engint=engint+0.5*wmass*ve**2
  ve=b1*ve+b2*gasdev(idum)
  engint=engint-0.5*wmass*ve**2
  do iatom=1,natoms
    c2=sqrt((1.0-c1**2)*temperature/masses(iatom))
      engint=engint+0.5*masses(iatom)*velocities(iatom)**2
      velocities(iatom)=c1*velocities(iatom)+c2*gasdev(idum)
      engint=engint-0.5*masses(iatom)*velocities(iatom)**2
  end do
end subroutine thermostat

subroutine compute_remove_vcm(natoms,masses,velocities)
! remove the velocity of center mass from the velocities
  implicit none
  integer, intent(in) :: natoms
  real,    intent(in) :: masses(natoms)
  real,    intent(inout) :: velocities(natoms)
  real :: vcm
  real :: mass
  integer :: iatom
  vcm=0.0
  mass=sum(masses(:))

  do iatom=1,natoms
    vcm=vcm+masses(iatom)*velocities(iatom)
  end do

  vcm=vcm/mass

  do iatom=1,natoms
    velocities(iatom)=velocities(iatom)-vcm
  end do
end subroutine compute_remove_vcm

subroutine compute_pressure(natoms,engkin,forces,positions,cell,virial)
! Instantaneous pressure from virial theorem
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: engkin
  real,    intent(in)    :: forces(natoms)
  real,    intent(in)    :: positions(natoms)
  real,    intent(in)    :: cell
  real,    intent(out)   :: virial
  integer :: iatom
  virial=0.0

  virial=virial+2.0*engkin

    do iatom=1,natoms
    virial=virial+forces(iatom)*positions(iatom)   ! potential energy contribution  
  end do

  virial=virial/cell
end subroutine compute_pressure

subroutine sinhx(x2,poly)
! series of tayle of sinh(x)/x
 implicit none
 real, intent(in)     :: x2
 real, intent(out)    :: poly
 real :: a2,a4,a6,a8,a10
 a2=1.0/(2*3);
 a4=a2/(4*5);
 a6=a4/(6*7);
 a8=a6/(8*9);
 a10=a8/(10*11);
 poly=1.0+x2*(a2+x2*(a4+x2*(a6+x2*(a8+x2*a10))))
end subroutine sinhx

subroutine temperature_system(temperature,Tstar,Tstop,istep,nstep)
 implicit none
 real,    intent(in)  :: Tstar
 real,    intent(in)  :: Tstop
 integer, intent(in)  :: istep
 integer, intent(in)  :: nstep
 real,    intent(out) :: temperature

 temperature=Tstar+(Tstop-Tstar)*istep/nstep

end subroutine temperature_system

subroutine pressure_system(pressure,Pstar,Pstop,istep,nstep)
 implicit none
 real,    intent(in)  :: Pstar
 real,    intent(in)  :: Pstop
 integer, intent(in)  :: istep
 integer, intent(in)  :: nstep
 real,    intent(out) :: pressure

 pressure=Pstar+(Pstop-Pstar)*istep/nstep

end subroutine pressure_system

subroutine mass_barostat(wmass,natoms,temperature,Tb)
 implicit none
 real,    intent(in)  :: temperature
 real,    intent(in)  :: Tb
 integer, intent(in)  :: natoms
 real,    intent(out) :: wmass

 wmass=natoms*temperature*Tb*Tb

end subroutine mass_barostat

end module routines




