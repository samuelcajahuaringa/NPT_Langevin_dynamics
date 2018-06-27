module routines
implicit none
contains
subroutine read_input(Tstart,Tstop,tstep,Tp, &
                      Pstart,Pstop,Tb, &  
                      forcecutoff,listcutoff,nstep,&
                      nequil,nconfig,nstat, &
                      inputfile,outputfile,trajfile,statfile, &
                      maxneighbours,ensemble,npt_style,idum)
  implicit none
  real,           intent(out) :: Tstart
  real,           intent(out) :: Tstop
  real,           intent(out) :: tstep
  real,           intent(out) :: Tp
  real,           intent(out) :: Pstart
  real,           intent(out) :: Pstop       
  real,           intent(out) :: Tb      
  real,           intent(out) :: forcecutoff
  real,           intent(out) :: listcutoff
  integer,        intent(out) :: nstep
  integer,        intent(out) :: nequil
  integer,        intent(out) :: nconfig
  integer,        intent(out) :: nstat
  integer,        intent(out) :: maxneighbours
  character(256), intent(out) :: ensemble 
  character(256), intent(out) :: npt_style
  character(256), intent(out) :: inputfile
  character(256), intent(out) :: outputfile
  character(256), intent(out) :: trajfile
  character(256), intent(out) :: statfile
  integer, intent(out) :: idum
  integer :: iostat
  character(256) :: line,keyword,keyword1
  integer :: i
  logical :: foundsharp
! default values
  Tstart=1.0
  Tstop=Tstart
  tstep=0.005
  Tp=1.0
  Pstart=0.0
  Pstop=Pstart
  Tb=5.0
  forcecutoff=2.5
  listcutoff=3.0
  nstep=1
  nconfig=10
  nstat=1
  nequil=1
  maxneighbours=1000
  idum=0
  statfile=""
  trajfile=""
  outputfile=""
  inputfile=""
  ensemble="nve"  
  npt_style="iso"
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
    case("ensemble")
      read(line,*) keyword1,ensemble
    case("npt_style")
      read(line,*) keyword1,npt_style
    case("temperature")
      read(line,*) keyword1,Tstart,Tstop
    case("tstep")
      read(line,*) keyword1,tstep
    case("Tp")
      read(line,*) keyword1,Tp
    case("pressure")
      read(line,*) keyword1,Pstart,Pstop  
    case("Tb")
      read(line,*) keyword1,Tb
!    case("wmass")
!      read(line,*) keyword1,wmass  
    case("forcecutoff")
      read(line,*) keyword1,forcecutoff
    case("listcutoff")
      read(line,*) keyword1,listcutoff
    case("nstep")
      read(line,*) keyword1,nstep
    case("nequil")
      read(line,*) keyword1,nequil
    case("nconfig")
      read(line,*) keyword1,nconfig,trajfile
    case("nstat")
      read(line,*) keyword1,nstat,statfile
    case("maxneighbours")
      read(line,*) keyword1,maxneighbours
    case("inputfile")
      read(line,*) keyword1,inputfile
    case("outputfile")
      read(line,*) keyword1,outputfile
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

subroutine triangular(a,n1,n2,natoms,positions,cell)
  implicit none
  real,    intent (in)  :: a
  integer, intent (in)  :: n1
  integer, intent (in)  :: n2
  integer, intent (in)  :: natoms
  real,    intent (out) :: positions(2,natoms)
  real,    intent (out) :: cell(2,2)
  integer :: i
  integer :: j
  integer :: k
  cell(1,1) = n1*a
  cell(2,2) = n2*sqrt(3.0)*a
  cell(1,2) = 0.0
  cell(2,1) = 0.0
  k=1

  do i=0,n1-1
     do j=0,n2-1
        positions(1,2*k-1) = i*a
        positions(2,2*k-1) = j*a*sqrt(3.0)
        positions(1,2*k) = (i+0.5)*a
        positions(2,2*k) = (j*sqrt(3.0)+sqrt(3.0)*0.5)*a
        k=k+1
     end do
  end do
end subroutine triangular

subroutine square(a,n1,n2,natoms,positions,cell)
  implicit none
  real,    intent (in)  :: a
  integer, intent (in)  :: n1
  integer, intent (in)  :: n2
  integer, intent (in)  :: natoms
  real,    intent (out) :: positions(2,natoms)
  real,    intent (out) :: cell(2,2)
  integer :: i
  integer :: j
  integer :: k

  cell(1,1) = n1*a
  cell(2,2) = n2*a
  cell(1,2) = 0.0
  cell(2,1) = 0.0

  k=1

  do i=0,n1-1
     do j=0,n2-1
        positions(1,k) = i*a
        positions(2,k) = j*a
        k=k+1
     end do
  end do
end subroutine square

subroutine read_positions(inputfile,natoms,positions,h)
! read positions and cell from a file called inputfile
! natoms (input variable) and number of atoms in the file should be consistent
  implicit none
  character(*), intent(in)  :: inputfile
  integer,      intent(in)  :: natoms
  real,         intent(out) :: positions(2,natoms)
  real,         intent(out) :: h(2,2)
  integer :: iatom
  character(100) :: atomname
  real :: pos(2)
  open(10,file=inputfile)
  
  read(10,*)
  read(10,*)h(1,1),h(1,2)
  read(10,*)h(2,1),h(2,2) 
  do iatom=1,natoms
    read(10,*)pos(1),pos(2)
    call put_in_box(h,pos,2)
    positions(1,iatom)=pos(1)
    positions(2,iatom)=pos(2)
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

subroutine write_positions(trajfile,natoms,positions,h,istep)
! write positions on xyzfile trajfile
! positions are appended at the end of the file
  implicit none
  character(*), intent(in) :: trajfile
  integer,      intent(in) :: natoms
  real,         intent(in) :: positions(2,natoms)
  real,         intent(in) :: h(2,2)
  integer,      intent(in) :: istep
  integer :: iatom
  real :: pos(2,natoms)
  real :: hh(2,2)
  real :: area
  real :: ri(2)
  logical, save :: first=.true.
  real :: a,b,ab
  real :: Ro(2,2)

! wrap atomics positions
  do iatom=1,natoms
    pos(1,iatom)=positions(1,iatom)
    pos(2,iatom)=positions(2,iatom)
    ri(1)=pos(1,iatom)
    ri(2)=pos(2,iatom)
    call put_in_box(h,ri,2)
    pos(1,iatom)=ri(1)
    pos(2,iatom)=ri(2)
  end do

  if(first) then
    open(10,file=trajfile)
    first=.false.
  else
    open(10,file=trajfile,position="append")
  end if
  ! configuration in extended XYZ format
1000 format(I4)
  write(10,1000)natoms
2000 format(A10,E12.5,2X,E12.5,2X,F7.3,2X,E12.5,2X,E12.5,2X,F7.3,2X,F7.3,2X,F7.3,2X,F7.3,2X,A26,I7)
  write(10,2000)'Lattice ="',h(1,1),h(2,1),0.0,h(1,2),h(2,2),0.0,1.0,0.0,0.0,'" Properties=pos:R:3 Time=',istep
3000 format(A2,2X,E12.5,2X,E12.5,2X,E12.5)
  do iatom=1,natoms
! usually, it is better not to apply pbc here, so that diffusion
! is more easily calculated from a trajectory file:
    !pos(1)=positions(1,iatom)
    !pos(2)=positions(2,iatom)
    !call put_in_box(cell,pos)
    write(10,3000)"Ar",pos(1,iatom),pos(2,iatom),0.5  
  end do
  close(10)
end subroutine write_positions

subroutine write_final_positions(outputfile,natoms,positions,cell)
! write positions on file outputfile
  character(*), intent(in) :: outputfile
  integer,      intent(in) :: natoms
  real,         intent(in) :: positions(2,natoms)
  real,         intent(in) :: cell(2,2)
  integer :: iatom
  real :: pos(2)
  open(10,file=outputfile)
  write(10,*) natoms
  write(10,*) cell(1,1),cell(1,2)
  write(10,*) cell(2,1),cell(2,2)
  do iatom=1,natoms
! usually, it is better not to apply pbc here, so that diffusion
! is more easily calculated from a trajectory file:
    pos(1)=positions(1,iatom)
    pos(2)=positions(2,iatom)
    call put_in_box(cell,pos,2)
    write(10,*) pos
  end do
  close(10)
end subroutine write_final_positions

subroutine write_statistics(statfile,istep,tstep,natoms,engkin,engconf,pressure,pinst,h,wmass,ve,engint)
! write statistics on file statfile
  character(*), intent(in) :: statfile
  integer,      intent(in) :: istep
  real,         intent(in) :: tstep
  integer,      intent(in) :: natoms
  real,         intent(in) :: engkin
  real,         intent(in) :: engconf
  real,         intent(in) :: pressure
  real,         intent(in) :: pinst(2,2) 
  real,         intent(in) :: h(2,2)
  real,         intent(in) :: ve(2,2)
  real,         intent(in) :: wmass
  real,         intent(in) :: engint
  real :: press                ! pressure
  real :: pxy                  ! non-diagonal pressure
  real :: area                 ! area of box
  real :: lx,ly                ! magnitud of box sides                
  real :: ab                   ! dot product of box vector 
  real :: angle                ! angle between lx and ly
  real :: energy               ! energy of system   
  real :: enthalpy             ! enthalpy energy 
  real :: Hpr                  ! Parrinello-Rahman energy
  real :: effective            ! Effective energy conservation
  real :: engkinbar 
  logical, save :: first=.true.
  integer, save :: last_time_reopened=0

  press=0.5*(pinst(1,1)+pinst(2,2))
  pxy=0.5*(pinst(1,2)+pinst(2,1))
  area=h(1,1)*h(2,2)-h(1,2)*h(2,1)
  lx=sqrt(h(1,1)*h(1,1)+h(2,1)*h(2,1))
  ly=sqrt(h(1,2)*h(1,2)+h(2,2)*h(2,2))
  ab = h(1,1)*h(2,1)+h(1,2)*h(2,2)
  angle=acos(ab/(lx*ly))*180.0*0.25/atan(1.0)
  energy=engkin+engconf
  enthalpy=energy+pressure*area
  engkinbar=0.5*wmass*(ve(1,1)*ve(1,1)+ve(2,2)*ve(2,2)+2.0*ve(1,2)*ve(1,2))
  Hpr=enthalpy+engkinbar
  effective=Hpr+engint

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
  
2000 format (I9,2X,F6.3,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5, &
     2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5) 
!3000 format (I9,2X,F6.3,2X,F12.8,2X,E12.5,2X,E12.5)
  ! step temp engconf energy area lx ly angle pressure pxx pyy pxy vexx veyy vexy
  write(666,2000) istep,2.0*engkin/(2.0*natoms),engconf,energy,area,lx,ly,angle,press,pxy,ve(1,1),ve(2,2),ve(1,2),enthalpy,effective
!  write(666,3000) istep,2.0*engkin/(2.0*natoms),pressure,area,press
end subroutine write_statistics

subroutine randomize_velocities(natoms,temperature,masses,velocities,idum)
! randomize the velocities according to the temperature
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: temperature
  real,    intent(in)    :: masses(natoms)
  real,    intent(out)   :: velocities(2,natoms)
  integer, intent(inout) :: idum
  real, external :: gasdev
  integer :: iatom,i
  do iatom=1,natoms
    do i=1,2
      velocities(i,iatom)=sqrt(temperature/masses(iatom))*gasdev(idum)
    end do
  end do
end subroutine randomize_velocities

subroutine compute_engkin(natoms,masses,velocities,engkin,engkin_tensor)
! calculate the kinetic energy from the velocities
  implicit none
  integer, intent(in)  :: natoms
  real,    intent(in)  :: masses(natoms)
  real,    intent(in)  :: velocities(2,natoms)
  real,    intent(out) :: engkin
  real,    intent(out) :: engkin_tensor(3)
  integer :: iatom
  engkin=0.0
  engkin_tensor=0.0
  do iatom=1,natoms
    engkin=engkin+0.5*masses(iatom)*sum(velocities(:,iatom)**2)
    engkin_tensor(1)=engkin_tensor(1)+0.5*masses(iatom)*velocities(1,iatom)**2
    engkin_tensor(2)=engkin_tensor(2)+0.5*masses(iatom)*velocities(2,iatom)**2
    engkin_tensor(3)=engkin_tensor(3)+0.5*masses(iatom)*velocities(1,iatom)*velocities(2,iatom)
  end do
end subroutine compute_engkin

subroutine compute_pressure(pinst,area,virial,engkin_tensor)
  implicit none
  real, intent(in)  :: area
  real, intent(in)  :: virial(2,2)
  real, intent(in)  :: engkin_tensor(3)
  real, intent(out) :: pinst(2,2)

  pinst(1,1)=(2.0*engkin_tensor(1)+virial(1,1))/area
  pinst(2,2)=(2.0*engkin_tensor(2)+virial(2,2))/area
  pinst(1,2)=(2.0*engkin_tensor(3)+virial(1,2))/area
  pinst(2,1)=(2.0*engkin_tensor(3)+virial(2,1))/area
end subroutine compute_pressure
  
subroutine pbc(h,ri,rj,rij,d)
! apply periodic boundary condition to a distance
  implicit none
  integer, intent(in) :: d
  real, intent(in)  :: h(d,d)
  real, intent(in)  :: ri(d)
  real, intent(in)  :: rj(d)
  real, intent(out) :: rij(d)
  real :: si(d)
  real :: sj(d)
  real :: sij(d)
  real :: h_inv(d,d)
  integer :: i
  real :: area
  
  area=h(1,1)*h(2,2)-h(2,1)*h(1,2)

  h_inv(1,1) = h(2,2)/area
  h_inv(2,2) = h(1,1)/area
  h_inv(1,2) = -h(1,2)/area
  h_inv(2,1) = -h(2,1)/area  
!  h_inv(1) = 1.0/h(1)
!  h_inv(2) = 1.0/h(2)
!  h_inv(3) = -h(3)*h_inv(1)*h_inv(2)

  si = ri

  call prod(h_inv,si,d)

  sj = rj

  call prod(h_inv,sj,d)


!  si(1) = ri(1)*h_inv(1) + ri(2)*h_inv(3)
!  si(2) = ri(2)*h_inv(2)

!  sj(1) = rj(1)*h_inv(1) + rj(2)*h_inv(3)
!  sj(2) = rj(2)*h_inv(2)

  do i=1,2
    sij(i) = si(i) - sj(i)
    sij(i) = sij(i)-nint(sij(i))
  end do

  call prod(h,sij,d)

  rij = sij

!  rij(1) = sij(1)*h(1) + sij(2)*h(3) 
!  rij(2) = sij(2)*h(2)

end subroutine pbc

subroutine put_in_box(h,ri,d)
  implicit none
  integer, intent(in) :: d
  real, intent(in)  :: h(d,d)
  real, intent(inout)  :: ri(d)
  real :: si(d)
  real :: h_inv(d,d)
  integer :: i
  real :: area
  
  area=h(1,1)*h(2,2)-h(2,1)*h(1,2)

  h_inv(1,1) = h(2,2)/area
  h_inv(2,2) = h(1,1)/area
  h_inv(1,2) = -h(1,2)/area
  h_inv(2,1) = -h(2,1)/area
 
  si = ri

  call prod(h_inv,si,d)

  si(1)=modulo(si(1),1.0)
  si(2)=modulo(si(2),1.0)

  call prod(h,si,d)

  ri = si 

end subroutine put_in_box 

subroutine check_list(natoms,positions,positions0,listcutoff,forcecutoff,recompute)
! check if the neighbour list have to be recomputed
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: positions(2,natoms)
  real,    intent(in)    :: positions0(2,natoms)
  real,    intent(in)    :: listcutoff
  real,    intent(in)    :: forcecutoff
  logical, intent(out)   :: recompute
  real :: displacement(2) ! displacement from positions0 to positions
  real :: delta2          ! square of the 'skin' thickness
  integer :: iatom
  recompute=.false.
  delta2=(0.5*(listcutoff-forcecutoff))**2
! if ANY atom moved more than half of the skin thickness, recompute is set to .true.
  do iatom=1,natoms
    displacement=positions(:,iatom)-positions0(:,iatom)
    if(sum(displacement**2)>delta2) recompute=.true.
  end do
end subroutine check_list

subroutine compute_list(natoms,listsize,positions,h,listcutoff,point,list)
  implicit none
  integer, intent(in)  :: natoms
  integer, intent(in)  :: listsize
  real,    intent(in)  :: positions(2,natoms)
  real,    intent(in)  :: h(2,2)
  real,    intent(in)  :: listcutoff
! see Allen-Tildesey for a definition of point and list
  integer, intent(out) :: point(natoms)
  integer, intent(out) :: list(listsize)
  integer :: iatom,jatom  ! indexes of the two involved atoms
  real :: rij(2)          ! distance of the two atoms
  real :: ri(2)
  real :: rj(2)
  real :: listcutoff2     ! squared list cutoff
  listcutoff2=listcutoff**2
  point(1)=1
  do iatom=1,natoms-1
    point(iatom+1)=point(iatom)
    do jatom=iatom+1,natoms
      ri=positions(:,iatom)
      rj=positions(:,jatom)    
      call pbc(h,ri,rj,rij,2)
! if the interparticle distance is larger than the cutoff, skip
      if(sum(rij**2)>listcutoff2) cycle
      if(point(iatom+1)>listsize) then
! too many neighbours
        write(0,*) "Verlet list size exceeded"
        write(0,*) "Increase maxneighbours"
        stop
      end if
      list(point(iatom+1))=jatom
      point(iatom+1)=point(iatom+1)+1
    end do
  end do
end subroutine compute_list

subroutine compute_forces(natoms,listsize,positions,h,c0,c2,forcecutoff,point,list,forces,engconf,virial)
  implicit none
  integer, intent(in)  :: natoms
  integer, intent(in)  :: listsize
  real,    intent(in)  :: positions(2,natoms)
  real,    intent(in)  :: h(2,2)
!  real,    intent(in)  :: eps
!  real,    intent(in)  :: sigma
  real,    intent(in)  :: c0
  real,    intent(in)  :: c2
  real,    intent(in)  :: forcecutoff
  integer, intent(in)  :: point(natoms)
  integer, intent(in)  :: list(listsize)
  real,    intent(out) :: forces(2,natoms)
  real,    intent(out) :: engconf
  real,    intent(out) :: virial(2,2)
  integer :: iatom,jatom  ! indexes of the two involved atoms
  integer :: jlist        ! counter for the neighbours of iatom
  real :: rij(2)          ! distance of the two atoms
  real :: ri(2)
  real :: rj(2)
  real :: rsq
  real :: forcecutoff2    ! squared force cutoff
  real :: Ulj             ! energy lennard-jones contribution
  real :: Ugauss          ! energy gauss contribution  
  real :: flj             ! force lennard-jones contribution
  real :: fgauss          ! force gauss contribution 
  real :: f(2)            ! force
  real :: w
  real :: lambda
  real :: ro
  w=5.0
  ro=1.5
  lambda=1.7
  forcecutoff2=forcecutoff**2
  engconf=0.0
  forces=0.0
  virial=0.0
  do iatom=1,natoms-1
    do jlist=point(iatom),point(iatom+1)-1
      jatom=list(jlist)
      ri=positions(:,iatom)
      rj=positions(:,jatom)
      call pbc(h,ri,rj,rij,2)
! distance_pbc is a vector
      rsq=sum(rij**2)
! if the interparticle distance is larger than the cutoff, skip
      if(rsq>forcecutoff2) cycle
      Ulj = 4.0*(1.0/rsq**6 - 1.0/rsq**3 + c2*rsq + c0)
      Ugauss = -lambda*exp(-(w*(sqrt(rsq)-ro))**2)
      engconf=engconf + Ulj + Ugauss
      flj = 4.0*(12.0/rsq**7 - 6.0/rsq**4 - 2.0*c2)
      fgauss = 2.0*w*w*(sqrt(rsq)-ro)*Ugauss
      f = rij*(flj + fgauss/sqrt(rsq))
! virial contribution to instantaneous pressure
      virial(1,1)=virial(1,1)+f(1)*rij(1)  
      virial(2,2)=virial(2,2)+f(2)*rij(2)
      virial(1,2)=virial(1,2)+f(1)*rij(2)
      virial(2,1)=virial(2,1)+f(2)*rij(1)
! same force on the two atoms, with opposite sign:
      forces(:,iatom)=forces(:,iatom)+f
      forces(:,jatom)=forces(:,jatom)-f
    end do
  end do
end subroutine compute_forces

subroutine thermostat_particles(natoms,masses,dt,friction,temperature,velocities,engint,idum)
! Langevin thermostat, implemented as decribed in Bussi and Parrinello, Phys. Rev. E (2007)
! it is a linear combination of old velocities and new, randomly chosen, velocity,
! with proper coefficients
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: masses(natoms)
  real,    intent(in)    :: dt
  real,    intent(in)    :: friction
  real,    intent(in)    :: temperature
  real,    intent(inout) :: velocities(2,natoms)
  real,    intent(inout) :: engint ! contribution of phase-space compression
                                   ! its increment is equal to minus the kinetic-energy increment due to the
                                   ! thermostat
  integer, intent(inout) :: idum
  real :: c1 ! coefficient for the old velocity
  real :: c2 ! coefficient for the new velocity
  real, external :: gasdev
  integer :: i,iatom

  c1=exp(-friction*dt)
  do iatom=1,natoms
    c2=sqrt((1.0-c1**2)*temperature/masses(iatom))
    do i=1,2
      engint=engint+0.5*masses(iatom)*velocities(i,iatom)**2
      velocities(i,iatom)=c1*velocities(i,iatom)+c2*gasdev(idum)
      engint=engint-0.5*masses(iatom)*velocities(i,iatom)**2
    end do
  end do

end subroutine thermostat_particles

subroutine thermostat_barostat_iso(wmass,dt,frictionL,temperature,ve,engint,idum)
  implicit none
  real,    intent(in)    :: wmass
  real,    intent(in)    :: dt
  real,    intent(in)    :: frictionL
  real,    intent(in)    :: temperature
  real,    intent(inout) :: ve
  real,    intent(inout) :: engint ! contribution of phase-space compression
                                   ! its increment is equal to minus the kinetic-energy increment due to the
                                   ! thermostat
  integer, intent(inout) :: idum
  real :: b1
  real :: b2
  real, external :: gasdev

  b1=exp(-frictionL*dt)
  b2=sqrt((1.0-b1**2)*temperature/wmass)
  engint=engint+0.5*wmass*ve*ve
  ve=b1*ve+b2*gasdev(idum)
  engint=engint-0.5*wmass*ve*ve

end subroutine thermostat_barostat_iso

subroutine thermostat_barostat_flex(wmass,dt,frictionL,temperature,ve,engint,idum)
  implicit none
  real,    intent(in)    :: wmass
  real,    intent(in)    :: dt
  real,    intent(in)    :: frictionL
  real,    intent(in)    :: temperature
  real,    intent(inout) :: ve(2,2)
  real,    intent(inout) :: engint ! contribution of phase-space compression
                                   ! its increment is equal to minus the kinetic-energy increment due to the
                                   ! thermostat
  integer, intent(inout) :: idum
  real :: b1
  real :: b2
  real, external :: gasdev

  b1=exp(-frictionL*dt)
  b2=sqrt((1.0-b1**2)*temperature/wmass)
  engint=engint+0.5*wmass*(ve(1,1)*ve(1,1)+ve(2,2)*ve(2,2)+2.0*ve(1,2)*ve(1,2))
  ve(1,1)=b1*ve(1,1)+b2*gasdev(idum)
  ve(2,2)=b1*ve(2,2)+b2*gasdev(idum)
  ve(1,2)=b1*ve(1,2)+b2*gasdev(idum)
  ve(2,1)=ve(1,2)
  engint=engint-0.5*wmass*(ve(1,1)*ve(1,1)+ve(2,2)*ve(2,2)+2.0*ve(1,2)*ve(1,2))

end subroutine thermostat_barostat_flex

subroutine velocity_barostat_iso(ve,area,pinst,pressure,natoms,engkin,tstep,wmass)
  implicit none
  real,    intent(in)    :: area
  real,    intent(in)    :: pinst(2,2)
  real,    intent(in)    :: pressure
  integer, intent(in)    :: natoms
  real,    intent(in)    :: engkin
  real,    intent(in)    :: tstep
  real,    intent(in)    :: wmass
  real,    intent(inout) :: ve
  real :: trz
  trz=0.5*(pinst(1,1)+pinst(2,2))

  ve=ve+(2.0*area*(trz-pressure)+2.0/(2.0*natoms)*2.0*engkin)*0.5*tstep/wmass
end subroutine velocity_barostat_iso

subroutine velocity_barostat_flex(ve,area,pinst,pressure,natoms,engkin,tstep,wmass)
  implicit none
  real,    intent(in)    :: area
  real,    intent(in)    :: pinst(2,2)
  real,    intent(in)    :: pressure
  integer, intent(in)    :: natoms
  real,    intent(in)    :: engkin
  real,    intent(in)    :: tstep
  real,    intent(in)    :: wmass
  real,    intent(inout) :: ve(2,2)
  real :: pxy
  
  pxy = 0.5*(pinst(1,2)+pinst(2,1))

  ve(1,1)=ve(1,1)+(area*(pinst(1,1)-pressure)+1.0/(2.0*natoms)*2.0*engkin)*0.5*tstep/wmass
  ve(2,2)=ve(2,2)+(area*(pinst(2,2)-pressure)+1.0/(2.0*natoms)*2.0*engkin)*0.5*tstep/wmass
  ve(1,2)=ve(1,2)+area*pxy*0.5*tstep/wmass
  ve(2,1)=ve(1,2)
end subroutine velocity_barostat_flex

subroutine velocities_particles_iso(natoms,velocities,masses,forces,ve,tstep)
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: masses(natoms)
  real,    intent(inout) :: velocities(2,natoms)
  real,    intent(in)    :: forces(2,natoms)
  real,    intent(in)    :: ve
  real,    intent(in)    :: tstep
  integer :: iatom
  real :: lambda
  real :: alpha
  real :: vel_scale
  real :: vel_poly

  alpha=2.0/(2.0*natoms)
  lambda=-alpha*ve*0.25*tstep 
  vel_scale=exp(lambda)
  call sinhx(lambda**2,vel_poly)

do iatom=1,natoms
   velocities(:,iatom)=vel_scale*(vel_scale*velocities(:,iatom)+vel_poly*forces(:,iatom)*0.5*tstep/masses(iatom))
end do
end subroutine velocities_particles_iso

subroutine velocities_particles_flex(natoms,velocities,masses,forces,O,OT,ve1,ve2,tstep)
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: masses(natoms)
  real,    intent(inout) :: velocities(2,natoms)
  real,    intent(in)    :: forces(2,natoms)
  real,    intent(in)    :: O(2,2)
  real,    intent(in)    :: OT(2,2)
  real,    intent(in)    :: ve1
  real,    intent(in)    :: ve2
  real,    intent(in)    :: tstep
  real :: v(2)
  real :: f(2)
  integer :: iatom
  real :: lambda1
  real :: lambda2
  real :: b
  real :: vel_scale1
  real :: vel_scale2
  real :: vel_poly1
  real :: vel_poly2
  real :: D1(2,2)
  real :: D2(2,2)
  real :: vector1(2)
  real :: vector2(2)
  
  b=1.0/(2.0*natoms)
  lambda1=(ve1+b*(ve1+ve2))*0.25*tstep
  lambda2=(ve2+b*(ve1+ve2))*0.25*tstep
  vel_scale1=exp(-lambda1)
  vel_scale2=exp(-lambda2)  
  call sinhx(lambda1**2,vel_poly1)
  call sinhx(lambda2**2,vel_poly2)

  D1(1,1)=vel_scale1**2
  D1(2,2)=vel_scale2**2
  D1(1,2)=0.0
  D1(2,1)=0.0
    
  D2(1,1)=vel_scale1*vel_poly1
  D2(2,2)=vel_scale2*vel_poly2
  D2(1,2)=0.0
  D2(2,1)=0.0

do iatom=1,natoms
  v(1)=velocities(1,iatom)
  v(2)=velocities(2,iatom)

  call prod(OT,v,2)
  call prod(D1,v,2)
  call prod(O,v,2)

  f(1)=forces(1,iatom)
  f(2)=forces(2,iatom)

  call prod(OT,f,2)
  call prod(D2,f,2)
  call prod(O,f,2)

  velocities(1,iatom)=v(1)+0.5*f(1)*tstep/masses(iatom)
  velocities(2,iatom)=v(2)+0.5*f(2)*tstep/masses(iatom) 
end do    

end subroutine velocities_particles_flex

subroutine positions_particles_iso(natoms,positions,velocities,ve,tstep)
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(inout) :: positions(2,natoms)
  real,    intent(in)    :: velocities(2,natoms)
  real,    intent(in)    :: ve
  real,    intent(in)    :: tstep
  integer :: iatom
  real :: lambda
  real :: pos_scale
  real :: pos_poly
  lambda=ve*0.5*tstep
  pos_scale=exp(lambda)
  call sinhx(lambda**2,pos_poly)

do iatom=1,natoms
    positions(:,iatom)=pos_scale*(pos_scale*positions(:,iatom)+pos_poly*velocities(:,iatom)*tstep)
end do
end subroutine positions_particles_iso

subroutine positions_particles_flex(natoms,positions,velocities,O,OT,ve1,ve2,tstep)
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(inout) :: positions(2,natoms)
  real,    intent(in)    :: velocities(2,natoms)
  real,    intent(in)    :: O(2,2)
  real,    intent(in)    :: OT(2,2)
  real,    intent(in)    :: ve1
  real,    intent(in)    :: ve2
  real,    intent(in)    :: tstep
  real :: v(2)
  real :: r(2)
  integer :: iatom
  real :: lambda1
  real :: lambda2
  real :: vel_scale1
  real :: vel_scale2
  real :: vel_poly1
  real :: vel_poly2
  real :: D1(2,2)
  real :: D2(2,2)
  real :: vector1(2)
  real :: vector2(2)

  lambda1=ve1*0.5*tstep
  lambda2=ve2*0.5*tstep
  vel_scale1=exp(lambda1)
  vel_scale2=exp(lambda2)
  call sinhx(lambda1**2,vel_poly1)
  call sinhx(lambda2**2,vel_poly2)

  D1(1,1)=vel_scale1**2
  D1(2,2)=vel_scale2**2
  D1(1,2)=0.0
  D1(2,1)=0.0

  D2(1,1)=vel_scale1*vel_poly1
  D2(2,2)=vel_scale2*vel_poly2
  D2(1,2)=0.0
  D2(2,1)=0.0


do iatom=1,natoms
  r(1)=positions(1,iatom)
  r(2)=positions(2,iatom)

  call prod(OT,r,2)
  call prod(D1,r,2)
  call prod(O,r,2)

  v(1)=velocities(1,iatom)
  v(2)=velocities(2,iatom)

  call prod(OT,v,2)
  call prod(D2,v,2)
  call prod(O,v,2)

  positions(1,iatom)=r(1)+v(1)*tstep
  positions(2,iatom)=r(2)+v(2)*tstep
end do

end subroutine positions_particles_flex

subroutine box_propagated_iso(h,ve,tstep)
  implicit none
  real, intent(inout) :: h(2,2)
  real, intent(in)    :: ve
  real, intent(in)    :: tstep

  h=exp(ve*tstep)*h
 
end subroutine box_propagated_iso

subroutine box_propagated_flex(h,Obox,OTbox,h1,h2,h3,h4,tstep)
  implicit none
  real, intent(inout) :: h(2,2)
  real, intent(in)    :: Obox(4,4)
  real, intent(in)    :: OTbox(4,4)
  real, intent(in)    :: h1,h2,h3,h4
  real, intent(in)    :: tstep
  real :: hg(4)
  real :: Dg(4,4)

  hg(1) = h(1,1)
  hg(2) = h(1,2)
  hg(3) = h(2,1)
  hg(4) = h(2,2)

  Dg = 0.0
  Dg(1,1) = exp(h1*tstep)
  Dg(2,2) = exp(h2*tstep)
  Dg(3,3) = exp(h3*tstep)  
  Dg(4,4) = exp(h4*tstep)
 
  call prod(OTbox,hg,4)
  call prod(Dg,hg,4) 
  call prod(Obox,hg,4)

  h(1,1) = hg(1)
  h(1,2) = hg(2)
  h(2,1) = hg(3)
  h(2,2) = hg(4)
end subroutine box_propagated_flex

subroutine orthogonal_array(ve,O,OT,ve1,ve2)
  implicit none
  real, intent(in)  :: ve(2,2)
  real, intent(out) :: O(2,2)
  real, intent(out) :: OT(2,2)
  real, intent(out) :: ve1
  real, intent(out) :: ve2
  real, parameter :: abserr=1.0e-10
  real :: vg(2,2)
  real :: u1(2)
  real :: u2(2)
  real :: v1(2)
  real :: v2(2)

  vg(1,1) = ve(1,1)
  vg(2,2) = ve(2,2)
  vg(1,2) = ve(1,2)
  vg(2,1) = ve(2,1)

  call jacobi(vg,O,abserr,2)
 
  ve1=vg(1,1)
  ve2=vg(2,2) 

  v1(1)=O(1,1)
  v1(2)=O(2,1)

  u1(1)=v1(1) !/(v1(1)**2 + v1(2)**2)
  u1(2)=v1(2) !/(v1(1)**2 + v1(2)**2)
  
  v2(1)=O(1,2)
  v2(2)=O(2,2)

  u2(1)=v2(1)-(v2(1)*u1(1)+v2(2)*u1(2))*u1(1)
  u2(2)=v2(2)-(v2(1)*u1(1)+v2(2)*u1(2))*u1(2)

  O(1,1)=u1(1)
  O(2,1)=u1(2)
  O(1,2)=u2(1)
  O(2,2)=u2(2)

  OT(1,1)=O(1,1)
  OT(2,2)=O(2,2)
  OT(1,2)=O(2,1)
  OT(2,1)=O(1,2)
end subroutine orthogonal_array

subroutine orthogonal_box_array(ve,Obox,OTbox,h1,h2,h3,h4)
  implicit none
  real, intent(in)  :: ve(2,2)
  real, intent(out) :: Obox(4,4)
  real, intent(out) :: OTbox(4,4)
  real, intent(out) :: h1
  real, intent(out) :: h2
  real, intent(out) :: h3
  real, intent(out) :: h4
  real, parameter :: abserr=1.0e-10
  real :: vg(4,4)
  real :: u1(4)
  real :: u2(4)
  real :: u3(4)
  real :: u4(4)
  real :: v1(4)
  real :: v2(4)
  real :: v3(4)
  real :: v4(4)
  integer :: i,j,k

  vg(1,1) = ve(1,1)
  vg(1,2) = 0.0
  vg(1,3) = ve(1,2)
  vg(1,4) = 0.0

  vg(2,1) = 0.0
  vg(2,2) = ve(1,1)
  vg(2,3) = 0.0
  vg(2,4) = ve(1,2)

  vg(3,1) = ve(2,1)
  vg(3,2) = 0.0
  vg(3,3) = ve(2,2)
  vg(3,4) = 0.0

  vg(4,1) = 0.0
  vg(4,2) = ve(2,1)
  vg(4,3) = 0.0
  vg(4,4) = ve(2,2)

  call jacobi(vg,Obox,abserr,4)

  h1 = vg(1,1)
  h2 = vg(2,2)
  h3 = vg(3,3)
  h4 = vg(4,4)

  do i=1,4
     v1(i)=Obox(i,1)
     v2(i)=Obox(i,2)
     v3(i)=Obox(i,3)
     v4(i)=Obox(i,4)
  end do

  do i=1,4
     u1(i) = v1(i) 
  end do

  do i=1,4
     u2(i) = v2(i)-(v2(1)*u1(1)+v2(2)*u1(2)+v2(3)*u1(3)+v2(4)*u1(4))*u1(i)
  end do

  do i=1,4
     u2(i) = u2(i)/sqrt(u2(1)**2 + u2(2)**2 + u2(3)**2 + u2(4)**2)
  end do

  do i=1,4
     u3(i) = v3(i)-(v3(1)*u1(1)+v3(2)*u1(2)+v3(3)*u1(3)+v3(4)*u1(4))*u1(i) -(v3(1)*u2(1)+v3(2)*u2(2)+v3(3)*u2(3)+v3(4)*u2(4))*u2(i) 
  end do

  do i=1,4
     u3(i) = u3(i)/sqrt(u3(1)**2 + u3(2)**2 + u3(3)**2 + u3(4)**2)
  end do

  do i=1,4
     u4(i) = v4(i)-(v4(1)*u1(1)+v4(2)*u1(2)+v4(3)*u1(3)+v4(4)*u1(4))*u1(i) -(v4(1)*u2(1)+v4(2)*u2(2)+v4(3)*u2(3)+v4(4)*u2(4))*u2(i) 
     u4(i) = u4(i)-(v4(1)*u3(1)+v4(2)*u3(2)+v4(3)*u3(3)+v4(4)*u3(4))*u3(i)
  end do

  do i=1,4
     u4(i) = u4(i)/sqrt(u4(1)**2 + u4(2)**2 + u4(3)**2 + u4(4)**2)
  end do

  do i=1,4
     Obox(i,1)=u1(i)
     Obox(i,2)=u2(i)
     Obox(i,3)=u3(i)
     Obox(i,4)=u4(i)
  end do

  do i=1,4  
     do j=1,4
        OTbox(i,j)=Obox(j,i)
     end do
  end do
end subroutine orthogonal_box_array

subroutine compute_vcm(natoms,masses,velocities,vcm)
  implicit none
  integer, intent(in) :: natoms
  real,    intent(in) :: masses(natoms)
  real,    intent(in) :: velocities(2,natoms)
  real,    intent(out) :: vcm(2)
  real :: mass
  integer :: iatom
  vcm=0.0
  mass=sum(masses(:))

  do iatom=1,natoms
    vcm(:)=vcm(:)+masses(iatom)*velocities(:,iatom)
  end do

  vcm(:)=vcm(:)/mass

end subroutine compute_vcm

subroutine compute_remove_vcm(natoms,masses,velocities)
! remove the velocity of center mass from the velocities
  implicit none
  integer, intent(in) :: natoms
  real,    intent(in) :: masses(natoms)
  real,    intent(inout) :: velocities(2,natoms)
  real :: vcm(2)
  real :: mass
  integer :: iatom
  vcm=0.0
  mass=sum(masses(:))
  
  do iatom=1,natoms
    vcm(:)=vcm(:)+masses(iatom)*velocities(:,iatom)
  end do
  
  vcm(:)=vcm(:)/mass
  
  do iatom=1,natoms
    velocities(:,iatom)=velocities(:,iatom)-vcm(:)
  end do 
end subroutine compute_remove_vcm

subroutine compute_verlet_velocities(natoms,velocities,masses,forces,tstep)
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: masses(natoms)
  real,    intent(inout) :: velocities(2,natoms)
  real,    intent(in)    :: forces(2,natoms)
  real,    intent(in)    :: tstep
  integer :: iatom
 
  do iatom=1,natoms
     velocities(:,iatom)=velocities(:,iatom)+forces(:,iatom)*0.5*tstep/masses(iatom)
  end do

end subroutine compute_verlet_velocities

subroutine compute_verlet_positions(natoms,positions,velocities,tstep)
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(in)    :: velocities(2,natoms)
  real,    intent(inout) :: positions(2,natoms)
  real,    intent(in)    :: tstep
  integer :: iatom

  do iatom=1,natoms
    positions(:,iatom)=positions(:,iatom)+velocities(:,iatom)*tstep
  end do

end subroutine compute_verlet_positions

subroutine put_in_rcm_frame(natoms,positions,cell,rcm)
! put the atoms in box across pass the boundary limits
  implicit none
  integer, intent(in) :: natoms
  real,    intent(in) :: cell(3)
  real,    intent(in) :: rcm(3)
  real,    intent(inout) ::  positions(3,natoms)
  integer :: i,j

  do i=1,3
     do j=1,natoms
        do
          if (positions(i,j) >= cell(i)-rcm(i)) then
            positions(i,j)=positions(i,j)-cell(i)
          else if (positions(i,j) < -rcm(i)) then
            positions(i,j)=positions(i,j)+cell(i)
          else
            exit
          end if
        end do
     end do
  end do
end subroutine put_in_rcm_frame

subroutine compute_center_mass(natoms,masses,positions,rcm)
! calculated the center mass of system
  implicit none
  integer, intent(in)  :: natoms 
  real,    intent(in)  :: positions(2,natoms)
  real,    intent(in)  :: masses(natoms)
  real,    intent(out) :: rcm(2)
  integer :: iatom
  real    :: mass

  mass=sum(masses(:))

  do iatom=1,natoms
    rcm(:)=rcm(:)+masses(iatom)*positions(:,iatom)
  end do

  rcm(:)=rcm(:)/mass
end subroutine compute_center_mass

subroutine compute_positions_rcm(natoms,masses,positions,rcm)
! calculated the position relative at center mass frame
  implicit none
  integer, intent(in)    :: natoms
  real,    intent(inout) :: positions(2,natoms)
  real,    intent(in)    :: masses(natoms)
  real,    intent(out)   :: rcm(2)
  integer :: iatom
  real    :: mass

  mass=sum(masses(:))

  do iatom=1,natoms
    rcm(:)=rcm(:)+masses(iatom)*positions(:,iatom)
  end do

  rcm(:)=rcm(:)/mass

  do iatom=1,natoms
    positions(:,iatom)=positions(:,iatom)-rcm(:)
  end do
end subroutine compute_positions_rcm

subroutine sinhx(x2,poly)
! series of tayle of sinh(x)/x
 implicit none
 real,    intent(in)     :: x2
 real,    intent(out)    :: poly
 real                    :: a2,a4,a6,a8
 a2=1.0/(2*3);
 a4=a2/(4*5);
 a6=a4/(6*7);
 a8=a6/(8*9); 

 poly=1.0+x2*(a2+x2*(a4+x2*(a6+x2*a8)))
 
end subroutine sinhx

subroutine coeff(sigma,rc,c0,c2)

 implicit none
 real, intent(in)  :: sigma
 real, intent(in)  :: rc 
 real, intent(out) :: c0
 real, intent(out) :: c2
 real :: a1
 a1 = sigma/rc
 c0 = 4.0*a1**6 - 7.0*a1**12
 c2 = -3.0*a1**8 + 6.0*a1**14
end subroutine coeff

subroutine jacobi(a,x,abserr,n)
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
implicit none
integer :: i, j, k
integer, intent(in) :: n
real, intent(inout) :: a(n,n),x(n,n)
real, intent(in) :: abserr
real :: b2, bar
real :: beta, coeff, c, s, cs, sc

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
x = 0.0
do i=1,n
  x(i,i) = 1.0
end do

! find the sum of all off-diagonal elements (squared)
b2 = 0.0
do i=1,n
  do j=1,n
    if (i.ne.j) b2 = b2 + a(i,j)**2
  end do
end do

if (b2 <= abserr) return

! average for off-diagonal elements /2
bar = 0.5*b2/float(n*n)

do while (b2.gt.abserr)
  do i=1,n-1
    do j=i+1,n
      if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
      b2 = b2 - 2.0*a(j,i)**2
      bar = 0.5*b2/float(n*n)
! calculate coefficient c and s for Givens matrix
      beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
      coeff = 0.5*beta/sqrt(1.0+beta**2)
      s = sqrt(max(0.5+coeff,0.0))
      c = sqrt(max(0.5-coeff,0.0))
! recalculate rows i and j
      do k=1,n
        cs =  c*a(i,k)+s*a(j,k)
        sc = -s*a(i,k)+c*a(j,k)
        a(i,k) = cs
        a(j,k) = sc
      end do
! new matrix a_{k+1} from a_{k}, and eigenvectors 
      do k=1,n
        cs =  c*a(k,i)+s*a(k,j)
        sc = -s*a(k,i)+c*a(k,j)
        a(k,i) = cs
        a(k,j) = sc
        cs =  c*x(k,i)+s*x(k,j)
        sc = -s*x(k,i)+c*x(k,j)
        x(k,i) = cs
        x(k,j) = sc
      end do
    end do
  end do
end do
return
end subroutine jacobi

subroutine prod(a,v,d)
 implicit none
 integer, intent(in)     :: d
 real,    intent(in)     :: a(d,d)
 real,    intent(inout)  :: v(d)
 real :: b(d)
 integer :: i,j
 
 do i=1,d 
    b(i) = 0.0
    do j=1,d
       b(i) = b(i) + a(i,j)*v(j)
    end do
 end do

 do i=1,d
    v(i)=b(i)
 end do
end subroutine prod

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

subroutine mass_barostat(wmass,natoms,temperature,Tb,d,style)
 implicit none
 real,    intent(in)  :: temperature
 real,    intent(in)  :: Tb
 integer, intent(in)  :: natoms
 integer, intent(in)  :: d
 integer, intent(in)  :: style
 real,    intent(out) :: wmass

 if (style == 1) then  
    wmass=(d*natoms+d)*temperature*Tb*Tb/d
 else 
    wmass=(d*natoms+d)*temperature*Tb*Tb
 end if
end subroutine mass_barostat
 
end module routines
