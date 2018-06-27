! particle in one-dimensional periodic potential in the NPT ensemble by the Langevin dynamics
! author: Oscar Samuel Cajahuaringa Macollunco and Alex Antonelli
! copyright 2018 Samuel Cajahuaringa <samuelif@ifi.unicamp.br>
! article: Stochastic sampling of the isothermal-isobaric ensemble: phase diagram of crystalline solids 
! this is free software: you can redistribute it and/or modify it under the terms of 
! the GNU General Public License as published by the Free Software Foundation.
! The code is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even 
! the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
! See the GNU General Public License for more details.

program nptlangevin
implicit none

real :: xi ! positions
real :: pi ! momentum
real :: fi ! force
real :: mi ! mass
real :: Li ! length of box
real :: ve ! velocity of box
real :: W  ! mass of barostat
real :: Ge
real :: kT ! temperature energy
real :: Pext  ! pressure
real :: Tp
real :: Tb
real :: gp
real :: gb
real :: a1
real :: a2
real :: b1
real :: b2
! barostat 
real :: vscale   ! barostat velocities scale
real :: xscale   ! barostat postions scale
real :: dt ! timestep
integer :: nrun  ! number of steps of thermalization
integer :: nequil ! number of steps of equilibration
integer :: nstat ! stride for output statistics
real :: Hi ! instantaneous energy of system
real :: Ho ! initial energy of system
real :: dH ! error in the energy of system
real :: tmp
real    :: engkin         ! kinetic energy
real    :: engpot         ! configurational energy
real    :: engint
real    :: Pinst
integer :: istep          ! step counter
integer :: seed
real, external :: gasdev

! define variables
call read_input(nrun,nequil,nstat,dt,xi,pi,mi,Li,ve,kT,Tp,Tb,Pext,seed)

gp=0.5/Tp
gb=0.1/Tb

W=2.0*kT*Tb*Tb

a1=exp(-0.5*gp*dt)
a2=sqrt((1.0-a1**2)*kT*mi)
b1=exp(-0.5*gb*dt)
b2=sqrt((1.0-b1**2)*kT/W)
fi=force(xi,Li)
engpot=potential(xi,Li)
engkin=0.5*pi*pi/mi
engint=0.0
! write the parameters in output so they can be checked
write(*,*) "#position               : ",xi
write(*,*) "#momentum               : ",pi
write(*,*) "#mass                   : ",mi
write(*,*) "#box length             : ",Li
write(*,*) "#box velocity           : ",ve
write(*,*) "#barostat mass          : ",W
write(*,*) "#Temperature            : ",kT
write(*,*) "#Pressure               : ",Pext
write(*,*) "#engkin                 : ",engkin
write(*,*) "#engpot                 : ",engpot
write(*,*) "#time step              : ",dt
write(*,*) "#Number of steps Run    : ",nrun
write(*,*) "#Number of steps Equil  : ",nequil
write(*,*) "#Stride for statistics  : ",nstat
write(*,*) "#-------------------------"
write(*,*) "#step qi pi Li ve Epot Etotal Enthalpy Effective_Enthalpy"

do istep=1,nrun

  ve=b1*ve+b2*gasdev(seed)

  pi=a1*pi+a2*gasdev(seed)

  engkin = 0.5*pi*pi/mi

  Pinst = (2.0*engkin-2.0*engpot)/Li

  Ge = Li*(Pinst-Pext) + 2.0*engkin

  ve = ve + 0.5*dt*Ge/W

  vscale = exp(-0.5*ve*dt)
  pi = pi*vscale**2 + 0.5*dt*fi*vscale*sinhx(-0.5*ve*dt)

  xscale = exp(0.5*ve*dt)
  xi = xi*xscale**2 + dt*pi/mi*xscale*sinhx(0.5*ve*dt)

  Li = Li*exp(ve*dt)

  if (xi > Li) then
     xi = xi - Li
  elseif (xi < 0.0) then
     xi = xi + Li
  end if

  fi = force(xi,Li)

  pi = pi*vscale**2 + 0.5*dt*fi*vscale*sinhx(-0.5*ve*dt)

  engkin = 0.5*pi*pi/mi

  engpot = potential(xi,Li)

  Pinst = (2.0*engkin-2.0*engpot)/Li

  Ge = Li*(Pinst-Pext) + 2.0*engkin

  ve = ve + 0.5*dt*Ge/W

  ve=b1*ve+b2*gasdev(seed)

  pi=a1*pi+a2*gasdev(seed)

end do

1000 format (F8.1,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)
2000 format (F12.4,2X,E12.5,2X,E12.5)

do istep=1,nequil

  engint=engint+0.5*pi*pi/mi+0.5*W*ve*ve

  ve=b1*ve+b2*gasdev(seed)

  pi=a1*pi+a2*gasdev(seed)

  engint=engint-0.5*pi*pi/mi-0.5*W*ve*ve

  engkin = 0.5*pi*pi/mi

  Pinst = (2.0*engkin-2.0*engpot)/Li
 
  Ge = Li*(Pinst-Pext) + 2.0*engkin
 
  ve = ve + 0.5*dt*Ge/W

  vscale = exp(-0.5*ve*dt)
  pi = pi*vscale**2 + 0.5*dt*fi*vscale*sinhx(-0.5*ve*dt)
  
  xscale = exp(0.5*ve*dt)
  xi = xi*xscale**2 + dt*pi/mi*xscale*sinhx(0.5*ve*dt)

  Li = Li*exp(ve*dt)
  
  if (xi > Li) then
     xi = xi - Li
  elseif (xi < 0.0) then
     xi = xi + Li
  end if
  
  fi = force(xi,Li)

  pi = pi*vscale**2 + 0.5*dt*fi*vscale*sinhx(-0.5*ve*dt)
 
  engkin = 0.5*pi*pi/mi

  engpot = potential(xi,Li) 

  Pinst = (2.0*engkin-2.0*engpot)/Li

  Ge = Li*(Pinst-Pext) + 2.0*engkin

  ve = ve + 0.5*dt*Ge/W

  engint=engint+0.5*pi*pi/mi+0.5*W*ve*ve

  ve=b1*ve+b2*gasdev(seed)

  pi=a1*pi+a2*gasdev(seed)

  engint=engint-0.5*pi*pi/mi-0.5*W*ve*ve

  engkin = 0.5*pi*pi/mi

! eventually, write positions and statistics
  if(modulo(istep,nstat)==0)  then 
    Hi = engkin + engpot + Pext*Li 
    tmp = Hi + 0.5*W*ve*ve + engint  
    write(*,1000)istep*dt,xi,pi,Li,ve,engpot,engpot+engkin,Hi,tmp
  end if 
 
end do

write(*,*)"#Execution completed"

contains

subroutine read_input(nrun,nequil,nstat,dt,xi,pi,mi,Li,ve,kT,Tp,Tb,Pext,seed)
  implicit none
  integer,        intent(out) :: nrun
  integer,        intent(out) :: nequil
  integer,        intent(out) :: nstat
  real,           intent(out) :: dt
  real,           intent(out) :: xi
  real,           intent(out) :: pi
  real,           intent(out) :: mi 
  real,           intent(out) :: Li
  real,           intent(out) :: ve
!  real,           intent(out) :: W
  real,           intent(out) :: kT
  real,           intent(out) :: Tp
  real,           intent(out) :: Tb
  real,           intent(out) :: Pext
  integer,        intent(out) :: seed!
!  character(256), intent(out) :: inputfile

  integer :: iostat
  character(256) :: line,keyword,keyword1
  integer :: i
  logical :: foundsharp
! default values
  nrun=1
  nequil=1
  nstat=1
  dt=0.01
  xi=0.0
  pi=0.0
  mi=1.0
  Li=8.0
  ve=0.0
  !W=18.0
  kT=1.0
  Tp=1.0
  Tb=3.0
  Pext=1.0
  seed=0
!  inputfile=""

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
    case("nrun")
      read(line,*) keyword1,nrun
    case("nequil")
      read(line,*) keyword1,nequil
    case("nstat")
      read(line,*) keyword1,nstat
    case("dt")
      read(line,*) keyword1,dt
    case("xi")
      read(line,*) keyword1,xi
    case("pi")
      read(line,*) keyword1,pi
    case("mi")
      read(line,*) keyword1,mi
    case("Li")
      read(line,*) keyword1,Li
    case("ve")
      read(line,*) keyword1,ve
    case("kT")
      read(line,*) keyword1,kT
    case("Tp")
      read(line,*) keyword1,Tp
    case("Tb")
      read(line,*) keyword1,Tb
    case("Pext")
      read(line,*) keyword1,Pext
!    case("inputfile")
!      read(line,*) keyword1,inputfile
    case("seed")
      read(line,*) keyword1,seed
      seed=-seed ! idum for ran1() needs to be negative
    case default
! an unknown word will stop the execution
      write(0,*) "Unknown keyword :",trim(keyword)
      stop
    end select
  end do
end subroutine read_input

function sinhx(x)
! series of tayle of sinh(x)/x
 implicit none
 real,  intent(in) :: x
 real :: sinhx
 real :: a2,a4,a6,a8,a10,x2
 a2=1.0/(2*3)
 a4=a2/(4*5)
 a6=a4/(6*7)
 a8=a6/(8*9)
 a10=a8/(10*11)
 x2=x**2
 sinhx=1.0+x2*(a2+x2*(a4+x2*(a6+x2*(a8+x2*a10))))
end function sinhx

function potential(q,L)
! potential energy
 implicit none
 real,  intent(in) :: q
 real,  intent(in) :: L
 real :: potential
 real :: pi

 pi=4.0*atan(1.0)

 potential=(1.0-cos(2.0*pi*q/L))*(0.5*L/pi)**2 
end function potential

function force(q,L)
! force
 implicit none
 real,  intent(in) :: q
 real,  intent(in) :: L
 real :: force
 real :: pi
 
 pi=4.0*atan(1.0)

 force=-0.5*L/pi*sin(2.0*pi*q/L)
end function force

function dUdV(q,L)
! derivate of potential energy with respect of length
 implicit none
 real,  intent(in) :: q
 real,  intent(in) :: L
 real :: dUdV
 real :: pi

 pi=3.141592653589793

 dUdV=0.5*L/(pi**2)*(1.0-cos(2.0*pi*q/L))-0.5*q/pi*sin(2.0*pi*q/L)
end function dUdV

end program nptlangevin
