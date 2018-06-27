module routines
implicit none
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
!    case("W")
!      read(line,*) keyword1,W
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
  !if(inputfile=="") then
  !  write(0,*) "Specify input file"
  !  stop
  !end if
end subroutine read_input

subroutine initialize(M,etai,vetai,Gi,Qi,kT,tau,Nf,seed)
  implicit none
  integer, intent(in) :: M
  real, intent(inout) :: etai(M)
  real, intent(inout) :: vetai(M)
  real, intent(inout) :: Gi(M)
  real, intent(inout) :: Qi(M)
  real,    intent(in) :: kT
  real,    intent(in) :: tau
  integer, intent(in) :: Nf
  integer, intent(in) :: seed 
  real, external :: gasdev
  integer :: k 

  Qi(1) = Nf*kT*tau**2

  do k=2,M 
     Qi(k) = kT*tau**2
  end do

  do k=1,M
     etai(k) = 0.0
     vetai(k) = 0.0 !sqrt(kT/Qi(k))*gasdev(seed)
     if (k > 1) then
     Gi(k) = (Qi(k-1)*vetai(k-1)**2 - kT)/Qi(k)
     end if
  end do
  Gi(1) = 0.0
end subroutine initialize

subroutine nhc_thermostat(dt,engkin,Nf,kT,M,etai,vetai,Gi,Qi,tscale)
  implicit none
  real, intent(in)    :: dt
  real, intent(in)    :: engkin
  integer, intent(in)    :: Nf
  real, intent(in)    :: kT
  integer, intent(in)    :: M
  real, intent(inout) :: etai(M)
  real, intent(inout) :: vetai(M)
  real, intent(inout) :: Gi(M)
  real, intent(in)    :: Qi(M)
  real, intent(out)   :: tscale 
  integer :: k
  real :: aa
  real :: ekin
  real :: dts
  real :: dt2
  real :: dt4
  real :: dt8
  ekin=2.0*engkin
  tscale=1.0
  dts=dt
  dt2=0.5*dt  ! 0.5*dt
  dt4=0.25*dt ! 0.25*dt
  dt8=0.125*dt ! 0.125*dt

  Gi(1) = (ekin-Nf*kT)/Qi(1)

  vetai(M)=vetai(M)+dt4*Gi(M)

  do k=M-1,2,-1
     aa=exp(-dt8*vetai(k+1))
     vetai(k)=vetai(k)*aa
     vetai(k)=vetai(k)+dt4*Gi(k)
     vetai(k)=vetai(k)*aa
  end do

  aa=exp(-dt8*vetai(2))
  vetai(1)=vetai(1)*aa
  vetai(1)=vetai(1)+dt4*Gi(1)
  vetai(1)=vetai(1)*aa

  tscale=tscale*exp(-dt2*vetai(1))
     
  Gi(1)=(tscale*tscale*ekin-Nf*kT)/Qi(1)

  do k=1,M
     etai(k)=etai(k)+dt2*vetai(k)
  end do

  vetai(1)=vetai(1)*aa
  vetai(1)=vetai(1)+dt4*Gi(1)
  vetai(1)=vetai(1)*aa   
                       
  do k=2,M-1 
     aa=exp(-dt8*vetai(k+1))
     vetai(k)=vetai(k)*aa
     Gi(k)=(Qi(k-1)*vetai(k-1)**2 - kT)/Qi(k)  
     vetai(k)=vetai(k)+dt4*Gi(k)
     vetai(k)=vetai(k)*aa
  end do
 
  Gi(M)=(Qi(M-1)*vetai(M-1)**2 - kT)/Qi(M)
  vetai(M)=vetai(M)+dt4*Gi(M)    
     
end subroutine nhc_thermostat

subroutine hamilton(engkin,engpot,P,Li,bengkin,Nf,kT,M1,Qi,etai,vetai,M2,Qb,etab,vetab,Hi)
! energy of system
  implicit none
  real,    intent(in)  :: engkin
  real,    intent(in)  :: engpot
  real,    intent(in)  :: P
  real,    intent(in)  :: Li
  real,    intent(in)  :: bengkin 
  integer, intent(in)  :: Nf
  real,    intent(in)  :: kT
  integer, intent(in)  :: M1
  real,    intent(in)  :: Qi(M1)
  real,    intent(in)  :: etai(M1)
  real,    intent(in)  :: vetai(M1)
  integer, intent(in)  :: M2
  real,    intent(in)  :: Qb(M2)
  real,    intent(in)  :: etab(M2)
  real,    intent(in)  :: vetab(M2)
  real,    intent(out) :: Hi 
  integer :: i
  Hi=0.0
  Hi=Hi + engkin + engpot + P*li + bengkin 
  Hi=Hi + real(Nf)*kT*etai(1)
  do i=2,M1
     Hi=Hi+kT*etai(i)
  end do

  do i=1,M1
     Hi=Hi+0.5*Qi(i)*vetai(i)**2
  end do
  
  do i=1,M2
     Hi=Hi+kT*etab(i)+0.5*Qb(i)*vetab(i)**2
  end do
  
end subroutine hamilton

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

end module routines
