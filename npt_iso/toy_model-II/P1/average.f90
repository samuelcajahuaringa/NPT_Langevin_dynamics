program true_values
implicit none
integer :: i,j,k,l,m,n,t,nb_i,nb_f,tb_i,tb_f
integer :: ioerr,i_repeat,blk,stp1,stp2,stp
integer :: nstep,trun,nblock,tblock
integer :: natoms
real(8) :: atoms
real(8), dimension(:), allocatable :: a ! stored data values (nstep)
real(8), dimension(:), allocatable ::sigma
real(8) :: average,variance,stddev,err,si
real(8) :: a_blk,a_blk2,a_run,a_avg,a_var,a_var_1,a_err
real(8) :: tcorr,true_rms,true_err
real(8) :: sigma_ave,sigma_sigma_ave,err_sigma_ave
!nstep = 4194304
!natoms = 1000
!atoms = 1.0*natoms
!allocate(a(1:nstep))

open(10,file="input_val",status="old",action="read")
read(10,*)nstep
read(10,*)natoms
read(10,*)tcorr
close(10)

atoms = 1.0*natoms
allocate(a(1:nstep))

open(10,file="cell.dat",status="old",action="read")
do i=1,nstep
   read(10,*)a(i)
end do
close(10)

! Simple calculation of average and variance

open(10,file="values.dat")
1000 FORMAT (A85)
write(10,1000)"#tcorr avg_per_atom err_per_atom sigma_per_atom var_sigma_per_atom err_sigma_per_atom"

a_avg   = SUM(a)/nstep             ! Sample average

2000 FORMAT (I4,2X,F10.6,2X,F13.6,2X,F13.6,2X,F13.6,2X,F13.6)

tblock = int(tcorr)
nblock = int(nstep/tblock)

allocate(sigma(tblock))

sigma = 0.0

trun = nblock*tblock

do t = 1,tblock
   a_blk = 0.0
   a_blk2 = 0.0

   do blk = 1,nblock
      stp = t + (blk-1)*tblock
      a_blk = a_blk + a(stp)
      a_blk2 = a_blk2 + a(stp)**2
   end do

   a_blk = a_blk/nblock
   a_var = a_blk2/nblock - a_blk**2  ! sigma2 
   sigma(t) = sqrt(a_var)
end do

a_blk = 0.0
a_blk2 = 0.0

sigma_ave = sum(sigma(1:tblock))/real(tblock)
true_rms = sigma_ave
true_err = true_rms/sqrt(1.0*tblock)

a_var = 0.0
do t = 1,tblock
   a_var = a_var + (sigma(t)-sigma_ave)**2
end do

sigma_sigma_ave = sqrt(a_var/tblock)
err_sigma_ave = sigma_sigma_ave/sqrt(1.0*tblock) 

write(10,2000)tblock,a_avg/atoms,true_err/atoms,sigma_ave/sqrt(atoms),sigma_sigma_ave/sqrt(atoms),err_sigma_ave/sqrt(atoms)

close(10)
deallocate(a)
deallocate(sigma)
end program true_values
