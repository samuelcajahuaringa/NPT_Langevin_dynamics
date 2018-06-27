program statistical_inefficiency
implicit none
integer :: i,j,k,l,m,n,nb_i,nb_f,tb_i,tb_f
integer :: ioerr,i_repeat,blk,stp1,stp2
integer :: nstep,trun,nblock,tblock
integer :: natoms
real(8) :: atoms
real(8), dimension(:), allocatable :: a ! stored data values (nstep)
real(8) :: average,variance,stddev,err,si
real(8) :: a_blk,a_run,a_avg,a_var,a_var_1,a_err
real(8) :: true_rms,true_err

nstep = 4194304
natoms = 1000
atoms = 1.0*natoms

allocate(a(1:nstep))

open(10,file="cell.dat",status="old",action="read")
do i=1,nstep
   read(10,*)a(i)
end do
close(10)

! Simple calculation of average and variance

open(10,file="si.dat")
1000 FORMAT (A40)
write(10,1000)"#tblock avg_per_atom rms_per_atom err_per_atom si"

a_avg   = SUM(a)/nstep             ! Sample average
a_var_1 = 0.0
do i=1,nstep
   a_var_1 = a_var_1 + (a(i)-a_avg)**2 
end do
a_var_1 = a_var_1/nstep         ! square variance
a_err   = SQRT(a_var_1 /nstep)  ! Error estimate neglecting any correlations

2000 FORMAT (I4,2X,F10.6,2X,F13.6,2X,F13.6,2X,F13.6)
write(10,2000)1,a_avg/atoms,sqrt(a_var_1/atoms),a_err/atoms,1.0

! We must take account of the correlations between successive values in time
! The two common methods which follow are actually very similar
! They differ in the choice of block lengths

! Traditional block analysis
! The rationale here is that 20 (say) independent block averages should be enough to get a reasonable
! estimate of the desired quantities, and that there is not a lot to gain by looking at more (shorter) blocks.
! Some workers just assume that the run may be divided into 10 or 20 blocks, which they hope will be independent.
! This is exactly what we do in our example programs, just for simplicity.
! We cannot recommend this in general, unless there is good reason to support the assumption of independence.
! If the 20 blocks are not independent, then attention should be focused on fewer (longer) blocks,
! rather than more (shorter) ones, and a plot of squared error estimate, or statistical inefficiency,
! vs 1/tblock carried out to extrapolate to tblock=nstep. The loop below provides the data for that plot.


!DO nblock = 20, 4, -1 ! Loop over number, and hence length, of blocks
tb_i = 2
tb_f = 2000
nb_i = nstep/tb_i
nb_f = nstep/tb_f

do tblock = tb_i,tb_f
   nblock = int(nstep/tblock)
   trun = nblock*tblock
   a_run = sum(a(1:trun))/trun
   a_var = 0.0

   do blk = 1,nblock
      stp1 = 1 + (blk-1)*tblock
      stp2 = blk*tblock
      a_blk = sum(a(stp1:stp2))/tblock
      a_var = a_var + (a_blk - a_run)**2 
   end do
   
   a_var = a_var/nblock
   si = tblock * a_var/ a_var_1
   a_err = sqrt(a_var/nblock)       

   write(10,2000)tblock,a_run/atoms,sqrt(a_var/atoms),a_err/atoms,si
end do

close(10)
deallocate(a)

end program statistical_inefficiency
