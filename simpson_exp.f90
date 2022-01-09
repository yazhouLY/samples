      program test1
      implicit none
      integer::i
      integer::nnu !Simpson积分小区间的数目
      real*8::xival !Simpson积分的积分区间长度
      real*8,allocatable,dimension(:)::rr !Simpson积分的格点
      real*8,allocatable,dimension(:)::rw !Simpson积分的权重系数
      real*8::f
      
      nnu=1000 !给出Simpson积分的小区间的数目
      xival=40.0 !给出Simpson积分的积分区间长度
      allocate(rr(0:nnu),rw(0:nnu))
      call simpson(nnu,xival,rr,rw)  !初始化高斯积分，输出积分格点与权重系数
      f=0.d0
      do i=0,nnu
            f=f+rr(i)**2*rw(i)
      end do
      write(*,*)f
      end program test1


 ! subroutine to generate the grid points and width for Simpson integration
       subroutine simpson(nnu,xival,rr,rw)
       implicit none 
       integer :: nxmx,nnu
       real*8 :: xival
       real*8,dimension(1:nnu) :: rr,rw
       real*8,dimension(1:nnu+1) :: sxx,wxx
       real*8 :: dx, d43,d23
       integer :: nxmx1, nxmx2,i 
       
        nxmx=nnu+1
        dx=xival/float(nxmx-1)
 
        d43=4.d0/3.d0
        d23=2.d0/3.d0
   
        wxx(1)=dx/3.d0
        wxx(nxmx)=dx/3.d0
        sxx(1)=0.d0
        sxx(nxmx)=float(nxmx-1)*dx

 
        nxmx1=nxmx-1
        nxmx2=nxmx-2
 
        do 50 i=2,nxmx1,2
        sxx(i)=float(i-1)*dx
  50    wxx(i)=d43*dx
 
        do 55 i=3,nxmx2,2
        sxx(i)=float(i-1)*dx
  55    wxx(i)=d23*dx
  
  
        do i=1,nnu
          rr(i)=sxx(i+1)
          rw(i)=wxx(i+1)
        end do 
       end subroutine 



