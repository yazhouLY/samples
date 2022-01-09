      program gaussinterpolation
      implicit none
      real*8,dimension(0:2000) :: F ! np波函数纵坐标的值
      real*8,dimension(0:2000) :: xv ! np波函数横坐标的值 
      real*8 ::b(0:2000,2)
      real*8 :: fival ! interpolation function
      real*8 :: r
      real*8 :: alpha=1
      real*8 :: hcm =0.01
      integer :: i,j


      open(11,file='fort.7',status='old') !文件fort.7中储存了np波函数
      do i=0,2000
      read(11,*)(b(i,j),j=1,2)
      xv(i)=b(i,1) !把np波函数横坐标的值赋值到数组xv中
      F(i)=b(i,2) !把np波函数纵坐标的值赋值到数组F中
      end do
      
      r=4.5*hcm
      write(*,*)fival(r,xv,F,2001,alpha)

      end program gaussinterpolation

  !       REAL 4-point lagrange interpolation routine.
  !       interpolates thr FUNCTION value fival at point r from an
  !       array of points stored in fdis(ndm). this array is assumed
  !       to be defined such that the first element fdis(1) CONTAINS
  !       the FUNCTION value at r=xv(1) and xv(2 .. ndm) are monotonically
  !       increasing.
  !          r        a value of self variable
  !          xv       point of self variable(an array of one demension)
  !          fdis     value of function(an array of one demension)
  !          ndm      number of grids
  !          alpha    coefficient 
  !   ************************************************************************
            FUNCTION fival(r,xv,fdis,ndm,alpha)
            IMPLICIT REAL*8(A-H,O-Z)
            REAL*8 fdis(ndm),y1,y2,y3,y4
            DIMENSION xv(ndm)
            IF(r.GT.xv(ndm)) go to 9
            DO 5 k=1,ndm-2
 5          IF(r.LT.xv(k)) go to 6
            k=ndm-2
 6          nst=MAX(k-1,1)
            x1=xv(nst)
            x2=xv(nst+1)
            x3=xv(nst+2)
            x4=xv(nst+3)
            y1=fdis(nst+0)
            y2=fdis(nst+1)
            y3=fdis(nst+2)
            y4=fdis(nst+3)
            pii1=(x1-x2)*(x1-x3)*(x1-x4)
            pii2=(x2-x1)*(x2-x3)*(x2-x4)
            pii3=(x3-x1)*(x3-x2)*(x3-x4)
            pii4=(x4-x1)*(x4-x2)*(x4-x3)
            xd1=r-x1
            xd2=r-x2
            xd3=r-x3
            xd4=r-x4
            pi1=xd2*xd3*xd4
            pi2=xd1*xd3*xd4
            pi3=xd1*xd2*xd4
            pi4=xd1*xd2*xd3
            fival=y1*pi1/pii1+y2*pi2/pii2+y3*pi3/pii3+y4*pi4/pii4
            RETURN
 9          fival=fdis(ndm) * EXP(alpha*(xv(ndm)-r))
            RETURN
            END function
