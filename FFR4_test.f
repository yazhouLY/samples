      program ffrtest
      implicit none
      real*8,dimension(0:2000) :: F ! np波函数
      real*8 ::b(0:2000,2)
      real*8 :: FFR4 ! interpolation function
      real*8 :: Y
      integer :: i,j


      open(11,file='fort.7',status='old')
      do i=0,2000
      read(11,*)(b(i,j),j=1,2)
      F(i)=b(i,2)
      end do
      Y=4.5
      write(*,*)FFR4(Y,F,2001)

      end program ffrtest


      FUNCTION FFR4(Y,F,N)
            IMPLICIT REAL*8(A-H,O-Z)
            REAL*8 F(N),P,P1,P2,Q,X,FFR4
            REAL*8 Y
            PARAMETER(X=.16666666666667)
            P=Y
            I=P
            write(*,*)I
            IF(I.LE.0) GO TO 2
            IF(I.GE.N-2) GO TO 4
  1         P=P-I
            P1=P-1.
            P2=P-2.
            Q=P+1.
            FFR4=(-P2*F(I)+Q*F(I+3))*P*P1*X+(P1*F(I+1)-P*F(I+2))*Q*P2*.5
            RETURN
  2         IF(I.LT.0) GO TO 3
            I=1
            GO TO 1
  3         FFR4=F(1)
            RETURN
  4         IF(I.GT.N-2) GO TO 5
            I=N-3
            GO TO 1
  5         FFR4=F(N)
            RETURN
            END function



