program finite
      implicit none
      interface 
      subroutine finitedifference(n,h,array1,array2)
            implicit none
            integer n
            real h
            real::array1(n),array2(n)
      end subroutine
      end interface
      
      
      real::b(2001,2)
      integer::i,j
      integer::n=2001
      real::array1(2001),array2(2001)

      open(11,file='fort.7',status='old')

      do i=1,2001
      read(11,*)(b(i,j),j=1,2)
      array1(i)=b(i,2)
      end do


     call finitedifference(n,0.01,array1,array2)


      do i=1,2001
      write(15,*)array2(i)
      end do

      end


!!!this subroutine is to obtain the second derivative
! with a forth-order accuracy
!!!caution:  there are arrays as parameter,so there should be
! an interface when calls it.
! n----------- number of grids
! h ---------- step length of small interval
! array1 ----- fundamental function
! array2 ----- second derivative of fundamental function
      subroutine finitedifference(n,h,array1,array2)
            implicit none
            integer::n
            real::h
            integer::i
            real::array1(n),array2(n)

            do i=1,1000
            array2(i)=(2.*array1(i)-5.*array1(i+1)+4.*array1(i+2)-1.*array1(i+3))/h**2
            end do
      
            do i=2001,1001,-1
            array2(i)=(2.*array1(i)-5.*array1(i-1)+4.*array1(i-2)-1.*array1(i-3))/h**2
            end do
      
       end subroutine
      
      
      




