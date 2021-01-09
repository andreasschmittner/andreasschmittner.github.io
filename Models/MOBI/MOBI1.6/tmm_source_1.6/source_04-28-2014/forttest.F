      subroutine test(m,itr,arr)
      
      integer m, itr
      real*8 arr(m)
      
      integer i, j
      
C      write(*,*)m,n
C      do j=1,n
C        do i=1,m
C          write(*,*)i,j,arr(i,j)
C        enddo
C      enddo  

      write(*,*)m,itr
      do j=1,m
		write(*,*)j,arr(j)
      enddo  
      
      end