       subroutine icedata
!=======================================================================
!     read and interpolate ice data
!     for a specific data set that gives ice thickness (relative to
!     present day) and area every 1000 years from -19000 to 2000
!     calendar years (21kbp to 0kbp)

!     based on code by: M. Eby
!=======================================================================
#if defined uvic_embm && defined uvic_embm_icedata

      implicit none

      integer i, iou, j, n, ln, ib(10), ic(10)

      logical first_time, intrp, exists, inqvardef

      real wt3, wt1, c100, yrl(3), iyr(3)

      real, allocatable :: time(:)

      save time, ln, yrl, first_time

# include "param.h"
# include "atm.h"
# include "cembm.h"
# include "ice.h"
# include "levind.h"
# include "tmngr.h"

      character(120) :: fname, name, new_file_name

      c100 = 100.

!     recalculate land masks to be consistent with aicel
      fname = new_file_name ("tmsk.nc")
      inquire (file=trim(fname), exist=exists)
      if (exists) then
        call openfile (fname, iou)
        exists = inqvardef('tmsk', iou)
        ib(:) = 1
        ic(:) = imt
        ic(2) = jmt
        if (exists) call getvara ('tmsk', iou, imt*jmt, ib, ic, tmsk
     &,   c1, c0)
        call closefile (iou)
      endif
      if (.not. exists) then
        tmsk(:,:) = 0.
        do j=1,jmtm1
          do i=1,imtm1
            if (kmt(i,j) .gt. 0.) tmsk(i,j) = 1.
          enddo
        enddo
      endif

      name = "icedata.nc"

      if (.not. allocated (time)) then
        fname = new_file_name (name)
        inquire (file=trim(fname), exist=exists)
        if (.not. exists) then
          print*, "Warning => ", trim(fname), " does not exist."
          ln = 3
          allocate ( time(ln) )
          time(:) = year0
          aicel(:,:,:) = 0.
          hicel(:,:,:) = 0.
          first_time = .false.
        else
          call openfile (fname, iou)
          call getdimlen ('time', iou, ln)
          allocate ( time(ln) )
          ib(:) = 1
          ic(:) = ln
          call getvara ('time', iou, ln, ib, ic, time, c1, c0)
          first_time = .true.
          call closefile (iou)
        endif
        iyr(:) = 0
        yrl(:) = 0.
      else
        first_time = .false.
      endif

# if defined uvic_embm_icedata_transient
      yrl(2) = min(time(ln), max(time(1), year0 + relyr))
# else
      yrl(2) = min(time(ln), max(time(1), ice_yr))
# endif

      intrp = .false.
      if (yrl(2) .gt. time(1) .and. yrl(2) .lt. time(ln)) intrp = .true.

      if (first_time .or. yrl(2) .gt. yrl(3)) then
!       read data
        ib(:) = 1
        ic(:) = imt
        ic(2) = jmt
        ic(3) = 1
        fname = new_file_name (name)
        if (intrp) then
          do n=2,ln
            if (time(n-1) .le. yrl(2) .and. time(n) .ge. yrl(2)) then
              yrl(1) = time(n-1)
              iyr(1) = n-1
              yrl(3) = time(n)
              iyr(3) = n
            endif
          enddo
          call openfile (fname, iou)
          ib(3) = iyr(1)
          print*, "=> reading ice data for year:",yrl(1)
          call getvara ('hicel', iou, imt*jmt, ib, ic
     &,     hicel(1,1,1), c100, c0)
          call getvara ('aicel', iou, imt*jmt, ib, ic
     &,     aicel(1,1,1), c1, c0)
          call closefile (iou)
          call embmbc (hicel(1,1,1))
          call embmbc (aicel(1,1,1))
          call openfile (fname, iou)
          ib(3) = iyr(3)
          print*, "=> reading ice data for year:",yrl(3)
          call getvara ('hicel', iou, imt*jmt, ib, ic
     &,     hicel(1,1,3), c100, c0)
          call getvara ('aicel', iou, imt*jmt, ib, ic
     &,     aicel(1,1,3), c1, c0)
          call closefile (iou)
          call embmbc (hicel(1,1,3))
          call embmbc (aicel(1,1,3))
        else
          if (yrl(2) .le. time(1)) then
            n = 1
            yrl(1) = time(1)
            yrl(3) = time(1)
            iyr(n) = 1
          else
            n = 3
            yrl(1) = time(ln)
            yrl(3) = time(ln)
            iyr(n) = ln
          endif
          call openfile (fname, iou)
          ib(3) = iyr(n)
          print*, "=> reading ice data for year:",yrl(2)
          call getvara ('hicel', iou, imt*jmt, ib, ic
     &,     hicel(1,1,2), c100, c0)
          call getvara ('aicel', iou, imt*jmt, ib, ic
     &,     aicel(1,1,2), c1, c0)
          call closefile (iou)
          call embmbc (hicel(1,1,2))
          call embmbc (aicel(1,1,2))
          hicel(:,:,1) = hicel(:,:,2)
          hicel(:,:,3) = hicel(:,:,2)
          aicel(:,:,1) = aicel(:,:,2)
          aicel(:,:,3) = aicel(:,:,2)
        endif
      endif

      if (intrp) then
!       interpolate data
        wt1 = 1.
        if (yrl(3) .ne. yrl(1)) wt1 = (yrl(3)-yrl(2))/(yrl(3)-yrl(1))
        wt1 = max(0., min(1., wt1))
        wt3 = 1. - wt1
        do j=1,jmt
          do i=1,imt
            hicel(i,j,2) = hicel(i,j,1)*wt1 + hicel(i,j,3)*wt3
            aicel(i,j,2) = aicel(i,j,1)*wt1 + aicel(i,j,3)*wt3
            if (aicel(i,j,2) .lt. 0.5) then
              aicel(i,j,2) = 0.
              hicel(i,j,2) = 0.
            else
              aicel(i,j,2) = 1.
            endif
          enddo
        enddo
      endif
      call embmbc (hicel(1,1,2))
      call embmbc (aicel(1,1,2))

      do j=1,jmtm1
        do i=1,imtm1
          if (aicel(i,j,2) .ge. 0.5) tmsk(i,j) = 0.
        enddo
      enddo
      call embmbc (tmsk)
#endif

      return
      end
