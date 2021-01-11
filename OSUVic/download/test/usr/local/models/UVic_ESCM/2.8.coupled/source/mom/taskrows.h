!====================== include file "taskrows.h" ======================

!     max_tasks  = maximum number of processors
!     num_processors = requested number of processors (each gets a task)
!     jstask         = first row accessed by a task
!     jetask         = last row accessed by a task
!     num_crpt       = number of calculated rows per task
!     pn             = processor number

      parameter (max_tasks = 2048)
      common /taskrows/ num_processors
      common /taskrows/ jstask(max_tasks), jetask(max_tasks)
      common /taskrows/ num_rcpt(max_tasks)
#if defined coarse_grained_parallelism

!     for CRAY

!DIR$ TASKCOMMON procsi
#endif
      integer pn
      common /procsi/ pn
