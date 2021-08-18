! llarray.f90
module M
   implicit none
!   integer, parameter :: dp = kind([double precision ::])
   type PTRllnode
      type(llnode), pointer :: p => NULL()
   end type PTRllnode
   type llnode
      integer col
!      complex(dp) payload
      type(PTRllnode) next
   end type llnode
end module M

program Q
   use M
   implicit none
   type(PTRllnode), pointer :: p
   type(PTRllnode), allocatable, target :: heads(:)
   integer nrow
   integer row
   integer col
   integer iunit, istat
!   complex(dp) payload

   open(newunit=iunit,file='csriinput.txt',status='old')
   read(iunit,*) nrow
   allocate(heads(nrow))
   do
      read(iunit,*,iostat=istat) row,col !,payload
      if(istat /= 0 .OR. row < lbound(heads,1) .OR. row > ubound(heads,1)) exit
      p => heads(row)
      do while(associated(p%p))
         p => p%p%next
      end do
      allocate(p%p)
      p%p%col = col
!      p%p%payload = payload
   end do

   do row = lbound(heads,1),ubound(heads,1)
      print '(*(g0))', 'row: ',row
      p => heads(row)
      do while(associated(p%p))
         print '(*(g0))', '   col: ',p%p%col
         p => p%p%next
      end do
   end do
end program Q