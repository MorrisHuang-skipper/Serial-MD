!> module for array of linked list
module llarray
    implicit none
    type container
        type(llnode), pointer :: p => NULL()
    end type container
    type llnode
        integer value
        type(container) next
    end type llnode
end module llarray

!> program that test for module llarray
! program lltest
!     use llarray
!     implicit none
!     integer :: i, j, n
!     type(container), pointer :: tmp, current, next
!     type(container), allocatable, target :: heads(:)
!     n = 0
!     do 
!         allocate(heads(1000))
!         do i = 1, 1000
!             do j = 1, 1000
!                 tmp => heads(i) 
!                 do while (associated(tmp%pnt))
!                     tmp => tmp%pnt%next 
!                 end do
!                 allocate(tmp%pnt)
!                 tmp%pnt%value=i*j
!             end do
!         end do

!         ! > delete linked list
!         do i = 1, 1000
!             current => heads(i)
!             do while (associated(current % pnt))
!                 next => current % pnt % next
!                 if (associated(next % pnt)) then
!                     deallocate(current % pnt)
!                     nullify(current % pnt)
!                 end if
!                 current => next
!             end do
!         end do
!         deallocate(heads)
!         call sleep(3)
!         print *, n
!         n = n + 1
!     end do
! end program lltest

!> test 2
! program lltest
!     use llarray
!     implicit none
!     integer :: i, j
!     type(container), pointer :: tmp, newhead
!     type(container), allocatable, target :: heads(:)
!     allocate(heads(10))

!     ! do i = 1, 1
!     !     do j = 1, 2
!     !         tmp => heads(i) 
!     !         do while (associated(tmp%pnt))
!     !             tmp => tmp%pnt%next 
!     !         end do
!     !         allocate(tmp%pnt)
!     !         tmp%pnt%value=i*j
!     !         tmp => tails(i)
!     !     end do
!     ! end do

!     ! do i = 1, 5
!     !     print *, 'i', i
!     !     tmp => heads(i)
!     !     do while(associated(tmp%pnt))
!     !         print *, tmp%pnt%value
!     !         tmp => tmp%pnt%next
!     !     end do
!     ! end do
! end program lltest

! program lltest
!     ! > simple linked list
!     integer :: i
!     type LLnode
!         integer :: value
!         type(LLnode), pointer :: next => NULL()
!     end type LLnode
!     type(LLnode), pointer :: head, tmp
!     do i = 1, 10
!         allocate(tmp)
!         tmp % value = i
!         tmp % next => head
!         head => tmp
!         nullify(tmp)
!     end do

!     do i = 1, 10
!         print*, head % value
!         head => head % next
!     end do
! end program lltest