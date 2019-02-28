
include "../fortran_lib/MailboxMod.f90"


program simulate_fortran

use MailboxMod

    implicit none
    integer :: i
    type(mbm_MailboxInfo) :: MI
    character(1024)  :: msg

    integer :: max_try, stat
    logical :: get_through

    max_try = 5
    stat = 0

    call mbm_setDefault(MI)

    call mbm_hello(MI, max_try)

    print *, "Hello finish."
    get_through = .true.

    do i = 1, 5

        print *, "CESM doing other stuff... ", i
        
        write (msg, "(A, I5)") "Step : ", i
        call sleep(3)

        print *, "Send message: ", msg
        call mbm_send(MI, msg, max_try, stat)
        if ( stat /= 0 ) then
            print *, "Cannot send message."
            get_through = .false.
            exit
        end if

        print *, "Receiving message ... "
        call mbm_recv(MI, msg, max_try, stat)
        if ( stat /= 0 ) then
            print *, "Cannot recv message."
            get_through = .false.
            exit
        end if


        print *, "Message received: ", trim(msg)

    end do


    msg = "<<END>>"
    call mbm_send(MI, msg)

    print *, "Program ends."

end program
