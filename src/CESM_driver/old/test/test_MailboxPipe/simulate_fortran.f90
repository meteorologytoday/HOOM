
include "../../fortran_lib/MailboxPipeMod.f90"


program simulate_fortran

use MailboxPipeMod

    implicit none
    integer :: i
    type(mbp_MailboxInfo) :: MI
    character(1024)  :: msg

    integer :: stat
    logical :: get_through

    stat = 0

    call mbp_setDefault(MI)
    call mbp_hello(MI)

    print *, "Hello finish."
    get_through = .true.

    do i = 1, 5

        print *, "CESM doing other stuff... ", i
        
        write (msg, "(A, I5)") "Step : ", i
        call sleep(3)

        print *, "Send message: [", trim(msg), "]"
        stat = mbp_send(MI, msg)
        if ( stat /= 0 ) then
            print *, "Cannot send message."
            get_through = .false.
            exit
        end if

        print *, "Receiving message ... "
        stat = mbp_recv(MI, msg)
        if ( stat /= 0 ) then
            print *, "Cannot recv message."
            get_through = .false.
            exit
        end if


        print *, "Message received: [", trim(msg), "]"

    end do


    msg = "<<END>>"
    stat = mbp_send(MI, msg)
    if ( stat /= 0 ) then
        print *, "Cannot send message."
    end if


    print *, "Program ends."

end program
