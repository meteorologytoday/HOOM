
include "../../fortran_lib/MailboxPipeMod2.f90"


program simulate_fortran

use MailboxPipeMod2

    implicit none
    integer :: i
    type(mbp_MailboxInfo) :: MI
    character(1024)  :: msg

    integer :: stat
    logical :: get_through

    real(8), pointer :: dat(:)
    integer :: n = 5

    allocate(dat(n))

    do i = 1, n
        dat(i) = i * 10
    end do

    stat = 0

    call mbp_setDefault(MI)
    !call mbp_hello(MI)

    write (msg, "(I5)") n
    stat = mbp_send_txt(MI, msg)

    print *, "Hello finish."
    get_through = .true.

    do i = 1, 5

        print *, "CESM doing other stuff... ", i
        
        write (msg, "(A, I5)") "Step : ", i
        call sleep(3)

        print *, "Send message: [", trim(msg), "]"
        stat = mbp_send_txt(MI, msg)
        if ( stat .gt. 0 ) then
            print *, "Cannot send message."
            get_through = .false.
            exit
        end if
        
        stat = mbp_send_bin(MI, dat, n)
        print *, "Sending binary: ", dat

        print *, "Receiving message ... "
        stat = mbp_recv_txt(MI, msg)
        if ( stat .gt. 0 ) then
            print *, "Cannot recv message."
            get_through = .false.
            exit
        end if

        print *, "Message received: [", trim(msg), "]"
        
        stat = mbp_recv_bin(MI, dat, n)
        print *, "Received binary: ", dat


    end do


    msg = "<<END>>"
    stat = mbp_send_txt(MI, msg)
    if ( stat .gt. 0 ) then
        print *, "Cannot send message."
    end if


    print *, "Program ends."

end program
