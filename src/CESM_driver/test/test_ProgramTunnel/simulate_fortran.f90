
include "../../fortran_lib/ProgramTunnelMod.f90"


program simulate_fortran

use ProgramTunnelMod

    implicit none
    integer :: i
    type(ptm_TunnelSet) :: TS
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

    call ptm_setDefaultTunnelSet(TS)

    write (msg, "(I5)") n
    stat = ptm_sendText(TS, msg)

    print *, "Hello finish."
    get_through = .true.

    do i = 1, 5

        print *, "CESM doing other stuff... ", i
        
        write (msg, "(A, I5)") "Step : ", i
        call sleep(3)

        print *, "Send message: [", trim(msg), "]"
        stat = ptm_sendText(TS, msg)
        if ( stat .gt. 0 ) then
            print *, "Cannot send message."
            get_through = .false.
            exit
        end if
        

        stat = ptm_sendBinary(TS, "FLX", dat, n)
        print *, "Binary sent: ", dat


        print *, "Receiving message ... "
        stat = ptm_recvText(TS, msg)
        if ( stat .gt. 0 ) then
            print *, "Cannot recv message."
            get_through = .false.
            exit
        end if

        print *, "Message received: [", trim(msg), "]"
        
        stat = ptm_recvBinary(TS, "SST", dat, n)
        print *, "Binary received: ", dat


    end do


    msg = "<<END>>"
    stat = ptm_sendText(TS, msg)
    if ( stat .gt. 0 ) then
        print *, "Cannot send message."
    end if


    print *, "Program ends."

end program
