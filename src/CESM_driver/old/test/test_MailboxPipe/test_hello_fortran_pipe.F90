
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

    print *, "Program ends."

end program
