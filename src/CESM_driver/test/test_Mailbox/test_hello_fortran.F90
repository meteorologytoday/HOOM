
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

    print *, "Program ends."

end program
