module MailboxPipeMod2
implicit none

type mbp_MailboxInfo

    Integer :: recv_txt_fd
    Integer :: send_txt_fd
    Integer :: recv_bin_fd
    Integer :: send_bin_fd
    
    character(len = 256) :: recv_txt_fn
    character(len = 256) :: send_txt_fn
    character(len = 256) :: recv_bin_fn
    character(len = 256) :: send_bin_fn


    character(len = 256) :: log_file
end type


contains

integer function mbp_get_file_unit()
    
    integer :: lu, iostat
    logical :: opened
      
    do lu = 999, 1,-1
       inquire (unit=lu, opened=opened, iostat=iostat)
       if (iostat.ne.0) cycle
       if (.not.opened) exit
    end do
    
    mbp_get_file_unit = lu
    return
end function 

subroutine mbp_setDefault(MI)
    implicit none
    type(mbp_MailboxInfo) :: MI

    MI%recv_txt_fn  = "_mymodel2cesm_txt.fifo"
    MI%send_txt_fn  = "_cesm2mymodel_txt.fifo"
 
    MI%recv_bin_fn  = "_mymodel2cesm_bin.fifo"
    MI%send_bin_fn  = "_cesm2mymodel_bin.fifo"
 
    MI%log_file = "log"

    MI%recv_txt_fd = mbp_get_file_unit()
    MI%send_txt_fd = mbp_get_file_unit()

    MI%recv_bin_fd = mbp_get_file_unit()
    MI%send_bin_fd = mbp_get_file_unit()
end subroutine 

subroutine mbp_appendPath(MI, path)
    implicit none
    type(mbp_MailboxInfo) :: MI
    character(len=256) :: path

    MI%recv_txt_fn  = path // "/" // MI%recv_txt_fn 
    MI%send_txt_fn  = path // "/" // MI%send_txt_fn
    MI%recv_bin_fn  = path // "/" // MI%recv_bin_fn 
    MI%send_bin_fn  = path // "/" // MI%send_bin_fn 
 
    MI%log_file = path // "/" // MI%log_file

end subroutine 

integer function mbp_recv_txt(MI, msg)
    implicit none
    type(mbp_MailboxInfo)  :: MI
    character(len=*)       :: msg
    
    logical :: file_exists

    mbp_recv_txt = 0
    open(unit=MI%recv_txt_fd, file=MI%recv_txt_fn, form="formatted", access="stream", action="read", iostat=mbp_recv_txt)
    if (mbp_recv_txt .gt. 0) then
        print *, "ERROR OPENING RECV TXT PIPE, errcode:", mbp_recv_txt
        close(MI%recv_txt_fd)
        return
    end if

    mbp_recv_txt = 0
    read (MI%recv_txt_fd, '(A)', iostat=mbp_recv_txt) msg
    if (mbp_recv_txt .gt. 0) then
        print *, msg
        print *, "ERROR READING RECV TXT PIPE, errcode:", mbp_recv_txt
        close(MI%recv_txt_fd)
        return
    end if

    close(MI%recv_txt_fd)
    
    msg = trim(msg)

end function


integer function mbp_send_txt(MI, msg)
    implicit none
    type(mbp_MailboxInfo)  :: MI
    character(len=*)       :: msg

    mbp_send_txt = 0
    open(unit=MI%send_txt_fd, file=MI%send_txt_fn, form="formatted", access="stream", action="write", iostat=mbp_send_txt)
    if (mbp_send_txt .gt. 0) then
        print *, "[mbp_send] Error during open."
        close(MI%send_txt_fd)
        return
    end if

    mbp_send_txt = 0
    write (MI%send_txt_fd, *, iostat=mbp_send_txt) msg
    if (mbp_send_txt .gt. 0) then
        print *, "[mbp_send] Error during write."
        close(MI%send_txt_fd)
        return
    end if
   
    close(MI%send_txt_fd)

end function



integer function mbp_recv_bin(MI, dat, n)
    implicit none
    type(mbp_MailboxInfo)  :: MI
    real(8), intent(inout) :: dat(n)
    integer, intent(in)    :: n
    integer                :: i


    mbp_recv_bin = 0
    open(unit=MI%recv_bin_fd, file=MI%recv_bin_fn, form="unformatted", &
         access="stream", action="read", iostat=mbp_recv_bin,  &
         convert='LITTLE_ENDIAN')

    if (mbp_recv_bin .gt. 0) then
        print *, "ERROR OPENING RECV BIN PIPE, errcode: ", mbp_recv_bin
        close(MI%recv_bin_fd)
        return
    end if

    mbp_recv_bin = 0
    read (MI%recv_bin_fd, iostat=mbp_recv_bin) (dat(i),i=1,n,1)
    if (mbp_recv_bin .gt. 0) then
        print *, "ERROR READING RECV BIN PIPE, errcode:", mbp_recv_bin
        close(MI%recv_bin_fd)
        return
    end if

    close(MI%recv_bin_fd)
    
end function


integer function mbp_send_bin(MI, dat, n)
    implicit none
    type(mbp_MailboxInfo)  :: MI
    real(8), intent(in)    :: dat(n)
    integer, intent(in)    :: n
    integer                :: i

    mbp_send_bin = 0
    open(unit=MI%send_bin_fd, file=MI%send_bin_fn, form="unformatted", &
         access="stream", action="write", iostat=mbp_send_bin,  &
         convert='LITTLE_ENDIAN')

    if (mbp_send_bin .gt. 0) then
        print *, "[mbp_send_bin] Error during open."
        close(MI%send_bin_fd)
        return
    end if

    mbp_send_bin = 0
    write (MI%send_bin_fd, iostat=mbp_send_bin) (dat(i), i=1,n,1)
    if (mbp_send_bin .gt. 0) then
        print *, "[mbp_send_bin] Error during write, err code: ", mbp_send_bin
        close(MI%send_bin_fd)
        return
    end if
   
    close(MI%send_bin_fd)

end function


subroutine mbp_hello(MI)
    implicit none
    type(mbp_MailboxInfo) :: MI
    character(256) :: msg

    integer :: stat

    stat = mbp_recv_txt(MI, msg)
    if (stat .gt. 0) then
        print *, "Something went wrong during recv stage. Error io stat: ", stat
        return
    end if

    if (mbp_messageCompare(msg, "<<TEST>>")) then
        print *, "Recv hello!"
    else
        print *, len(msg), " : ", len("<<TEST>>")
        print *, "Weird msg: [", msg, "]"
    end if

    stat = mbp_send_txt(MI, "<<TEST>>")
    if (stat .gt. 0) then
        print *, "Something went wrong during send stage. Error io stat: ", stat
        return
    end if

end subroutine

logical function mbp_messageCompare(msg1, msg2)
    implicit none
    character(*) :: msg1, msg2

    if (msg1 .eq. msg2) then
        mbp_messageCompare = .true.
    else
        mbp_messageCompare = .false.
    end if

end function


end module MailboxPipeMod2
