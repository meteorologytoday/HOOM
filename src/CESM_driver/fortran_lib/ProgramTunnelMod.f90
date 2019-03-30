module ProgramTunnelMod
implicit none

integer, parameter :: c_send_txt = 1, c_recv_txt = 2, c_send_bin = 3, c_recv_bin = 4
character(len=256), parameter :: keys(4) = (/"X2Y_txt", "Y2X_txt", "X2Y_bin", "Y2X_bin"/)
integer, parameter :: fd_beg = 50


type ptm_Tunnel
    Integer :: next_idx
    Integer :: fds(2)
    character(len=256) :: fns(2)
end type


type ptm_TunnelSet
    type(ptm_Tunnel) :: tnls(4)
end type


contains

integer function ptm_get_file_unit()
    integer :: lu, iostat
    logical :: opened
      
    do lu = 999, 1,-1
       inquire (unit=lu, opened=opened, iostat=iostat)
       if (iostat.ne.0) cycle
       if (.not.opened) exit
    end do
    
    ptm_get_file_unit = lu
    return
end function 

subroutine ptm_makeFilename(filename, id, n)
    implicit none
    character(len=*) :: filename, id
    integer          :: n

    write (filename, '(A, A, A, I1, A)')  "_", trim(id), "_", n, ".fifo"

end subroutine

subroutine ptm_setDefault(TS)
    implicit none
    type(ptm_TunnelSet) :: TS
    integer :: i, j
    do i = 1, 4
        do j = 1, 2
            call ptm_makeFilename(TS%tnls(i)%fns(j), keys(i), j)
            TS%tnls(i)%fds(j) = fd_beg + (i-1) * 2 + (j-1)
            TS%tnls(i)%next_idx = 1
        end do
    end do
end subroutine 


subroutine ptm_printSummary(TS)
    implicit none
    type(ptm_TunnelSet) :: TS
    integer :: i, j

    do i = 1, 4
        print *, "keys(", i, ") => ", trim(keys(i))
    end do

    do i = 1, 4
        do j = 1, 2
            print *, trim(keys(i)), "(", TS%tnls(i)%fds(j), ") =>", trim(TS%tnls(i)%fns(j))

        end do
    end do

    print *, "Next sendText   uses idx: ", TS%tnls(c_send_txt)%next_idx
    print *, "Next sendBinary uses idx: ", TS%tnls(c_recv_txt)%next_idx
    print *, "Next recvText   uses idx: ", TS%tnls(c_send_bin)%next_idx
    print *, "Next recvBinary uses idx: ", TS%tnls(c_recv_bin)%next_idx

end subroutine

subroutine ptm_iterTunnel(TS, n)
    implicit none
    type(ptm_TunnelSet) :: TS
    integer :: tmp, n 

    !tmp = TS%tnls(n)%next_idx 
    TS%tnls(n)%next_idx = mod(TS%tnls(n)%next_idx, 2) + 1
    !print *, "What is the tunnel file you get? ", trim(TS%tnls(n)%fns(tmp))

end subroutine

integer function ptm_sendText(TS, msg)
    implicit none
    type(ptm_TunnelSet), target  :: TS
    character(len=*)     :: msg

    character(len=256), pointer :: fn
    type(ptm_Tunnel), pointer :: tnl
    integer :: idx, fd

    tnl => TS%tnls(c_send_txt)
    idx = tnl%next_idx
 
    fd =  tnl%fds(idx)   
    fn => tnl%fns(idx)

    call ptm_iterTunnel(TS, c_send_txt)

    ptm_sendText = 0
    open(unit=fd, file=fn, form="formatted", access="stream", action="write", iostat=ptm_sendText)
    if (ptm_sendText .gt. 0) then
        print *, "[ptm_sendText] Error during open."
        close(fd)
        return
    end if

    ptm_sendText = 0
    write (fd, *, iostat=ptm_sendText) msg
    if (ptm_sendText .gt. 0) then
        print *, "[ptm_sendText] Error during write."
        close(fd)
        return
    end if
   
    close(fd)
    
end function


logical function ptm_messageCompare(msg1, msg2)
    implicit none
    character(*) :: msg1, msg2

    if (msg1 .eq. msg2) then
        ptm_messageCompare = .true.
    else
        ptm_messageCompare = .false.
    end if

end function




end module 
