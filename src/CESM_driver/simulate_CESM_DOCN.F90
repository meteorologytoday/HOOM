include "./fortran_lib/field_tools.f90"
include "./fortran_lib/MailboxPipeMod.f90"

program CESM_SIMULATE
    use field_tools
    use MailboxPipeMod

    implicit none

    character(1024)       :: x_msg, x_fn, x_datetime_str
    type(mbp_MailboxInfo) :: x_MI
    integer :: x_w_fd, x_r_fd, x_curr_ymd, x_curr_tod
    integer :: x_stat, x_max_try

    real(8), pointer     :: x_hflx(:), x_swflx(:), x_taux(:), x_tauy(:)

    character(*), parameter :: x_F00 = "(a, '.ssm.', a, '.', a)" 

    real(8), pointer :: somtp(:)

    integer :: max_time = 3
    integer :: t
    integer :: lsize = 11600 

    allocate(somtp(lsize))

    x_curr_ymd = 10101
    x_curr_tod = 1800


    call docn_comp_init()
    do t = 1, max_time
        write (x_datetime_str, '(i0.8,A,i0.8)') x_curr_ymd, "-", x_curr_tod
        call docn_comp_run(x_datetime_str, t, 86400.0d0)
    end do
    call docn_comp_final()


CONTAINS

subroutine docn_comp_run(clock, t, dt)
    implicit none
    character(*) :: clock
    integer :: t
    real(8) :: dt

    integer :: n, stat
    logical :: firstcall

    
    if (t == 1) then
        firstcall = .true.
    else
        firstcall = .false.
    end if
 
    print *, "# WHAT TIME IS IT? ", trim(x_datetime_str)

    if (firstcall) then

        x_w_fd = mbp_get_file_unit()
        x_r_fd = mbp_get_file_unit()

        allocate(x_hflx(lsize))
        allocate(x_swflx(lsize))
        allocate(x_taux(lsize))
        allocate(x_tauy(lsize))

        do n = 1,lsize
            x_hflx(n)  = 0.0
            x_swflx(n) = 0.0
            x_taux(n)  = 0.0
            x_tauy(n)  = 0.0
        end do
         
        call mbp_setDefault(x_MI)

        x_msg = "MSG:INIT;CESMTIME:"//trim(x_datetime_str)//";"
        
        x_fn = "init_sst.bin"
        call write_1Dfield(x_w_fd, x_fn, x_hflx, lsize)
 
        x_fn = "SST_NEW.bin"
        call write_1Dfield(x_w_fd, x_fn, x_hflx, lsize)
        
        x_fn = "init_sst.bin"
        x_msg = trim(x_msg)//"SST:"//trim(x_fn)//";"

         
        call stop_if_bad(mbp_send(x_MI, x_msg), "INIT_SEND")
        
        print *, "Init msg sent: ", trim(x_msg), "."
        print *, "Now receiving..."
        call stop_if_bad(mbp_recv(x_MI, x_msg), "INIT_RECV")

        if (mbp_messageCompare(x_msg, x_fn) .neqv. .true.) then
            print *, "SSM init failed. Recive message: ", x_msg
            call shr_sys_abort ('SSM init failed.')
        end if
        
        call read_1Dfield(x_r_fd, trim(x_fn), somtp, lsize)
        !call mbp_delFile(trim(x_fn), x_r_fd)
        
        !do n = 1, lsize
        !  o2x%rAttr(kt,n) = somtp(n)
        !  o2x%rAttr(kq,n) = 0.0_R8
        !end do
      else

        x_msg = "MSG:RUN;CESMTIME:"//trim(x_datetime_str)//";"
        print *, "Not first call."

        ! ... Extract information from coupler ...

        x_fn = "HFLX.bin"
        x_msg = trim(x_msg)//"HFLX:"//trim(x_fn)//";"
        call write_1Dfield(x_w_fd, x_fn, x_hflx, lsize)

        x_fn = "SWFLX.bin"
        x_msg = trim(x_msg)//"SWFLX:"//trim(x_fn)//";"
        call write_1Dfield(x_w_fd, x_fn, x_swflx, lsize)

        x_fn = "TAUX.bin"
        x_msg = trim(x_msg)//"TAUX:"//trim(x_fn)//";"
        call write_1Dfield(x_w_fd, x_fn, x_taux, lsize)

        x_fn = "TAUY.bin"
        x_msg = trim(x_msg)//"TAUY:"//trim(x_fn)//";"
        call write_1Dfield(x_w_fd, x_fn, x_tauy, lsize)


        x_msg = trim(x_msg)//"SST_NEW:SST_NEW.bin;"
        write (x_msg, "(A, A, F10.2, A)") trim(x_msg), ";DT:", dt, ";"
        call stop_if_bad(mbp_send(x_MI, x_msg), "RUN_SEND")
 
        ! SSM is doing some MAGICAL calculation...
        
        call stop_if_bad(mbp_recv(x_MI, x_msg), "RUN_RECV")

        call read_1Dfield(x_r_fd, trim(x_msg), somtp, lsize)
        !call mbp_delFile(trim(x_msg), x_r_fd)
 
        ! ... doing ice formation stuff ... 

      endif

      print *, "SSM_AQUAP done."

end subroutine docn_comp_run

subroutine docn_comp_init()
implicit none


    print *, "DOCN_COMP_INIT"

end subroutine docn_comp_init


subroutine docn_comp_final()
implicit none

    print *, "DOCN_COMP_FINAL"
    x_msg = "MSG:END"
    call stop_if_bad(mbp_send(x_MI, x_msg), "FINAL")

end subroutine docn_comp_final

subroutine stop_if_bad(stat, stage)
    integer      :: stat
    character(*) :: stage

    if (stat > 0) then
          print *, 'MailBox error during stage ['//trim(stage)//']. Error state: ', stat
          call shr_sys_abort('MailBox error during stage ['//trim(stage)//']')
    end if
end subroutine stop_if_bad

subroutine shr_sys_abort(msg)
    implicit none
    character(*) :: msg

    print *, "ERROR: ", msg
    stop
end subroutine
end program CESM_SIMULATE
