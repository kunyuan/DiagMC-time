
  !==============Trace elapsed time ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE set_time_elapse
    implicit none
    !-- read and calculate time (in seconds) -------------------------
    call date_and_time(date, time, zone, tval)
    t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
    h_curr = tval(5)
    t_prev = t_curr
    h_prev = h_curr
    return
  END SUBROUTINE set_time_elapse
    

  !==============Trace elapsed time ==================================
  !! THIS IS PROJECT-INDEPENDENT 
  SUBROUTINE time_elapse
    implicit none
    
    !-- read and calculate time (in seconds) -------------------------
    call date_and_time(date, time, zone, tval)
    t_curr = tval(5)*3600.d0+tval(6)*60.d0+tval(7)+tval(8)*0.001d0 
    h_curr = tval(5)

    t_elap = t_curr-t_prev
    if(h_curr<h_prev) t_elap = t_elap+24*3600.d0
    t_prev = t_curr
    h_prev = h_curr 
    return
  END SUBROUTINE time_elapse
  !===================================================================
