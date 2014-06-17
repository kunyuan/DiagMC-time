
SUBROUTINE test_subroutine
    implicit none
    !integer :: isamp

    !========  test drawing subroutine =====================
    call initialize_markov
    call print_config

    !======== analytic_integration =========================
    call read_GWGamma
    call calculate_Gam1
    call output_Gam1

    return
END SUBROUTINE
