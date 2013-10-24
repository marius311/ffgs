program test

    use ffgs

    type(ffgs_params) params

    params%aps_100 = 200
    params%aps_143 = 100
    params%aps_217 = 100
    params%rps_100_143 = 1
    params%rps_100_143 = 1
    params%rps_100_143 = 1

    params%acib_100 = 0
    params%acib_143 = 20
    params%acib_217 = 50
    params%rcib_100_143 = 1
    params%rcib_100_217 = 1
    params%rcib_143_217 = 1
    params%ncib = 0.8
    params%rszcib = -0.3

    params%atsz = 3
    params%aksz = 3

    params%fpol_100 = 0.05
    params%fpol_143 = 0.05
    params%fpol_217 = 0.05

    call init_ffgs("templates")

    print *, get_ffgs(params,143,143,'T','T',2,10)


end program
