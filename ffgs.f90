module ffgs

    implicit none

    integer, parameter :: ffreal = 8
    integer, parameter :: norm_ell = 3000
    real(ffreal), parameter :: norm_tsz_fr = 143
    real(ffreal), parameter :: pi = 3.14159265359
    real(ffreal), dimension(2:10000) :: tsz_template, ksz_template

    type ffgs_params

        real(ffreal) :: aps_100, aps_143, aps_217, rps_100_143, rps_100_217, rps_143_217
        real(ffreal) :: acib_100, acib_143, acib_217, rcib_100_143, rcib_100_217, rcib_143_217, ncib
        real(ffreal) :: rszcib
        real(ffreal) :: atsz, aksz
        real(ffreal) :: fpol_100, fpol_143, fpol_217

    end type

contains

    subroutine init_ffgs(template_folder)

        character(len=*) :: template_folder
        real(ffreal), dimension(2:10000,2) :: temp
        integer :: l

        open(1,file=template_folder//'/tsz_143_eps0.50.dat')
        read(1,*) (temp(l,:), l=2, 10000)
        close(1)
        tsz_template = (/ (temp(l,2)/temp(norm_ell,2)*norm_ell*(norm_ell+1)/l/(l+1), l=2, 10000) /)

        open(1,file=template_folder//'/cl_ksz_148_trac.dat')
        read(1,*) (temp(l,:), l=2, 9999)
        close(1)
        ksz_template = (/ (temp(l,2)/temp(norm_ell,2)*norm_ell*(norm_ell+1)/l/(l+1), l=2, 10000) /)

    end subroutine

    function get_ffgs(params,fr1,fr2,x1,x2,lmin,lmax) result (fg)

        type(ffgs_params) :: params
        integer :: fr1, fr2, l
        character(len=1) :: x1, x2
        integer :: lmin, lmax
        real(ffreal), dimension(lmin:lmax) :: fg, ps, cib_template, tszcib_template

        fg = 0

        !Poisson
        if (fr1==100 .and. fr2==100) then
            ps = params%aps_100
        else if (fr1==143 .and. fr2==143) then
            ps = params%aps_143
        else if (fr1==217 .and. fr2==217) then
            ps = params%aps_217
        else if (fr1==100 .and. fr2==143) then
            ps = sqrt(params%aps_100 * params%aps_143) * params%rps_100_143
        else if (fr1==100 .and. fr2==217) then
            ps = sqrt(params%aps_100 * params%aps_217) * params%rps_100_217
        else if (fr1==143 .and. fr2==217) then
            ps = sqrt(params%aps_143 * params%aps_217) * params%rps_143_217
        end if
        if (x1=='E') then
            if (fr1==100) then
                ps = ps*params%fpol_100
            else if (fr1==143) then
                ps = ps*params%fpol_143
            else if (fr1==217) then
                ps = ps*params%fpol_217
            end if
        end if
        if (x2=='E') then
            if (fr2==100) then
                ps = ps*params%fpol_100
            else if (fr2==143) then
                ps = ps*params%fpol_143
            else if (fr2==217) then
                ps = ps*params%fpol_217
            end if
        end if
        fg = fg + ps


        if (x1=='T' .and. x2=='T') then
            !Clustered
            cib_template = (/((l/3000.)**params%ncib , l=lmin, lmax)/)
            if (fr1==100 .and. fr2==100) then
                fg = fg + params%acib_100 * cib_template
            else if (fr1==143 .and. fr2==143) then
                fg = fg + params%acib_143 * cib_template
            else if (fr1==217 .and. fr2==217) then
                fg = fg + params%acib_217 * cib_template
            else if (fr1==100 .and. fr2==143) then
                fg = fg + sqrt(params%acib_100 * params%acib_143) * params%rcib_100_143 * cib_template
            else if (fr1==100 .and. fr2==217) then
                fg = fg + sqrt(params%acib_100 * params%acib_217) * params%rcib_100_217 * cib_template
            else if (fr1==143 .and. fr2==217) then
                fg = fg + sqrt(params%acib_143 * params%acib_217) * params%rcib_143_217 * cib_template
            endif

            !tSZ
            if (fr1==100 .and. fr2==100) then
                fg = fg + params%atsz * tsz_template(lmin:lmax) * tszdep(100._ffreal,100._ffreal)
            else if (fr1==143 .and. fr2==143) then
                fg = fg + params%atsz * tsz_template(lmin:lmax) * tszdep(143._ffreal,143._ffreal)
            else if (fr1==100 .and. fr2==143) then
                fg = fg + params%atsz * tsz_template(lmin:lmax) * tszdep(100._ffreal,143._ffreal)
            end if

            !tSZ-CIB
            tszcib_template = sqrt(tsz_template(lmin:lmax) * cib_template)
            if (fr1==100 .and. fr2==100) then
                fg = fg + 2 * params%rszcib * sqrt((params%atsz * tszdep(100._ffreal,100._ffreal)) * (params%acib_100)) * tszcib_template
            else if (fr1==143 .and. fr2==143) then
                fg = fg + 2 * params%rszcib * sqrt((params%atsz * tszdep(143._ffreal,143._ffreal)) * (params%acib_143)) * tszcib_template
            else if (fr1==100 .and. fr2==143) then
                fg = fg + params%rszcib * ( sqrt((params%atsz * tszdep(100._ffreal,100._ffreal)) * (params%acib_143)) + sqrt((params%atsz * tszdep(143._ffreal,143._ffreal)) * (params%acib_100)) ) * tszcib_template
            else if (fr1==100 .and. fr2==217) then
                fg = fg + params%rszcib * sqrt((params%atsz * tszdep(100._ffreal,100._ffreal)) * (params%acib_217)) * tszcib_template
            else if (fr1==143 .and. fr2==217) then
                fg = fg + params%rszcib * sqrt((params%atsz * tszdep(100._ffreal,100._ffreal)) * (params%acib_143)) * tszcib_template
            end if

            !kSZ
            fg = fg + params%aksz * ksz_template(lmin:lmax)

        end if

        fg = fg * 2.*pi/norm_ell/(norm_ell+1)

    end function


    function tszdep(fr1,fr2)
        real(ffreal) :: fr1,fr2,x1,x2,x0,tszdep
        x1 = fr1/56.78
        x2 = fr2/56.78
        x0 = norm_tsz_fr/56.78
        tszdep =  (x1*(exp(x1)+1)/(exp(x1)-1) - 4) * (x2*(exp(x2)+1)/(exp(x2)-1) - 4) / (x0*(exp(x0)+1)/(exp(x0)-1) - 4)**2
    end function



end module
