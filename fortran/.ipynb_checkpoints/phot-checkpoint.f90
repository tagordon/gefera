module phot
use iso_c_binding
use ellip

implicit none

!real*8, parameter :: pi = 4.d0 * Atan(1.d0), pihalf = 2.d0 * Atan(1.d0)
real*8, parameter :: twopithree = 8.d0 * Atan(1.d0) / 3.d0, twopi = 8.d0 * Atan(1.d0)
real*8, parameter :: o3 = 1.d0 / 3.d0, o9 = 1.d0 / 9.d0
real*8, parameter :: pithird = 4.d0 * Atan(1.d0) / 3.d0, pisixth = 4.d0 * Atan(1.d0) / 6.d0

contains

subroutine phis(rp, rm, bp, bm, bpm, cth, sth, pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)

    real*8 :: rp, rm, bp, bm, bpm, cth, sth
    
    ! intersection angle from planet center, intersection from moon center, same 
    ! values relative to bp and bm vectors respectively
    real*8 :: pp, pm, pp1, pm1, pp2, pm2
    
    ! derivatives 
    real*8 :: pp_rp, pp_rm, pp_bpm ! pp_theta = 1
    real*8 :: pm_rp, pm_rm, pm_bpm, thetam_theta, thetam_bp, thetam_bpm
    
    ! Four times the area of the triangle formed by rm, rp, and bpm
    real*8 :: delta
    ! Variables used in sorting the sides of the triangle
    real*8 :: a, b, c, tmp
    
    ! angle between bpm vector and bp vector, 
    ! angle between bpm vector and bm vector
    real*8 :: theta, thetam
    
    ! for avoiding divisions
    real*8 :: denom, obm
    
    thetam = Atan2(bp * sth, bpm - bp * cth)
    theta = Atan2(sth, cth)
    obm = 1.d0 / bm
    
    ! find 4 * area of triangle using modified Heron's formula 
    a = rp
    b = rm
    c = bpm
    if (c .gt. b) then
        tmp = c
        c = b
        b = tmp
    end if
    if (b .gt. a) then
        tmp = b
        b = a
        a = tmp
    end if
    delta = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    denom = 1.d0 / (delta * bpm * rm * rp)
    
    pm = Atan2(delta, (rm - rp) * (rm + rp) + bpm * bpm)   
    pm_bpm = ((rm + bpm) * (rm - bpm) - rp * rp) * denom * rm * rp
    pm_rp = 2 * rp * denom * bpm * rm * rp
    pm_rm = ((bpm - rm) * (bpm + rm) - rp * rp) * denom * rp * bpm
    
    pp = Atan2(delta, (rp - rm) * (rp + rm) + bpm * bpm)
    !pp = Asin(rm * Sin(pm) / rp)
    pp_bpm = ((rp - rm) * (rp + rm) - bpm * bpm) * denom * rm * rp
    pp_rp = ((bpm - rp) * (bpm + rp) - rm * rm) * denom * rm * bpm
    pp_rm = 2 * rm * denom * rm * rp * bpm
    
    thetam_bp = bpm * sth * obm * obm
    thetam_theta = ((bpm - bm) * (bpm + bm) - bp * bp) * 0.5 * obm * obm
    thetam_bpm = -bp * sth * obm * obm

    pm1 = thetam + pm
    pm2 = thetam - pm
    pp1 = theta + pp
    pp2 = theta - pp
    
    if (pm1 .gt. pi) then
        pm1 = pm1 - twopi
    end if
    if (pp1 .gt. pi) then
        pp1 = pp1 - twopi
    end if

end

subroutine kappas_p(rp, bp, kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)

    ! kp = angle to intersection from center of planet, 
    ! kps = angle to intersection from center of star 
    real*8 :: rp, bp, kp, kps
    
    ! derivatives 
    real*8 :: kp_rp, kp_bp, kps_rp, kps_bp
    
    ! variables used in sorting sides of triangle
    real*8 :: a, b, c
    
    ! four times the area of the triangle with sides rp, bp, and 1
    real*8 :: delta
    real*8 :: denom
    
    if (bp .gt. 1.d0) then
        a = bp
        b = 1.d0
        c = rp
    else
        a = 1.d0
        b = bp
        c = rp
        if (rp .gt. bp) then
            b = rp
            c = bp
        end if
    end if
    delta = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    denom = 1.d0 / (delta * bp * rp)
    
    kps = Atan2(delta, (1.d0 - rp) * (1.d0 + rp) + bp * bp)
    kps_bp = ((1.d0 - bp) * (1.d0 + bp) - rp * rp) * rp * denom
    kps_rp = 2 * rp * rp * bp * denom
    
    kp = Atan2(delta, (rp - 1.d0) * (rp + 1.d0) + bp * bp)
    kp_bp = ((rp + bp) * (rp - bp) - 1.d0) * rp * denom
    kp_rp = ((bp + rp) * (bp - rp) - 1.d0) * bp * denom
end 

subroutine kappas_m(rm, bp, bm, bpm, cth, sth, km, kms, km_rm, km_bp, km_bpm, km_theta, &
                    kms_rm, kms_bp, kms_bpm, kms_theta)

    ! km = angle to interection from center of moon, 
    ! kms = angle to intersection from center of planet
    real*8 :: rm, bp, bm, bpm, cth, sth, km, kms
    
    ! derivatives
    real*8 :: km_rm, km_bp, km_bpm, km_theta, kms_rm, kms_bp, kms_bpm, kms_theta
    
    ! variables used in sorting sides of triangle
    real*8 :: a, b, c
    
    ! four times the area of the triangle with sides rm, bm, and 1
    real*8 :: delta
    
    ! some useful quantities
    real*8 :: denom, xs, xm, yp, ypm, ytheta
    
    xs = (1.d0 - bm) * (1.d0 + bm) - rm * rm
    xm = (rm - bm) * (rm + bm) - 1.d0
    yp = bp - bpm * cth
    ypm = bpm - bp * cth
    ytheta = bp * bpm * sth
    
    if (bm .ge. 1.d0) then
        a = bm
        b = 1.d0
        c = rm
    else
        a = 1.d0
        b = bm
        c = rm
        if (rm .ge. bm) then
            b = rm
            c = bm
        end if
    end if
    delta = Sqrt((a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c)))
    denom = 1.d0 / (delta * bm * bm)
    
    km = Atan2(delta, (rm - 1.d0) * (rm + 1.d0) + bm * bm)
    kms = Atan2(delta, (1.d0 - rm) * (1.d0 + rm) + bm * bm)
    
    km_rm = ((bm + rm) * (bm - rm) - 1.d0) / (rm * delta)
    kms_rm = 2 * rm * bm * bm * denom
    
    km_theta = ytheta * xm * denom
    kms_theta = ytheta * xs * denom
    
    km_bpm = ypm * xm * denom
    kms_bpm = ypm * xs * denom
    
    km_bp = yp * xm * denom
    kms_bp = yp * xs * denom
end 

subroutine bm_x(bp, bm, bpm, cth, sth, dbm)

    real*8 :: bp, bm, bpm, cth, sth
    real*8, dimension(3) :: dbm
    real*8 :: obm 
    
    obm = 1.d0 / bm
    dbm(1) = (bp - bpm * cth) * obm
    dbm(2) = (bpm - bp * cth) * obm
    dbm(3) = bp * bpm * sth * obm
end 

! main loop to compute the flux at each timestep by finding the correct geometry and
! calling the integration routines 
subroutine flux(c1, c2, rp, rm, bp, bpm, cth, sth, lc, j) bind(C, name="flux")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: rp, rm
    real (c_double), bind(C), dimension(j) :: bp, cth, sth, bpm
    real (c_double), bind(C), intent(out), dimension(8, j) :: lc
    real*8, dimension(8) :: f0
    real*8, dimension(3) :: ld
    real (c_double), bind(C) :: c1, c2
    real*8 :: of0
    
    ! half angles: from planet to planet/star intersection, moon to moon/star
    ! intersection, spanned by planet on limb of star, spanned by moon on 
    ! limb of star
    real*8 :: kp, km, kps, kms
        
    ! derivatives of above angles
    real*8 :: kp_rp, kp_bp, kps_rp, kps_bp
    real*8 :: km_rm, km_bp, km_bpm, km_theta
    real*8 :: kms_rm, kms_bp, kms_bpm, kms_theta
    
    ! angles to planet-moon intersection from planet center 
    !relative to bp vector (and derivatives)
    real*8 :: pp1, pp2
    real*8 :: pp_rp, pp_rm, pp_bpm ! pp_theta = 1
    
    ! angles to planet-moon intersection from moon center 
    ! relative to bm vector (and derivatves)
    real*8 :: pm1, pm2
    real*8 :: pm_rp, pm_rm, thetam_bp, pm_bpm, thetam_theta
    ! derivative of angle between bpm and bm vector with respect to bpm
    real*8 :: thetam_bpm
    
    ! used to determine cases for three body overlaps, might not be needed. 
    ! Check if some of these (costheta, cosphi) can be removed when optimizing things later 
    real*8 :: phi, phi_bpm, phi_bp, phi_theta, d1, d2, delta, a, b, c, tmp
    
    ! For chain rule stuff
    real*8, dimension(j) :: bm
    real*8, dimension(3) :: dbm, dbm0
    
    bm = Sqrt((bp - bpm)**2.d0 + 2 * bp * bpm * (1.d0 - cth))
    dbm0 = 0.d0
    
    ld(1) = 1.d0 - c1 - 2 * c2
    ld(2) = c1 + 2 * c2
    ld(3) = c2
    
    ! normalization factors 
    f0(1) = ld(1) * pi + ld(2) * twopithree + ld(3) * pihalf
    f0(2) = 0.d0
    f0(3) = 0.d0
    f0(4) = 0.d0
    f0(5) = 0.d0
    f0(6) = 0.d0
    f0(7) = -pi + twopithree 
    f0(8) = -2 * pi + 2 * twopithree + pihalf
    
    of0 = 1.d0 / f0(1)
    
    do i=1,j,1
    
        lc(:, i) = f0 * of0
        
        if (rp .gt. 1.d0 + bp(i)) then
            lc(:, i) = 0.d0
        else if (rm .gt. 1.d0 + bm(i)) then 
            lc(:, i) = 0.d0
        else if (bpm(i) .gt. rp + rm) then
            ! moon and planet don't overlap each other 
            if (bp(i) .lt. 1.d0 - rp) then
                ! planet completely inside star
                lc(:, i) = lc(:, i) - 2 * Fcomplete(ld, rp, bp(i), dbm0, .TRUE.) * of0
            else if (bp(i) .lt. 1.d0 + rp) then
                ! planet partially overlaps star
                call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                lc(:, i) = lc(:, i) - 2 * (Fstar(ld, kps, kps_rp, 0.d0, kps_bp, 0.d0, 0.d0) &
                                    + F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                            dbm0, .TRUE., .TRUE.)) * of0
            end if
            if (bm(i) .lt. 1.d0 - rm) then
                ! moon completely inside star
                call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                lc(:, i) = lc(:, i) - 2 * Fcomplete(ld, rm, bm(i), dbm, .FALSE.) * of0
            else if (bm(i) .lt. 1.d0 + rm) then
                ! moon partially overlaps star
                call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                              km_rm, km_bp, km_bpm, km_theta, &
                              kms_rm, kms_bp, kms_bpm, kms_theta)
                call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                lc(:, i) = lc(:, i) - 2 * (Fstar(ld, kms, 0.d0, kms_rm, kms_bp, kms_bpm, kms_theta) &
                                    + F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                            dbm, .FALSE., .TRUE.)) * of0
            end if
        else
            ! moon and planet do overlap each other 
            if (bp(i) .gt. rp + 1.d0) then
                if (bm(i) .gt. rm + 1.d0) then
                    ! neither moon nor planet overlap star
                    lc(:, i) = f0 * of0
                else
                    ! moon partially overlaps star, planet does not overlap star
                    call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                    call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                    lc(:, i) = 2 * (Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                    - F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                            dbm, .FALSE., .TRUE.)) * of0

                end if
            else
                if (bm(i) .gt. rm + 1.d0) then
                    ! planet partially overlaps star, moon is outside of star
                    call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                    lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                    - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                        dbm0, .TRUE., .TRUE.)) * of0
                else
                    if (bp(i) + rp .le. 1.d0) then
                        if (bm(i) + rm .le. 1.d0) then
                            if (bpm(i) + rm .le. rp) then
                                ! moon and planet both overlap star, moon fully overlapped by planet
                                lc(:, i) = (f0 - 2 * Fcomplete(ld, rp, bp(i), dbm0, .TRUE.)) * of0
                            else
                                ! Case E
                                ! bookmark
                                ! moon and planet both overlap star, moon and planet partially overlap each other 
                                call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                call phis(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                          pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                lc(:, i) = (f0 - Arc(ld, pp1, pp2, rp, bp(i), &
                                                     pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                     -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                     dbm0, .TRUE., .FALSE., .FALSE.) &
                                                - Arc(ld, pm1, pm2, rm, bm(i), &
                                                     pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                     -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                     dbm, .FALSE., .FALSE., .FALSE.)) * of0
                            end if
                        else
                            ! Case F
                            ! planet fully overlaps star, moon partially overlaps star, both overlap each other 
                            call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                            call phis(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                      pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                            call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                            lc(:, i) = (2 * Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                        - Arc(ld, -km, pm2, rm, bm(i), &
                                              0.d0, -km_rm, -km_bp, -km_bpm, -km_theta, &
                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                              dbm, .FALSE., .TRUE., .FALSE.) &
                                        - Arc(ld, pm1, km, rm, bm(i), &
                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                              0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                              dbm, .FALSE., .FALSE., .TRUE.) &
                                        - Arc(ld, pp1, pp2, rp, bp(i), &
                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                              dbm0, .TRUE., .FALSE., .FALSE.)) * of0
                        end if
                    else
                        if (bm(i) + rm .le. 1.d0) then 
                            if (bpm(i) + rm .le. rp) then
                                ! planet partially overlaps star, moon fully overlaps star but is completely overlapped by planet 
                                call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                    dbm0, .TRUE., .TRUE.)) * of0
                            else
                                ! planet partially overlaps star, moon fully overlaps star and only partially overlaps planet
                                call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                call phis(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                      pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                lc(:, i) = (2 * Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                            - Arc(ld, -kp, pp2, rp, bp(i), &
                                                  -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                  -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                  dbm0, .TRUE., .TRUE., .FALSE.) &
                                            - Arc(ld, pp1, kp, rp, bp(i), &
                                                  pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                  kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                  dbm0, .TRUE., .FALSE., .TRUE.) &
                                            - Arc(ld, pm1, pm2, rm, bm(i), &
                                                  pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                  -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                  dbm, .FALSE., .FALSE., .FALSE.)) * of0
                            end if
                        else
                            if (bpm(i) + rm .le. rp) then
                                ! planet and moon both partially overlap star but moon is fully overlapped by the planet
                                call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                    dbm0, .TRUE., .TRUE.)) * of0
                            else
                                ! bookmark
                                !call compute_theta(rp,  bp(i), theta, phip, theta_bp, theta_rp, phip_bp, phip_rp)
                                !call compute_theta(rm,  bm(i), theta, phim, theta_bm, theta_rm, phim_bm, phim_rm)
                                call kappas_p(rp, bp(i), kp, kps, kp_rp, kp_bp, kps_rp, kps_bp)
                                call kappas_m(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms, &
                                      km_rm, km_bp, km_bpm, km_theta, &
                                      kms_rm, kms_bp, kms_bpm, kms_theta)
                                
                                a = bm(i)
                                b = bp(i)
                                c = bpm(i)
                                if (b .gt. a) then
                                    tmp = b
                                    b = a
                                    a = tmp
                                end if
                                if (c .gt. b) then
                                    tmp = c
                                    c = b
                                    b = tmp
                                end if
                                if (b .gt. a) then
                                    tmp = b
                                    b = a
                                    a = tmp
                                end if
                                delta = (a + (b + c)) * (c - (a - b)) * (c + (a - b)) * (a + (b - c))
                                if (delta .lt. 0.d0) then
                                    delta = 0.d0
                                else
                                    delta = Sqrt(delta)
                                end if
                                phi = Atan2(delta, (bm(i) - bpm(i)) * (bm(i) + bpm(i)) + bp(i) * bp(i))   
                                
                                call phis(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2, pp_rp, pp_rm, pp_bpm, &
                                          pm_rp, pm_rm, pm_bpm, thetam_bp, thetam_bpm, thetam_theta)
                                
                                if (phi + kms .le. kps) then  
                                        if (abs(pp2) .gt. kp) then
                                            ! planet and moon both partially overlap the star and each other but the 
                                            ! moon-star overlap is contained within the planet-star overlap
                                            lc(:, i) = 2 * (Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                            - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                                dbm0, .TRUE., .TRUE.)) * of0
                                        else
                                            ! planet and moon both partially overlap star and each other but the 
                                            ! moon-star intersections are overlapped by the planet
                                            call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                            lc(:, i) = (2 * Fstar(ld, pi - kps, -kps_rp, 0.d0, -kps_bp, 0.d0, 0.d0) &
                                                        - Arc(ld, -kp, pp2, rp, bp(i), &
                                                              -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              dbm0, .TRUE., .TRUE., .FALSE.) &
                                                        - Arc(ld, pp1, kp, rp, bp(i), &
                                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                              kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                              dbm0, .TRUE., .FALSE., .TRUE.) &
                                                        - Arc(ld, pm1, pm2, rm, bm(i), &
                                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                              dbm, .FALSE., .FALSE., .FALSE.)) * of0
                                        end if
    
                                else if (phi + kps .le. kms) then
                                    call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                    if ((bp(i) - rp) .le. (bm(i) - rm)) then
                                        ! Case L
                                        ! planet and moon both partially overlap the star and each other but the 
                                        ! planet-star intersections are overlapped by the moon
                                        lc(:, i) = (2 * Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                                        - Arc(ld, -km, pm2, rm, bm(i), &
                                                              0.d0, -km_rm, -km_bp, -km_bpm, -km_theta, &
                                                              -pm_rp, -pm_rm, thetam_bp, -pm_bpm + thetam_bpm, thetam_theta, &
                                                              dbm, .FALSE., .TRUE., .FALSE.) &
                                                        - Arc(ld, pm1, km, rm, bm(i), &
                                                              pm_rp, pm_rm, thetam_bp, pm_bpm + thetam_bpm, thetam_theta, &
                                                              0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                              dbm, .FALSE., .FALSE., .TRUE.) &
                                                        - Arc(ld, pp1, pp2, rp, bp(i), &
                                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              dbm0, .TRUE., .FALSE., .FALSE.)) * of0
                                    else
                                        ! planet and moon both partially overlap the star and each other but 
                                        ! the planet-star overlap is  entirely within the moon-star overlap
                                        lc(:, i) = 2 * (Fstar(ld, pi - kms, 0.d0, -kms_rm, -kms_bp, -kms_bpm, -kms_theta) &
                                                        - F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                            dbm, .FALSE., .TRUE.)) * of0
                                    end if
                                else
                                    ! bookmark
                                    d1 = rm * rm + bm(i) * bm(i) - 2 * rm * bm(i) * Cos(pm2)
                                    d2 = rm * rm + bm(i) * bm(i) - 2 * rm * bm(i) * Cos(pm1)
                                    call bm_x(bp(i), bm(i), bpm(i), cth(i), sth(i), dbm)
                                    if ((d1 .gt. 1.d0) .AND. (d2 .gt. 1.d0)) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! but the planet/moon overlap does not overlap the star
                                        if (pp1 * pp2 .lt. 0.d0) then
                                            lc(:, i) = 0.d0
                                        else
                                            lc(:, i) = 2 * (Fstar(ld, pi - (kps + kms), -kps_rp, -kms_rm, &
                                                              -(kps_bp + kms_bp), -kms_bpm, -kms_theta) &
                                                        - F(ld, kp, rp, bp(i), kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                            dbm0, .TRUE., .TRUE.) &
                                                        - F(ld, km, rm, bm(i), 0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                            dbm, .FALSE., .TRUE.)) * of0
                                        end if
                                    else if ((d1 .le. 1.d0) .AND. (d2 .le. 1.d0)) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap fully overlapping the star
                                        lc(:, i) = (2 * Fstar(ld, pi - (kps + kms), -kps_rp, -kms_rm, &
                                                              -(kps_bp + kms_bp), -kms_bpm, -kms_theta) &
                                                        - Arc(ld, -km, -pm1, rm, bm(i), &
                                                              0.d0, -km_rm, -km_bp, -km_bpm, -km_theta, &
                                                              -pm_rp, -pm_rm, -thetam_bp, -pm_bpm - thetam_bpm, -thetam_theta, &
                                                              dbm, .FALSE., .TRUE., .FALSE.) &
                                                        - Arc(ld, -pm2, km, rm, bm(i), &
                                                              pm_rp, pm_rm, -thetam_bp, pm_bpm - thetam_bpm, -thetam_theta, &
                                                              0.d0, km_rm, km_bp, km_bpm, km_theta, &
                                                              dbm, .FALSE., .FALSE., .TRUE.) &
                                                        - Arc(ld, pp1, kp, rp, bp(i), &
                                                              pp_rp, pp_rm, 0.d0, pp_bpm, 1.d0, &
                                                              kp_rp, 0.d0, kp_bp, 0.d0, 0.d0, &
                                                              dbm0, .TRUE., .FALSE., .TRUE.) &
                                                        - Arc(ld, -kp, pp2, rp, bp(i), &
                                                              -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              dbm0, .TRUE., .TRUE., .FALSE.)) * of0
                                    else
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap partially overlapping the star
                                        
                                        ! there might be a mistake somewhere in here... 
                                        phi_bpm = bp(i) * sth(i) / (bm(i) * bm(i))
                                        phi_theta = bpm(i) * (bp(i) * cth(i) - bpm(i)) / (bm(i) * bm(i))
                                        phi_bp = - bpm(i) * sth(i) / (bm(i) * bm(i))
                                        
                                        lc(:, i) = (2 * Fstar(ld, pi - 0.5 * (kps + kms + phi), &
                                                              -0.5 * kps_rp, -0.5 * kms_rm, &
                                                              -0.5 * (kps_bp + kms_bp + phi_bp), -0.5 * (kms_bpm + phi_bpm), &
                                                              -0.5 * (kms_theta + phi_theta)) &
                                                        - Arc(ld, -pm2, km, rm, bm(i), &
                                                              pm_rp, pm_rm, -thetam_bp, pm_bpm - thetam_bpm, -thetam_theta, &
                                                              0.d0, km_rm, km_bp, km_bpm, km_theta,  &
                                                              dbm, .FALSE., .FALSE., .TRUE.) &
                                                        - Arc(ld, -kp, pp2, rp, bp(i), &
                                                              -kp_rp, 0.d0, -kp_bp, 0.d0, 0.d0, &
                                                              -pp_rp, -pp_rm, 0.d0, -pp_bpm, 1.d0, &
                                                              dbm0, .TRUE., .TRUE., .FALSE.)) * of0
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if  
        end if
        lc(7, i) = lc(7, i) + (lc(1, i) - 1.d0) * pithird * of0
        lc(8, i) = lc(8, i) + (lc(1, i) - 1.d0) * 0.5d0 * pithird * of0
        lc(:, i) = lc(:, i) - f0 * of0
    end do
    return
    
end

! work out the right sign and order of the integration and call the integration routine 
! to integrate along an arbitrary arc of the planet or moon 
function Arc(ld, phi1, phi2, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta, &
            phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta, dbm, pflag, limbflag1, limbflag2)
                    
    real*8, dimension(8) :: Arc

    logical :: pflag, limbflag1, limbflag2
    real*8 :: phi1, phi2, r, b
    real*8, dimension(3) :: ld
    real*8, dimension(3) :: dbm
    real*8 :: phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta
    real*8 :: phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta
        
    if (phi1 < 0) then
        if (phi2 > 0) then
            Arc = F(ld, phi2, r, b, phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta, &
                    dbm, pflag, limbflag2) &
                + F(ld, -phi1, r, b, -phi1_rp, -phi1_rm, -phi1_bp, -phi1_bpm, -phi1_theta, &
                    dbm, pflag, limbflag1)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * Fcomplete(ld, r, b, dbm, pflag) &
                    + F(ld, -phi1, r, b, -phi1_rp, -phi1_rm, -phi1_bp, -phi1_bpm, -phi1_theta, &
                        dbm, pflag, limbflag1) &
                    - F(ld, -phi2, r, b, -phi2_rp, -phi2_rm, -phi2_bp, -phi2_bpm, -phi2_theta, &
                        dbm, pflag, limbflag2)
                return
            else
                Arc = - F(ld, -phi2, r, b, -phi2_rp, -phi2_rm, -phi2_bp, -phi2_bpm, -phi2_theta, &
                          dbm, pflag, limbflag2) &
                      + F(ld, -phi1, r, b, -phi1_rp, -phi1_rm, -phi1_bp, -phi1_bpm, -phi1_theta, &
                          dbm, pflag, limbflag1)
                return
            end if
        end if
    else
        if (phi2 < 0) then
            Arc = 2 * Fcomplete(ld, r, b, dbm, pflag) &
                - F(ld, phi1, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta, &
                    dbm, pflag, limbflag1) &
                - F(ld, -phi2, r, b, -phi2_rp, -phi2_rm, -phi2_bp, -phi2_bpm, -phi2_theta, &
                    dbm, pflag, limbflag2)
            return
        else
            if (phi2 < phi1) then
                Arc = 2 * Fcomplete(ld, r, b, dbm, pflag) &
                    + F(ld, phi2, r, b, phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta, &
                        dbm, pflag, limbflag2) &
                    - F(ld, phi1, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta, &
                        dbm, pflag, limbflag1)
                return
            else
                Arc = F(ld, phi2, r, b, phi2_rp, phi2_rm, phi2_bp, phi2_bpm, phi2_theta, &
                        dbm, pflag, limbflag2) &
                    - F(ld, phi1, r, b, phi1_rp, phi1_rm, phi1_bp, phi1_bpm, phi1_theta, &
                        dbm, pflag, limbflag1)
                return
            end if
        end if
    end if
    
    return
    
end function

! integrate along the limb of the star
function Fstar(ld, phi, phi_rp, phi_rm, phi_bp, phi_bpm, phi_theta)

    real*8, dimension(8) :: Fstar
    real*8, dimension(3) :: F_, F_rp, F_rm, F_bp, F_bpm, F_theta

    real*8 :: phi
    real*8 :: Fc_phi, Fq_phi, Fl_phi
    real*8 :: phi_bp, phi_rp, phi_rm, phi_bpm, phi_theta
    real*8, dimension(3) :: ld
        
    F_(1) = 0.5 * phi
    Fc_phi = 0.5
    F_rp(1) = Fc_phi * phi_rp
    F_rm(1) = Fc_phi * phi_rm
    F_bp(1) = Fc_phi * phi_bp
    F_bpm(1) = Fc_phi * phi_bpm
    F_theta(1) = Fc_phi * phi_theta
    
    F_(2) = phi * o3
    Fl_phi = o3
    F_rp(2) = Fl_phi * phi_rp
    F_rm(2) = Fl_phi * phi_rm
    F_bp(2) = Fl_phi * phi_bp
    F_bpm(2) = Fl_phi * phi_bpm
    F_theta(2) = Fl_phi * phi_theta

    F_(3) = 0.25 * phi
    Fq_phi = 0.25d0
    F_rp(3) = Fq_phi * phi_rp
    F_rm(3) = Fq_phi * phi_rm
    F_bp(3) = Fq_phi * phi_bp
    F_bpm(3) = Fq_phi * phi_bpm
    F_theta(3) = Fq_phi * phi_theta
    
    Fstar(1) = Sum(ld * F_)
    Fstar(2) = Sum(ld * F_rp)
    Fstar(3) = Sum(ld * F_rm)
    Fstar(4) = Sum(ld * F_bp)
    Fstar(5) = Sum(ld * F_bpm)
    Fstar(6) = Sum(ld * F_theta)
    Fstar(7) = - F_(1) + F_(2)
    Fstar(8) = -2 * F_(1) + 2 * F_(2) + F_(3)
    return
    
end function

! integrate around the entire planet/moon 
function Fcomplete(ld, r, b, dbm, pflag)

    real*8, dimension(8) :: Fcomplete
    real*8, dimension(3) :: F_, F_r, F_b
    real*8 :: sumF_b
    
    ! Are we integrating along the edge of the moon or the planet? 
    logical :: pflag
    
    ! Limb darkening params
    real*8, dimension(3) :: ld
    
    ! self explanatory
    real*8 :: r, b
    
    ! derivatives of input parameters 
    real*8, dimension(3) :: dbm
    
    ! convenient parameters
    real*8 :: r2, b2, br, bmr, bpr, obmr
    real*8 :: x, ox, o
    
    ! For the integral
    real*8 :: alpha, beta, gamma, n, m
    real*8 :: sgn
    real*8 :: ur, ub, vb, sqomm
    
    ! Elliptic integrals
    real*8 :: ellippi, ellipe, ellipf
    real*8 :: eplusf, eplusf_r, eplusf_b
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    obmr = 1.d0 / bmr
    bpr = b + r
    br = b * r
    
    F_(1) = r2 * pihalf  
    F_r(1) = r * pi
    F_b(1) = 0.d0
    
    F_(3) = pihalf * r2 * (b2 + 0.5 * r2) 
    F_r(3) = pi * r * (b2 + r2)
    F_b(3) = pi * r2 * b

    x = Sqrt((1.d0 - bmr) * (1.d0 + bmr))
    ox = 1.d0 / x
            
    alpha = (7 * r2 + b2 - 4.d0) * o9 * x
    beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) * o9 * ox        
            
    n = -4 * br * obmr * obmr
    m = 4 * br * ox * ox
            
    ur = 2 * r * x
    ub = x * ((b + 1.d0) * (b - 1.d0) + r2) / (3 * b)
    vb = ((-1.d0 + b2)**2.d0 - 2 * (1.d0 + b2) * r2 + r2 * r2) * o3 * ox / b
            
    sqomm = Sqrt(1.d0 - m)
    o = 1.d0
            
    if (b .eq. r) then 
        ellippi = 0.d0
        sgn = 0.d0
        gamma = 0.d0
    else
        ellippi = cel((sqomm), (1.d0 - n), (o), (o))
        sgn = Sign(1.d0, bmr)
        gamma = bpr * ox * o3 * obmr
    end if 
            
    ellipe = cel((sqomm), (o), (o), 1.d0 - m)
    ellipf = cel((sqomm), (o), (o), (o))
             
    eplusf = alpha * ellipe + beta * ellipf
    eplusf_r = ur * ellipe
    eplusf_b = ub * ellipe + vb * ellipf
            
    F_(2) = eplusf + gamma * ellippi + pisixth * (1.d0 - sgn)
    F_r(2) = eplusf_r
    F_b(2) = eplusf_b

    Fcomplete = 0.d0
    sumF_b = Sum(ld * F_b)
    
    if (pflag) then
        Fcomplete(1) = Sum(ld * F_)
        Fcomplete(2) = Sum(ld * F_r)
        Fcomplete(4) = Sum(ld * F_b)
        Fcomplete(7) = - F_(1) + F_(2)
        Fcomplete(8) = -2 * F_(1) + 2 * F_(2) + F_(3)
    else
        Fcomplete(1) = Sum(ld * F_)
        Fcomplete(3) = Sum(ld * F_r)
        Fcomplete(4) = sumF_b * dbm(1)
        Fcomplete(5) = sumF_b * dbm(2)
        Fcomplete(6) = sumF_b * dbm(3)
        Fcomplete(7) = - F_(1) + F_(2)
        Fcomplete(8) = -2 * F_(1) + 2 * F_(2) + F_(3)
    end if
    return

end function

! evaluate the integral at one arbitrary limit along the planet or moon's boundary 
function F(ld, phi, r, b, phi_rp, phi_rm, phi_bp, phi_bpm, phi_theta, dbm, pflag, limbflag)

    real*8, dimension(8) :: F
    real*8, dimension(3) :: F_, F_rp, F_rm, F_bp, F_bpm, F_theta
    
    ! Are we integrating along the edge of the moon or the planet? 
    logical :: pflag, limbflag
    
    ! Limb darkening params
    real*8, dimension(3) :: ld
    
    ! self explanatory
    real*8 :: phi, r, b
    
    ! derivatives of input parameters 
    real*8 :: phi_bp, phi_rp, phi_rm, phi_bpm, phi_theta
    real*8, dimension(3) :: dbm
    
    ! convenient parameters
    real*8 :: sphi, cphi, tphihalf, sphihalf, cphihalf
    real*8 :: r2, b2, br, bmr, bpr, obmr
    real*8 :: x, y, z, ox, oy, oz, tans, o
    
    ! components of flux and their derivatives
    real*8 :: Fc_phi, Fc_r, Fc_b
    real*8 :: Fq_phi, Fq_r, Fq_b
    real*8 :: Fl_phi, Fl_r, Fl_b
    
    ! For the integral
    real*8 :: alpha, beta, gamma, d, n, m
    real*8 :: ur, vr, ub, vb, pr, pb, sqomm
    
    ! Elliptic integrals
    real*8 :: ellippi, ellipe, ellipf
    real*8 :: eplusf, eplusf_r, eplusf_b
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    obmr = 1.d0 / bmr
    bpr = b + r
    br = b * r
    
    cphi = Cos(phi)
    sphi = Sin(phi)
    sphihalf = Sin(phi * 0.5)
    cphihalf = Cos(phi * 0.5)
    tphihalf = sphihalf / cphihalf
    
    if (phi .eq. pi) then
        F = Fcomplete(ld, r, b, dbm, pflag)
        return
    end if
        
    F_(1) = 0.5 * (r2 * phi - br * sphi)
    Fc_phi = 0.5 * (r2 - br * cphi)
    Fc_b = -0.5 * r * sphi
    Fc_r = 0.5 * (2 * r * phi - b * sphi)
    
    F_rp(1) = Fc_phi * phi_rp
    F_rm(1) = Fc_phi * phi_rm
    F_bp(1) = Fc_phi * phi_bp
    F_bpm(1) = Fc_phi * phi_bpm
    F_theta(1) = Fc_phi * phi_theta
    
    F_(3) = 0.25 * (r * (r * (2 * b2 + r2) * phi &
          + b * (-b2 - 3 * r2 + br * cphi) * sphi))
    Fq_phi = 0.25 * (r * (2 * b2 * r + r2 * r + b * (-((b2 + 3 * r2) * cphi) &
           + br * (1.d0 - 2 * sphi * sphi))))
    Fq_r = r * (b2 + r2) * phi - 0.25 * (b * (b2 + 9 * r2 - 2 * br * cphi) * sphi)
    Fq_b = 0.25 * (r * (-3 * (b2 + r2) * sphi + br * (4 * phi + 2 * sphi * cphi)))
                
    F_rp(3) = Fq_phi * phi_rp
    F_rm(3) = Fq_phi * phi_rm
    F_bp(3) = Fq_phi * phi_bp
    F_bpm(3) = Fq_phi * phi_bpm
    F_theta(3) = Fq_phi * phi_theta
    
    y = Sqrt(br)
    oy = 1.d0 / y
    x = b2 + r2 - 2 * br * cphi
    ox = 1.d0 / x

    pr = - b * sphi * o3 * ox
    pb = r * sphi * o3 * ox
    
    if (bpr .gt. 1.d0) then
            
        alpha = 2 * y * (7 * r2 + b2 - 4.d0) * o9
        beta = -(3.d0 + 2*r*(b2 * b + 5 * b2 * r + 3*r*(-2.d0 + r2) + b*(-4.d0 + 7*r2))) * o9 * 0.5 * oy
        gamma = bpr * o3 * 0.5 * oy * obmr
            
        m = (1.d0 - bmr) * (1.d0 + bmr) * 0.25 * oy * oy
        n = ((bmr + 1.d0) * (bmr - 1.d0)) * obmr * obmr
        sqomm = Sqrt(1.d0 - m)

        ur = 4 * r * y
        vr = - r * (bpr + 1.d0) * (bpr - 1.d0) * oy
        ub = 2 * r * (b2 + r2 - 1.d0) * o3 * oy
        vb = vr * o3
            
        if (limbflag) then
            
            d = phi * o3 * 0.5 - Atan(bpr * tphihalf * obmr) * o3
            Fl_phi = (r2 - br * cphi) * o3 * ox
                
            o = 1.d0
            ellippi = cel((sqomm), 1.d0 - n, (o), (o))
            ellipe = cel((sqomm), (o), (o), (1.d0 - m))
            ellipf = cel((sqomm), (o), (o), (o))
                
            eplusf = alpha * ellipe + beta * ellipf
            eplusf_r = ur * ellipe + vr * ellipf
            eplusf_b = ub * ellipe + vb * ellipf
                
        else   
            
            z = Sqrt(1.d0 - x)
            oz = 1.d0 / z
            d = o3 * (phi * 0.5 - Atan(bpr * tphihalf * obmr)) &
                - (2 * br * o9) * sphi * z

            Fl_phi = o3 * (1.d0 - z ** 3.d0) * (r2 - br * cphi) * ox

            pr = - ((1.d0 - x)**(1.5d0) - 1.d0) * pr
            pb = (1.d0 - (1.d0 + x) * Sqrt(1.d0 - x)) * pb
                                        
            tans = 1.d0 / Sqrt(m / (sphihalf * sphihalf) - 1.d0)
            o = 1.d0
            ellippi = el3((tans), (sqomm), 1.d0 - n)
            ellipe = el2((tans), (sqomm), (o), (1.d0 - m))
            ellipf = el2((tans), (sqomm), (o), (o))
                
            eplusf = alpha * ellipe + beta * ellipf
            eplusf_r = ur * ellipe + vr * ellipf
            eplusf_b = ub * ellipe + vb * ellipf

        end if
        F_(2) = eplusf + gamma * ellippi + d
        Fl_r = eplusf_r + pr
        Fl_b = eplusf_b + pb
            
    else
        
        y = Sqrt((1.d0 - bmr) * (1.d0 + bmr))
        oy = 1.d0 / y
            
        alpha = (7 * r2 + b2 - 4.d0) * y * o9
        beta = (r2 * r2 + b2 * b2 + r2 - b2 * (5.d0 + 2 * r2) + 1.d0) * o9 * oy
        gamma = bpr * o3 * oy * obmr
                        
        n = -4 * br * obmr * obmr
        m = 4 * br * oy * oy
            
        sqomm = Sqrt(1.d0 - m)
            
        pr = - ((1.d0 - x)**(1.5d0) - 1.d0) * pr
        pb = (1.d0 - (1.d0 + x) * Sqrt(1.d0 - x)) * pb
            
        ur = 2 * r * y
        ub = y * ((b + 1.d0) * (b - 1.d0) + r2) / (3 * b)
        vb = (b2 * b2 + ((r - 1.d0) * (r + 1.d0))**2.d0 - 2 * b2 * (1.d0 + r2)) / (3 * b * y) 
            
        d = o3 * (phi * 0.5 - Atan(bpr * tphihalf * obmr)) &
            - 2 * br * o9 * sphi * Sqrt(1.d0 - x)
             
        o = 1.d0
        ellipe = el2((tphihalf), (sqomm), (o), (1.d0 - m))
        ellipf = el2((tphihalf), (sqomm), (o), (o))
            
        if (b .eq. r) then
            Fl_phi = o3 * (1.d0 - (1.d0 - x) ** 1.5d0) * (r2 - br * cphi) * ox
            ellippi = 0.d0
            gamma = 0.d0
        else
            Fl_phi = o3 * (1.d0 - (1.d0 - x) ** 1.5d0) * (r2 - br * cphi) * ox
            ellippi = el3((tphihalf), (sqomm), (1.d0 - n))
        end if
                
        eplusf = alpha * ellipe + beta * ellipf
        eplusf_r = ur * ellipe
        eplusf_b = ub * ellipe + vb * ellipf 
            
        F_(2) = eplusf + gamma * ellippi + d
        Fl_r = eplusf_r + pr
        Fl_b = eplusf_b + pb
    end if
                
    F_rp(2) = Fl_phi * phi_rp
    F_rm(2) = Fl_phi * phi_rm
    F_bp(2) = Fl_phi * phi_bp
    F_bpm(2) = Fl_phi * phi_bpm
    F_theta(2) = Fl_phi * phi_theta
    
    if (pflag) then
        F_bp(3) = F_bp(3) + Fq_b
        F_rp(3) = F_rp(3) + Fq_r
        
        F_bp(1) = F_bp(1) + Fc_b
        F_rp(1) = F_rp(1) + Fc_r
        
        F_bp(2) = F_bp(2) + Fl_b
        F_rp(2) = F_rp(2) + Fl_r
    else
        F_theta(3) = F_theta(3) + Fq_b * dbm(3)
        F_bpm(3) = F_bpm(3) + Fq_b * dbm(2)
        F_bp(3) = F_bp(3) + Fq_b * dbm(1)
        F_rm(3) = F_rm(3) + Fq_r
        
        F_theta(1) = F_theta(1) + Fc_b * dbm(3)
        F_bpm(1) = F_bpm(1) + Fc_b * dbm(2)
        F_bp(1) = F_bp(1) + Fc_b * dbm(1)
        F_rm(1) = F_rm(1) + Fc_r
        
        F_theta(2) = F_theta(2) + Fl_b * dbm(3)
        F_bpm(2) = F_bpm(2) + Fl_b * dbm(2)
        F_bp(2) = F_bp(2) + Fl_b * dbm(1)
        F_rm(2) = F_rm(2) + Fl_r
    end if
    
    F(1) = Sum(ld * F_)
    F(2) = Sum(ld * F_rp)
    F(3) = Sum(ld * F_rm)
    F(4) = Sum(ld * F_bp)
    F(5) = Sum(ld * F_bpm)
    F(6) = Sum(ld * F_theta)
    F(7) = - F_(1) + F_(2)
    F(8) = -2 * F_(1) + 2 * F_(2) + F_(3)

    return
end function

end module phot