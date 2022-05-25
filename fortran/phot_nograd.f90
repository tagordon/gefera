module phot_nograd
use iso_c_binding
use ellip

implicit none

!real*8, parameter :: pi = 4.d0 * Atan(1.d0), pihalf = 2.d0 * Atan(1.d0)
real*8, parameter :: twopithree = 8.d0 * Atan(1.d0) / 3.d0, twopi = 8.d0 * Atan(1.d0)
real*8, parameter :: o3 = 1.d0 / 3.d0, o9 = 1.d0 / 9.d0
real*8, parameter :: pithird = 4.d0 * Atan(1.d0) / 3.d0, pisixth = 4.d0 * Atan(1.d0) / 6.d0

contains

subroutine phis_ng(rp, rm, bp, bm, bpm, cth, sth, pp1, pp2, pm1, pm2)

    real*8 :: rp, rm, bp, bm, bpm, cth, sth
    
    ! intersection angle from planet center, intersection from moon center, same 
    ! values relative to bp and bm vectors respectively
    real*8 :: pp, pm, pp1, pm1, pp2, pm2
    
    ! Four times the area of the triangle formed by rm, rp, and bpm
    real*8 :: delta
    ! Variables used in sorting the sides of the triangle
    real*8 :: a, b, c, tmp
    
    ! angle between bpm vector and bp vector, 
    ! angle between bpm vector and bm vector
    real*8 :: theta, thetam
    
    ! for avoiding divisions
    real*8 :: obm
    
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
    pm = Atan2(delta, (rm - rp) * (rm + rp) + bpm * bpm)   
    pp = Atan2(delta, (rp - rm) * (rp + rm) + bpm * bpm)

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

subroutine kappas_p_ng(rp, bp, kp, kps)

    ! kp = angle to intersection from center of planet, 
    ! kps = angle to intersection from center of star 
    real*8 :: rp, bp, kp, kps
    
    ! variables used in sorting sides of triangle
    real*8 :: a, b, c
    
    ! four times the area of the triangle with sides rp, bp, and 1
    real*8 :: delta
    
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
    kps = Atan2(delta, (1.d0 - rp) * (1.d0 + rp) + bp * bp)
    kp = Atan2(delta, (rp - 1.d0) * (rp + 1.d0) + bp * bp)
end 

subroutine kappas_m_ng(rm, bp, bm, bpm, cth, sth, km, kms)

    ! km = angle to interection from center of moon, 
    ! kms = angle to intersection from center of planet
    real*8 :: rm, bp, bm, bpm, cth, sth, km, kms
    
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
end 


! main loop to compute the flux at each timestep by finding the correct geometry and
! calling the integration routines 
subroutine flux_ng(c1, c2, rp, rm, bp, bpm, cth, sth, lc, j) bind(C, name="flux_ng")

    integer (c_int), bind(C) :: j
    integer :: i
    real (c_double), bind(C) :: rp, rm
    real (c_double), bind(C), dimension(j) :: bp, cth, sth, bpm
    real (c_double), bind(C), intent(out), dimension(j) :: lc
    real*8 :: f0
    real*8, dimension(3) :: ld
    real (c_double), bind(C) :: c1, c2
    real*8 :: of0
    
    ! half angles: from planet to planet/star intersection, moon to moon/star
    ! intersection, spanned by planet on limb of star, spanned by moon on 
    ! limb of star
    real*8 :: kp, km, kps, kms
    
    ! angles to planet-moon intersection from planet center 
    !relative to bp vector (and derivatives)
    real*8 :: pp1, pp2
    
    ! angles to planet-moon intersection from moon center 
    ! relative to bm vector (and derivatves)
    real*8 :: pm1, pm2
    
    ! used to determine cases for three body overlaps, might not be needed. 
    ! Check if some of these (costheta, cosphi) can be removed when optimizing things later 
    real*8 :: phi, d1, d2, delta, a, b, c, tmp
    
    ! For chain rule stuff
    real*8, dimension(j) :: bm
    
    bm = Sqrt((bp - bpm)**2.d0 + 2 * bp * bpm * (1.d0 - cth))
    
    ld(1) = 1.d0 - c1 - 2 * c2
    ld(2) = c1 + 2 * c2
    ld(3) = c2
    
    ! normalization factors 
    f0 = ld(1) * pi + ld(2) * twopithree + ld(3) * pihalf
    
    of0 = 1.d0 / f0
    
    do i=1,j,1
    
        lc(i) = f0 * of0
        
        if (rp .gt. 1.d0 + bp(i)) then
            lc(i) = 0.d0
        else if (rm .gt. 1.d0 + bm(i)) then 
            lc(i) = 0.d0
        else if (bpm(i) .gt. rp + rm) then
            ! moon and planet don't overlap each other 
            if (bp(i) .lt. 1.d0 - rp) then
                ! planet completely inside star
                lc(i) = lc(i) - 2 * Fcomplete_ng(ld, rp, bp(i)) * of0
            else if (bp(i) .lt. 1.d0 + rp) then
                ! planet partially overlaps star
                call kappas_p_ng(rp, bp(i), kp, kps)
                lc(i) = lc(i) - 2 * (Fstar_ng(ld, kps) &
                                    + F_ng(ld, kp, rp, bp(i), .TRUE.)) * of0
            end if
            if (bm(i) .lt. 1.d0 - rm) then
                ! moon completely inside star
                lc(i) = lc(i) - 2 * Fcomplete_ng(ld, rm, bm(i)) * of0
            else if (bm(i) .lt. 1.d0 + rm) then
                ! moon partially overlaps star
                call kappas_m_ng(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms)
                lc(i) = lc(i) - 2 * (Fstar_ng(ld, kms) &
                                    + F_ng(ld, km, rm, bm(i), .TRUE.)) * of0
            end if
        else
            ! moon and planet do overlap each other 
            if (bp(i) .gt. rp + 1.d0) then
                if (bm(i) .gt. rm + 1.d0) then
                    ! neither moon nor planet overlap star
                    lc(i) = f0 * of0
                else
                    ! moon partially overlaps star, planet does not overlap star
                    call kappas_m_ng(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms)
                    lc(i) = 2 * (Fstar_ng(ld, pi - kms) &
                                    - F_ng(ld, km, rm, bm(i), .TRUE.)) * of0

                end if
            else
                if (bm(i) .gt. rm + 1.d0) then
                    ! planet partially overlaps star, moon is outside of star
                    call kappas_p_ng(rp, bp(i), kp, kps)
                    lc(i) = 2 * (Fstar_ng(ld, pi - kps) &
                                    - F_ng(ld, kp, rp, bp(i), .TRUE.)) * of0
                else
                    if (bp(i) + rp .le. 1.d0) then
                        if (bm(i) + rm .le. 1.d0) then
                            if (bpm(i) + rm .le. rp) then
                                ! moon and planet both overlap star, moon is completely overlapped by planet
                                lc(i) = (f0 - 2 * Fcomplete_ng(ld, rp, bp(i))) * of0
                            else
                                ! Case E
                                ! moon and planet both overlap star, moon and planet partially overlap each other 
                                call phis_ng(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2)
                                lc(i) = (f0 - Arc_ng(ld, pp1, pp2, rp, bp(i), .FALSE., .FALSE.) &
                                                - Arc_ng(ld, pm1, pm2, rm, bm(i), .FALSE., .FALSE.)) * of0
                            end if
                        else
                            ! planet fully overlaps star, moon partially overlaps star, both overlap each other 
                            call phis_ng(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2)
                            call kappas_m_ng(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms)
                            lc(i) = (2 * Fstar_ng(ld, pi - kms) &
                                        - Arc_ng(ld, -km, pm2, rm, bm(i), .TRUE., .FALSE.) &
                                        - Arc_ng(ld, pm1, km, rm, bm(i), .FALSE., .TRUE.) &
                                        - Arc_ng(ld, pp1, pp2, rp, bp(i), .FALSE., .FALSE.)) * of0
                        end if
                    else
                        if (bm(i) + rm .le. 1.d0) then 
                            if (bpm(i) + rm .le. rp) then
                                ! planet partially overlaps star, moon fully overlaps star but is completely overlapped by planet 
                                call kappas_p_ng(rp, bp(i), kp, kps)
                                lc(i) = 2 * (Fstar_ng(ld, pi - kps) &
                                                - F_ng(ld, kp, rp, bp(i), .TRUE.)) * of0
                            else
                                ! planet partially overlaps star, moon fully overlaps star and only partially overlaps planet
                                call phis_ng(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2)
                                call kappas_p_ng(rp, bp(i), kp, kps)
                                lc(i) = (2 * Fstar_ng(ld, pi - kps) &
                                            - Arc_ng(ld, -kp, pp2, rp, bp(i), .TRUE., .FALSE.) &
                                            - Arc_ng(ld, pp1, kp, rp, bp(i), .FALSE., .TRUE.) &
                                            - Arc_ng(ld, pm1, pm2, rm, bm(i), .FALSE., .FALSE.)) * of0
                            end if
                        else
                            if (bpm(i) + rm .le. rp) then
                                ! planet and moon both partially overlap star but moon is fully overlapped by the planet
                                call kappas_p_ng(rp, bp(i), kp, kps)
                                lc(i) = 2 * (Fstar_ng(ld, pi - kps) &
                                                - F_ng(ld, kp, rp, bp(i), .TRUE.)) * of0
                            else
                                call kappas_p_ng(rp, bp(i), kp, kps)
                                call kappas_m_ng(rm, bp(i), bm(i), bpm(i), cth(i), sth(i), km, kms)
                                
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
                                
                                call phis_ng(rp, rm, bp(i), bm(i), bpm(i), cth(i), sth(i), pp1, pp2, pm1, pm2)
                                
                                if (phi + kms .le. kps) then 
                                        if (abs(pp2) .gt. kp) then
                                            ! planet and moon both partially overlap the star and each other but the 
                                            ! moon-star overlap is contained within the planet-star overlap
                                            lc(i) = 2 * (Fstar_ng(ld, pi - kps) &
                                                            - F_ng(ld, kp, rp, bp(i), .TRUE.)) * of0
                                        else
                                            ! planet and moon both partially overlap star and each other but the 
                                            ! moon-star intersections are overlapped by the planet
                                            lc(i) = (2 * Fstar_ng(ld, pi - kps) &
                                                        - Arc_ng(ld, -kp, pp2, rp, bp(i), .TRUE., .FALSE.) &
                                                        - Arc_ng(ld, pp1, kp, rp, bp(i), .FALSE., .TRUE.) &
                                                        - Arc_ng(ld, pm1, pm2, rm, bm(i), .FALSE., .FALSE.)) * of0
                                        end if
                                else if (phi + kps .le. kms) then
                                    if ((bp(i) - rp) .le. (bm(i) - rm)) then
                                        ! planet and moon both partially overlap the star and each other but the 
                                        ! planet-star intersections are overlapped by the moon
                                        lc(i) = (2 * Fstar_ng(ld, pi - kms) &
                                                        - Arc_ng(ld, -km, pm2, rm, bm(i), .TRUE., .FALSE.) &
                                                        - Arc_ng(ld, pm1, km, rm, bm(i), .FALSE., .TRUE.) &
                                                        - Arc_ng(ld, pp1, pp2, rp, bp(i), .FALSE., .FALSE.)) * of0
                                    else
                                        ! planet and moon both partially overlap the star and each other but 
                                        ! the planet-star overlap is  entirely within the moon-star overlap
                                        lc(i) = 2 * (Fstar_ng(ld, pi - kms) &
                                                        - F_ng(ld, km, rm, bm(i), .TRUE.)) * of0
                                    end if
                                else
                                    ! bookmark
                                    d1 = rm * rm + bm(i) * bm(i) - 2 * rm * bm(i) * Cos(pm2)
                                    d2 = rm * rm + bm(i) * bm(i) - 2 * rm * bm(i) * Cos(pm1)
                                    if ((d1 .gt. 1.d0) .AND. (d2 .gt. 1.d0)) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! but the planet/moon overlap does not overlap the star
                                        lc(i) = 2 * (Fstar_ng(ld, pi - (kps + kms)) &
                                                        - F_ng(ld, kp, rp, bp(i), .TRUE.) &
                                                        - F_ng(ld, km, rm, bm(i), .TRUE.)) * of0
                                    else if ((d1 .le. 1.d0) .AND. (d2 .le. 1.d0)) then
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap fully overlapping the star
                                        lc(i) = (2 * Fstar_ng(ld, pi - (kps + kms)) &
                                                        - Arc_ng(ld, -km, -pm1, rm, bm(i), .TRUE., .FALSE.) &
                                                        - Arc_ng(ld, -pm2, km, rm, bm(i), .FALSE., .TRUE.) &
                                                        - Arc_ng(ld, pp1, kp, rp, bp(i), .FALSE., .TRUE.) &
                                                        - Arc_ng(ld, -kp, pp2, rp, bp(i), .TRUE., .FALSE.)) * of0
                                    else
                                        ! planet and moon both partially overlap star and each other, 
                                        ! with the planet/moon overlap partially overlapping the star
                                        
                                        lc(i) = (2 * Fstar_ng(ld, pi - 0.5 * (kps + kms + phi)) &
                                                        - Arc_ng(ld, -pm2, km, rm, bm(i), .FALSE., .TRUE.) &
                                                        - Arc_ng(ld, -kp, pp2, rp, bp(i), .TRUE., .FALSE.)) * of0
                                    end if
                                end if
                            end if
                        end if
                    end if
                end if
            end if  
        end if
        lc(i) = lc(i) - f0 * of0
    end do
    return
    
end

! work out the right sign and order of the integration and call the integration routine 
! to integrate along an arbitrary arc of the planet or moon 
function Arc_ng(ld, phi1, phi2, r, b, limbflag1, limbflag2)
                    
    real*8 :: Arc_ng

    logical :: limbflag1, limbflag2
    real*8 :: phi1, phi2, r, b
    real*8, dimension(3) :: ld
        
    if (phi1 < 0) then
        if (phi2 > 0) then
            Arc_ng = F_ng(ld, phi2, r, b, limbflag2) &
                + F_ng(ld, -phi1, r, b, limbflag1)
            return
        else
            if (phi2 < phi1) then
                Arc_ng = 2 * Fcomplete_ng(ld, r, b) &
                    + F_ng(ld, -phi1, r, b, limbflag1) &
                    - F_ng(ld, -phi2, r, b, limbflag2)
                return
            else
                Arc_ng = - F_ng(ld, -phi2, r, b, limbflag2) &
                      + F_ng(ld, -phi1, r, b, limbflag1)
                return
            end if
        end if
    else
        if (phi2 < 0) then
            Arc_ng = 2 * Fcomplete_ng(ld, r, b) &
                - F_ng(ld, phi1, r, b, limbflag1) &
                - F_ng(ld, -phi2, r, b, limbflag2)
            return
        else
            if (phi2 < phi1) then
                Arc_ng = 2 * Fcomplete_ng(ld, r, b) &
                    + F_ng(ld, phi2, r, b, limbflag2) &
                    - F_ng(ld, phi1, r, b, limbflag1)
                return
            else
                Arc_ng = F_ng(ld, phi2, r, b, limbflag2) &
                    - F_ng(ld, phi1, r, b, limbflag1)
                return
            end if
        end if
    end if
    
    return
    
end function

! integrate along the limb of the star
function Fstar_ng(ld, phi)

    real*8 :: Fstar_ng
    real*8, dimension(3) :: F_

    real*8 :: phi
    real*8, dimension(3) :: ld
        
    F_(1) = 0.5 * phi
    F_(2) = phi * o3
    F_(3) = 0.25 * phi
    
    Fstar_ng = Sum(ld * F_)
    return
    
end function

! integrate around the entire planet/moon 
function Fcomplete_ng(ld, r, b)

    real*8 :: Fcomplete_ng
    real*8, dimension(3) :: F_
    
    ! Limb darkening params
    real*8, dimension(3) :: ld
    
    ! self explanatory
    real*8 :: r, b
    
    ! convenient parameters
    real*8 :: r2, b2, br, bmr, bpr, obmr
    real*8 :: x, ox, o

    
    ! For the integral
    real*8 :: alpha, beta, gamma, n, m
    real*8 :: sgn
    real*8 :: ur, ub, vb, sqomm
    
    ! Elliptic integrals
    real*8 :: ellippi
    real*8 :: eplusf
    
    r2 = r * r
    b2 = b * b
    bmr = b - r
    obmr = 1.d0 / bmr
    bpr = b + r
    br = b * r
    
    F_(1) = r2 * pihalf  
    F_(3) = pihalf * r2 * (b2 + 0.5 * r2) 
    
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

    eplusf = cel((sqomm), (o), alpha + beta, beta + alpha * (1.d0 - m))
            
    F_(2) = eplusf + gamma * ellippi + pisixth * (1.d0 - sgn)
    
    Fcomplete_ng = Sum(ld * F_)
    return

end function

! evaluate the integral at one arbitrary limit along the planet or moon's boundary 
function F_ng(ld, phi, r, b, limbflag)

    real*8 :: F_ng
    real*8, dimension(3) :: F_
    
    ! Are we integrating along the edge of the moon or the planet? 
    logical :: limbflag
    
    ! Limb darkening params
    real*8, dimension(3) :: ld
    
    ! self explanatory
    real*8 :: phi, r, b
    
    ! convenient parameters
    real*8 :: sphi, cphi, tphihalf, sphihalf, cphihalf
    real*8 :: r2, b2, br, bmr, bpr, obmr
    real*8 :: x, y, z, ox, oy, oz, tans, o
    
    ! For the integral
    real*8 :: alpha, beta, gamma, d, n, m
    real*8 :: ur, vr, ub, vb, pr, pb, sqomm
    
    ! Elliptic integrals
    real*8 :: ellippi
    real*8 :: eplusf
    
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
        F_ng = Fcomplete_ng(ld, r, b)
        return
    end if
        
    F_(1) = 0.5 * (r2 * phi - br * sphi)
    F_(3) = 0.25 * (r * (r * (2 * b2 + r2) * phi &
          + b * (-b2 - 3 * r2 + br * cphi) * sphi))
          
    y = Sqrt(br)
    oy = 1.d0 / y
    x = b2 + r2 - 2 * br * cphi
    ox = 1.d0 / x
          
    pr = - b * sphi * o3 * ox
    pb = r * sphi * o3 * ox
                
    if (bpr .gt. 1.d0) then
        
        y = Sqrt(br)
        oy = 1.d0 / y
        x = b2 + r2 - 2 * br * cphi
        ox = 1.d0 / x
            
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
                
            o = 1.d0
            ellippi = cel((sqomm), 1.d0 - n, (o), (o))
            eplusf = cel((sqomm), (o), alpha + beta, beta + alpha * (1.d0 - m))
        else                
            z = Sqrt((1.d0 - b) * (1.d0 + b) - r2 + 2 * br * cphi)
            oz = 1.d0 / z
            d = o3 * (phi * 0.5 - Atan(bpr * tphihalf * obmr)) &
                - (2 * br * o9) * sphi * z
                    
            pr = - ((1.d0 - x)**(1.5d0) - 1.d0) * pr
            pb = (1.d0 - (1.d0 + x) * Sqrt(1.d0 - x)) * pb 
             
            tans = 1.d0 / Sqrt(m / (sphihalf * sphihalf) - 1.d0)
                
            o = 1.d0
            ellippi = el3((tans), (sqomm), 1.d0 - n)
            eplusf = el2((tans), (sqomm), alpha + beta, beta + alpha * (1.d0 - m))

        end if
        F_(2) = eplusf + gamma * ellippi + d
            
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
            
        if (b .eq. r) then
            ellippi = 0.d0
            gamma = 0.d0
        else
            ellippi = el3((tphihalf), (sqomm), (1.d0 - n))
        end if
                
        eplusf = el2((tphihalf), (sqomm), alpha + beta, beta + alpha * (1.d0 - m))
            
        F_(2) = eplusf + gamma * ellippi + d
    end if
            
    F_ng = Sum(ld * F_)

    return
end function

end module phot_nograd