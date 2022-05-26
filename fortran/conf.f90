module conf
use iso_c_binding
use kepler
implicit none

contains

subroutine coords(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, j, &
    &xp, yp, zp, xm, ym, zm) bind(C, name="coords")
    
    integer (c_int), bind(C) :: j
    real*8 :: np, nm
    real*8 :: comegam, somegam
    real*8 :: cwp, swp, cip, cwm, swm, cim
    real*8, dimension(j) :: r, cosf, sinf, cosfw, sinfw
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip 
    real (c_double), bind(C) :: am, t0m, em, Pm, Om, wm, im
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), dimension(j), intent(out) :: xp, yp, zp, xm, ym, zm
            
    np = 2 * pi / Pp    
    nm = 2 * pi / Pm

    cwp = Cos(wp)
    swp = Sin(wp)
    cip = Cos(ip)
    
    comegam = Cos(Om)
    somegam = Sin(Om)
    cwm = Cos(wm)
    swm = Sin(wm)
    cim = Cos(im)
    
    call kepler_solve_RPP(np * (t - t0p), ep, cosf, sinf, j)
    
    r = ap * (1.d0 - ep * ep) / (1.d0 + ep * cosf)
    cosfw = cwp * cosf - swp * sinf
    sinfw = swp * cosf + sinf * cwp
    xp = -r * cosfw
    yp = -r * sinfw * cip
    zp = r * sinfw * Sin(ip)
    
    call kepler_solve_RPP(nm * (t - t0m), em, cosf, sinf, j)
    
    r = am * (1.d0 - em * em) / (1.d0 + em * cosf)
    cosfw = cwm * cosf - swm * sinf
    sinfw = swm * cosf + sinf * cwm
    xm = -r * (comegam * cosfw - somegam * sinfw * cim)
    ym = -r * (somegam * cosfw + comegam * sinfw * cim)
    zm = r * sinfw * Sin(im)
    
end

subroutine grad_coords(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, j, &
    &xp, yp, xm, ym, dxp, dyp, dxm, dym) bind(C, name="grad_coords")

    integer (c_int), bind(C) :: j
    real*8 :: np, nm
    real*8 :: comegam, somegam
    real*8 :: cwp, swp, cip, cwm, swm, cim, sip, sim
    real*8, dimension(j) :: r, cosf, sinf, cosfw, sinfw, denom
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip 
    real (c_double), bind(C) :: am, t0m, em, Pm, Om, wm, im
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), dimension(j), intent(out) :: xp, yp, xm, ym
    real (c_double), bind(C), dimension(j, 13), intent(out) :: dxp, dyp, dxm, dym
    
    real*8 :: np_Pp, nm_Pm
    real*8 :: np_ms, np_mp, np_mm, nm_mp, nm_mm
    real*8, dimension(j) :: r_t0, r_e, r_cosf
    real*8, dimension(j) :: f_M, f_e
    real*8, dimension(j) :: ccmss, cspsc, scpcs, ssmcc
    real*8, dimension(j) :: r_Pp, r_Pm
        
    np = 2 * pi / Pp
    np_Pp = -np / Pp
    
    nm = 2 * pi / Pm
    nm_Pm = -nm / Pm
    
    cwp = Cos(wp)
    swp = Sin(wp)
    cip = Cos(ip)
    sip = Sin(ip)
    
    comegam = Cos(Om)
    somegam = Sin(Om)
    cwm = Cos(wm)
    swm = Sin(wm)
    cim = Cos(im)
    sim = Sin(im)
    
    call grad_kepler_solve_RPP(np * (t - t0p), ep, cosf, sinf, f_e, f_M, j)
    
    denom = 1.d0 / (1.d0 + ep * cosf)
    r = ap * (1.d0 - ep * ep) * denom
    r_cosf = - r * ep * denom
    r_t0 = r_cosf * sinf * f_M * np
    r_e = -r_cosf * sinf * f_e - ap * (cosf + 2 * ep + cosf * ep * ep) * denom * denom
    r_Pp = -r_cosf * sinf * f_M * (t - t0p) * np_Pp
     
    cosfw = cwp * cosf - swp * sinf
    sinfw = swp * cosf + sinf * cwp
    
    ccmss = cosfw
    cspsc = sinfw
    scpcs = sinfw * cip
    ssmcc = - cosfw * cip
    xp = -r * ccmss
    yp = -r * scpcs
    
    dxp(:, 1) = - r * ccmss / ap
    dxm(:, 1) = 0.d0
    dxp(:, 2) = - r_t0 * ccmss - r * f_M * np * cspsc
    dxm(:, 2) = 0.d0
    dxp(:, 3) = - r_e * ccmss + r * f_e * cspsc
    dxm(:, 3) = 0.d0
    dxp(:, 4) = - r_Pp * ccmss + r * f_M * (t - t0p) * np_Pp * cspsc
    dxm(:, 4) = 0.d0
    dxp(:, 5) = r * cspsc
    dxm(:, 5) = 0.d0
    dxp(:, 6) = 0.d0
    dxm(:, 6) = 0.d0
        
    ! ms
    dyp(:, 1) = - r * scpcs / ap
    dym(:, 1) = 0.d0
    ! t0p
    dyp(:, 2) = - r_t0 * scpcs - r * f_M * np * ssmcc
    dym(:, 2) = 0.d0
    ! ep
    dyp(:, 3) = - r_e * scpcs + r * f_e * ssmcc
    dym(:, 3) = 0.d0
    ! Pp
    dyp(:, 4) = - r_Pp * scpcs + r * f_M * (t - t0p) * np_Pp * ssmcc
    dym(:, 4) = 0.d0
    ! wp
    dyp(:, 5) = r * ssmcc
    dym(:, 5) = 0.d0
    ! ip
    dyp(:, 6) = r * sinfw * sip
    dym(:, 6) = 0.d0
    
    call grad_kepler_solve_RPP(nm * (t - t0m), em, cosf, sinf, f_e, f_M, j)
    
    denom = 1.d0 / (1.d0 + em * cosf)
    r = am * (1.d0 - em * em) * denom
    r_cosf = - em * r * denom
    r_t0 = r_cosf * sinf * f_M * nm
    r_e = -r_cosf * sinf * f_e - am * (cosf + 2 * em + cosf * em * em) * denom * denom
    r_Pm = -r_cosf * sinf * f_M * (t - t0m) * nm_Pm
     
    cosfw = cwm * cosf - swm * sinf
    sinfw = swm * cosf + sinf * cwm
    
    ccmss = comegam * cosfw - somegam * sinfw * cim
    cspsc = comegam * sinfw + somegam * cosfw * cim
    scpcs = somegam * cosfw + comegam * sinfw * cim
    ssmcc = somegam * sinfw - comegam * cosfw * cim
    xm = -r * ccmss
    ym = -r * scpcs
        
    ! ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, wm, Om, im, mm
    
    dxp(:, 7) = 0.d0
    dxp(:, 8) = 0.d0
    dxp(:, 9) = 0.d0
    dxp(:, 10) = 0.d0
    dxp(:, 11) = 0.d0
    dxp(:, 12) = 0.d0
    dxp(:, 13) = 0.d0
    
    dyp(:, 7) = 0.d0
    dyp(:, 8) = 0.d0
    dyp(:, 9) = 0.d0
    dyp(:, 10) = 0.d0
    dyp(:, 11) = 0.d0
    dyp(:, 12) = 0.d0
    dyp(:, 13) = 0.d0

    dxm(:, 7) = -r * ccmss / am
    dxm(:, 8) = - r_t0 * ccmss - r * f_M * nm * cspsc
    dxm(:, 9) = - r_e * ccmss + r * f_e * cspsc
    dxm(:, 10) = - r_Pm * ccmss + r * f_M * (t - t0m) * nm_Pm * cspsc
    dxm(:, 11) = r * scpcs
    dxm(:, 12) = r * cspsc
    dxm(:, 13) = -r * somegam * sinfw * sim
    
    dym(:, 7) = - r * scpcs / am
    dym(:, 8) = - r_t0 * scpcs - r * f_M * nm * ssmcc
    dym(:, 9) = - r_e * scpcs + r * f_e * ssmcc
    dym(:, 10) = - r_Pm * scpcs + r * f_M * (t - t0m) * nm_Pm * ssmcc
    dym(:, 11) = -r * ccmss
    dym(:, 12) = r * ssmcc
    dym(:, 13) = r * comegam * sinfw * sim
    
end

subroutine grad_impacts(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, j, &
    &bp, bpm, theta, dbp, dbpm, dtheta) bind(C, name="grad_impacts")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm, bm2, bp2, bpm2, sth
    real*8, dimension(j) :: xp, yp, xm, ym
    real*8, dimension(j) :: bm_xm, bm_ym, bp_xp, bp_yp, bpm_xm, bpm_ym, bpm_xp, bpm_yp
    real*8, dimension(j) :: theta_bp, theta_bpm, theta_bm, denom
    real*8, dimension(j, 13) :: dxp, dyp, dxm, dym, dbm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    real (c_double), bind(C), intent(out), dimension(j, 13) :: dbp, dbpm, dtheta
    
    call grad_coords(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, &
                Pm, Om, wm, im, j, xp, yp, xm, ym, &
                dxp, dyp, dxm, dym)
        
    bm2 = xm**2.d0 + ym**2.d0
    bm = Sqrt(bm2)
    bm_xm = xm / bm
    bm_ym = ym / bm
    
    bp2 = xp**2.d0 + yp**2.d0
    bp = Sqrt(bp2)
    bp_xp = xp / bp
    bp_yp = yp / bp
    
    bpm2 = (xm - xp)**2.d0 + (ym - yp)**2.d0
    bpm = Sqrt(bpm2)
    bpm_xm = (xm - xp) / bpm
    bpm_ym = (ym - yp) / bpm
    bpm_xp = - bpm_xm
    bpm_yp = - bpm_ym
    
    do i=1,j,1
    
        a = bp(i) 
        b = bpm(i)
        c = bm(i)
        
        if (a .GT. b) then
            tmp = b
            b = a
            a = tmp
        end if
        if (b .GT. c) then
            mu = c - (a - b)
        else
            mu = b - (a - c)
        end if
            
        if ((bpm(i) .eq. 0.d0) .or. ((a - b + c) .lt. 1.d-15)) then
            theta(i) = 0.d0
        else if ((a - c + b) .lt. 1.d-15) then 
            theta(i) = pi
        else
            theta(i) = 2 * Atan(Sqrt(((a - b) + c) * mu / ((a + (b + c)) * ((a - c) + b))))
        end if
    
        sth(i) = Sin(theta(i))
        if (abs(sth(i)) .lt. 1.d-15) then 
            theta_bm(i) = 0.d0
            theta_bp(i) = 0.d0
            theta_bpm(i) = 0.d0
        else
            denom(i) = 1.d0 / (bp2(i) * bpm2(i) * sth(i))
            theta_bm(i) = bm(i) * bp(i) * bpm(i) * denom(i)
            theta_bp(i) = ((bpm(i) - bm(i)) * (bpm(i) + bm(i)) - bp2(i)) * 0.5 * bpm(i) * denom(i)
            theta_bpm(i) = ((bp(i) - bm(i)) * (bp(i) + bm(i)) - bpm2(i)) * 0.5 * bp(i) * denom(i)
        end if 
    end do
    
    do i=1,13,1
        dbm(:, i) = bm_xm * dxm(:, i) + bm_ym * dym(:, i)
        dbp(:, i) = bp_xp * dxp(:, i) + bp_yp * dyp(:, i)
        dbpm(:, i) = bpm_xm * dxm(:, i) + bpm_ym * dym(:, i) &
                   + bpm_xp * dxp(:, i) + bpm_yp * dyp(:, i)
        dtheta(:, i) = theta_bm * dbm(:, i) &
                     + theta_bp * dbp(:, i) &
                     + theta_bpm * dbpm(:, i)
    end do

    bp = bp
    bpm = bpm
    dbp = dbp
    dbpm = dbpm
            
end

subroutine impacts(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, j, bp, bpm, theta) bind(C, name="impacts")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm, bm2, bp2, bpm2
    real*8, dimension(j) :: xp, yp, xm, ym, zp, zm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    
    call coords(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, &
                Pm, Om, wm, im, j, xp, yp, zp, xm, ym, zm)
    
    bm2 = xm**2.d0 + ym**2.d0
    bm = Sqrt(bm2)
    
    bp2 = xp**2.d0 + yp**2.d0
    bp = Sqrt(bp2)
    
    bpm2 = (xm - xp)**2.d0 + (ym - yp)**2.d0
    bpm = Sqrt(bpm2)
    
    do i=1,j,1
    
        a = bp(i) 
        b = bpm(i)
        c = bm(i)
        
        if (a .gt. b) then
            tmp = b
            b = a
            a = tmp
        end if
        if (b .gt. c) then
            mu = c - (a - b)
        else
            mu = b - (a - c)
        end if
        
        if ((bpm(i) .lt. 1.d-15) .or. ((a - b + c) .lt. 1.d-15)) then
            theta(i) = 0.d0
        else if ((a - c + b) .lt. 1.d-15) then 
            theta(i) = pi
        else
            theta(i) = 2 * Atan(Sqrt(((a - b) + c) * mu / ((a + (b + c)) * ((a - c) + b))))
        end if
    end do

    bp = bp
    bpm = bpm
            
end

end