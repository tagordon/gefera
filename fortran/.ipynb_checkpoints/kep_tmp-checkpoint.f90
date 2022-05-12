module kepler_tmp
use iso_c_binding
implicit none

real*8, parameter :: pi = 4.d0 * Atan(1.d0), G = 8.887692445125634d-10
real*8, parameter :: o3 = 1.d0 / 3.d0, o6 = 1.d0 / 6.d0, pi_o12 = Atan(1.d0) / 3.d0
real*8, parameter :: au_rsun = 215.03215567054764
real*8, parameter :: if3 = 1.d0 / 6.d0, if5 = 1.d0 / (6.d0 * 20), if7 = 1.d0 / (6.d0 * 20 * 42)
real*8, parameter :: if9 = 1.d0 / (6.d0 * 20 * 42 * 72), if11 = 1.d0 / (6.d0 * 20 * 42 * 72 * 110)
real*8, parameter :: if13 = 1.d0 / (6.d0 * 20 * 42 * 72 * 110 * 156), if15 = 1.d0 / (6.d0 * 20 * 42 * 72 * 110 * 156 * 210)

contains

subroutine kepler_solve(M, ecc, cosf, sinf, j) bind(C, name="kepler_solve")

    integer :: j, i
    real*8 :: tol, err, x, ecc
    real (c_double), bind(C), dimension(j) :: M
    real (c_double), bind(C), dimension(j), intent(out) :: cosf, sinf
    real*8, dimension(j) :: tanfhalf, tanfhalf2, denom, E, sE, cE
    
    if (ecc .lt. 1.d-10) then
        cE = Cos(M)
        sE = Sin(M)
        
        tanfhalf = sE / (1.d0 + cE)
        tanfhalf2 = tanfhalf * tanfhalf
        denom = 1.d0 / (tanfhalf2 + 1.d0)
        cosf = (1.d0 - tanfhalf2) * denom
        sinf = 2 * tanfhalf * denom
    else
        tol = 1.d-15
        
        E = M
        x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
        
        do i=1,j,1
            err = 1.d0
            do while (abs(err) .gt. tol)
                err = - (E(i) - ecc * Sin(E(i)) - M(i)) / (1.d0 - ecc * Cos(E(i)))
                E(i) = E(i) + err
            end do
        end do    
        
        cE = Cos(E)
        tanfhalf = x * Sin(E) / (1.d0 + cE)
        tanfhalf2 = tanfhalf * tanfhalf
        denom = 1.d0 / (tanfhalf2 + 1.d0)
        cosf = (1.d0 - tanfhalf2) * denom
        sinf = 2 * tanfhalf * denom
    end if
end

function sine(x)

    real*8 :: sine
    real*8 :: x, x2
    
    x2 = x * x
    
    sine = x * (1.d0 - x2 * (if3 - &
                    x2 * (if5 - x2 * (if7 - x2 * (if9 - x2 * (if11 - x2 * (if13 - x2 * if15)))))))

end 

! copied from DFM's exoplanet code
function EAstart(M, ecc)

    real*8 :: EAstart
    real*8 :: M, ecc
    real*8 :: ome, sqrt_ome, chi, Lam, S, sigma, s2, s4, denom, E
    
    ome = 1.d0 - ecc
    sqrt_ome = Sqrt(ome)
    chi = M / (sqrt_ome * ome)
    Lam = Sqrt(8.d0 + 9 * chi * chi)
    S = (Lam + 3 * chi) ** (1.d0 / 3.d0)
    sigma = 6 * chi / (2.d0 + S * S + 4.d0 / (S * S))
    s2 = sigma * sigma
    s4 = s2 * s2
    
    denom = 1.d0 / (s2 + 2.d0)
    E = sigma * (1.d0 + s2 * ome * denom * ((s2 + 20.d0) / 60.d0 &
        + s2 * ome * denom * denom * (s2 * s4 + 25 * s4 + 340 * s2 + 840) / 1400))
    EAstart = E * sqrt_ome
end

! pretty much verbatim from DFM's exoplanet code 
subroutine kepler_solve_RPP(M, ecc, cosf, sinf, j) bind(C, name="kepler_solve_RPP")

    integer, bind(C) :: j
    integer :: i, k
    real*8, bind(C) :: ecc
    real*8, dimension(j), bind(C) :: M, cosf, sinf 
    real*8, dimension(j) :: tanfhalf, tanfhalf2, d
    real*8 :: g2s_e, g3s_e, g4s_e, g5s_e, g6s_e
    real*8, dimension(0:12) :: bounds
    real*8, dimension(0:8) :: EA_tab
    real*8 :: MA, EA, sE, cE, x, y
    real*8 :: B0, B1, B2, dx, idx
    integer :: MAsign, sig
    real*8 :: over_ecc
    real*8 :: num, denom, dEA
    
    g2s_e = 0.2588190451025207623489 * ecc
    g3s_e = 0.5 * ecc
    g4s_e = 0.7071067811865475244008 * ecc
    g5s_e = 0.8660254037844386467637 * ecc
    g6s_e = 0.9659258262890682867497 * ecc
    
    bounds(0) = 0.d0
    bounds(1) = pi_o12 - g2s_e
    bounds(2) = pi * o6 - g3s_e
    bounds(3) = pi * 0.25 - g4s_e
    bounds(4) = pi * o3 - g5s_e
    bounds(5) = 5 * pi_o12 - g6s_e
    bounds(6) = pi * 0.5 - ecc
    bounds(7) = 7 * pi_o12 - g6s_e
    bounds(8) = 2 * pi * o3 - g5s_e
    bounds(9) = 3 * pi * 0.25 - g4s_e
    bounds(10) = 5 * pi * o6 - g3s_e
    bounds(11) = 11 * pi_o12 - g2s_e
    bounds(12) = pi
    
    over_ecc = 1.d17
    if (ecc .gt. 1.d-17) then
        over_ecc = 1.d0 / ecc
    end if 
    
    do i=1,j,1
        MAsign = 1
        MA = Modulo(M(i), 2 * pi)
        if (MA .gt. pi) then 
            MAsign = -1
            MA = 2 * pi - MA
        end if
    
        if (2 * MA + 1.d0 - ecc .lt. 0.2d0) then 
            EA = EAstart(MA, ecc)
        else
            do k = 11, 1 , -1
                if (MA .gt. bounds(k)) then
                    exit
                end if
            end do
            EA_tab(0) = k * pi_o12
            EA_tab(6) = (k + 1.d0) * pi_o12
            if (k .ge. 6) then
                sig = 1
            else
                sig = -1
            end if
            
            x = 1.d0 / (1.d0 - ((6.d0 - k) * pi_o12 + sig * bounds(abs(6 - k))))
            y = -0.5 * (k * pi_o12 - bounds(k))
            EA_tab(1) = x
            EA_tab(2) = y * x * x * x
            
            x = 1.d0 / (1.d0 - ((5.d0 - k) * pi_o12 + sig * bounds(abs(5 - k))))
            y = -0.5 * ((k + 1.d0) * pi_o12 - bounds(k + 1))
            EA_tab(7) = x
            EA_tab(8) = y * x * x * x
            
            idx = 1.d0 / (bounds(k + 1) - bounds(k))
            
            B0 = idx * (-EA_tab(2) - idx * (EA_tab(1) - idx * pi_o12))
            B1 = idx * (-2 * EA_tab(2) - idx * (EA_tab(1) - EA_tab(7)))
            B2 = idx * (EA_tab(8) - EA_tab(2))
            
            EA_tab(3) = B2 - 4 * B1 + 10 * B0
            EA_tab(4) = (-2 * B2 + 7 * B1 - 15 * B0) * idx
            EA_tab(5) = (B2 - 3 * B1 + 6 * B0) * idx * idx
            
            dx = MA - bounds(k)
            EA = EA_tab(0) &
                 + dx * (EA_tab(1) + dx * (EA_tab(2) + dx * (EA_tab(3) + dx * (EA_tab(4) + dx * EA_tab(5)))))
                 
        end if
                 
        if (EA .lt. pi * 0.25) then
            sE = sine(EA)
            cE = Sqrt(1.d0 - sE * sE)
        else if (EA .gt. 3 * pi * 0.25) then 
            sE = sine(pi - EA)
            cE = - Sqrt(1.d0 - sE * sE)
        else
            cE = sine(pi * 0.5 - EA)
            sE = Sqrt(1.d0 - cE * cE)
        end if
            
        num = (MA - EA) * over_ecc + sE
        denom = over_ecc - cE
        dEA = num * denom / (denom * denom + 0.5 * sE * num)
            
        if ((ecc .lt. 0.78d0) .OR. (MA .gt. 0.4d0)) then 
            sinf(i) = MAsign * (sE * (1.d0 - 0.5 * dEA * dEA) + dEA * cE)
            cosf(i) = cE * (1.d0 - 0.5 * dEA * dEA) - dEA * sE
        else
            dEA = num / (denom + dEA * (0.5 * sE + o6 * cE * dEA))
            sinf(i) = MAsign * (sE * (1.d0 - 0.5 * dEA * dEA) + dEA * cE * (1.d0 - dEA * dEA * o6))
            cosf(i) = cE * (1.d0 - 0.5 * dEA * dEA) - dEA * sE * (1.d0 - dEA * dEA * o6)
        end if
    end do
    
    x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
    tanfhalf = x * sinf / (1.d0 + cosf)
    tanfhalf2 = tanfhalf * tanfhalf
    d = 1.d0 / (tanfhalf2 + 1.d0)
    cosf = (1.d0 - tanfhalf2) * d
    sinf = 2 * tanfhalf * d
    
end

subroutine grad_kepler_solve(M, ecc, cosf, sinf, f_e, f_M, j) bind(C, name="grad_kepler_solve")

    integer :: j, i
    real*8 :: tol, err, x, ecc
    real (c_double), bind(C), dimension(j) :: M
    real (c_double), bind(C), dimension(j), intent(out) :: cosf, sinf, f_e, f_M
    real*8, dimension(j) :: tanfhalf, tanfhalf2, denom, E, sE, cE
    
    if (ecc .lt. 1.d-10) then
        cE = Cos(M)
        sE = Sin(M)
        
        tanfhalf = sE / (1.d0 + cE)
        tanfhalf2 = tanfhalf * tanfhalf
        denom = 1.d0 / (tanfhalf2 + 1.d0)
        cosf = (1.d0 - tanfhalf2) * denom
        sinf = 2 * tanfhalf * denom
        
        f_M = 1.d0
        f_e = 2.d0 * sinf
    else
        tol = 1.d-10
        
        E = M
        x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
        
        do i=1,j,1
            err = 10.d0
            do while (abs(err) .gt. tol)
                err = - (E(i) - ecc * Sin(E(i)) - M(i)) / (1.d0 - ecc * Cos(E(i)))
                E(i) = E(i) + err
            end do
        end do    
        
        cE = Cos(E)
        tanfhalf = x * Sin(E) / (1.d0 + cE)
        tanfhalf2 = tanfhalf * tanfhalf
        denom = 1.d0 / (tanfhalf2 + 1.d0)
        cosf = (1.d0 - tanfhalf2) * denom
        sinf = 2 * tanfhalf * denom
        
        x = 1.d0 - ecc**2.d0
        f_M = (1.d0 + ecc * cosf)**2.d0 / x**1.5d0
        f_e = (2.d0 + ecc * cosf) * sinf / x
    end if
end

! pretty much verbatim from DFM's exoplanet code 
subroutine grad_kepler_solve_RPP(M, ecc, cosf, sinf, f_e, f_M, j)

    integer :: j, i, k
    real*8 :: ecc
    real*8, dimension(j) :: M, cosf, sinf, tanfhalf, tanfhalf2, d
    real*8, dimension(j) :: f_e, f_M
    real*8 g2s_e, g3s_e, g4s_e, g5s_e, g6s_e
    real*8, dimension(0:12) :: bounds
    real*8, dimension(0:8) :: EA_tab
    real*8 :: MA, EA, sE, cE, x, y
    real*8 :: B0, B1, B2, dx, idx
    integer :: MAsign, sig
    real*8 :: over_ecc
    real*8 :: num, denom, dEA
    
    g2s_e = 0.2588190451025207623489 * ecc
    g3s_e = 0.5 * ecc
    g4s_e = 0.7071067811865475244008 * ecc
    g5s_e = 0.8660254037844386467637 * ecc
    g6s_e = 0.9659258262890682867497 * ecc
    
    bounds(0) = 0.d0
    bounds(1) = pi_o12 - g2s_e
    bounds(2) = pi * o6 - g3s_e
    bounds(3) = pi * 0.25 - g4s_e
    bounds(4) = pi * o3 - g5s_e
    bounds(5) = 5 * pi_o12 - g6s_e
    bounds(6) = pi * 0.5 - ecc
    bounds(7) = 7 * pi_o12 - g6s_e
    bounds(8) = 2 * pi * o3 - g5s_e
    bounds(9) = 3 * pi * 0.25 - g4s_e
    bounds(10) = 5 * pi * o6 - g3s_e
    bounds(11) = 11 * pi_o12 - g2s_e
    bounds(12) = pi
    
    over_ecc = 1.d17
    if (ecc .gt. 1.d-17) then
        over_ecc = 1.d0 / ecc
    end if 
    
    do i=1,j,1
        MAsign = 1
        MA = Modulo(M(i), 2 * pi)
        if (MA .gt. pi) then 
            MAsign = -1
            MA = 2 * pi - MA
        end if
    
        if (2 * MA + 1.d0 - ecc .lt. 0.2d0) then 
            EA = EAstart(MA, ecc)
        else
            do k = 11, 1 , -1
                if (MA .gt. bounds(k)) then
                    exit
                end if
            end do
            EA_tab(0) = k * pi_o12
            EA_tab(6) = (k + 1.d0) * pi_o12
            if (k .ge. 6) then
                sig = 1
            else
                sig = -1
            end if
            
            x = 1.d0 / (1.d0 - ((6.d0 - k) * pi_o12 + sig * bounds(abs(6 - k))))
            y = -0.5 * (k * pi_o12 - bounds(k))
            EA_tab(1) = x
            EA_tab(2) = y * x * x * x
            
            x = 1.d0 / (1.d0 - ((5.d0 - k) * pi_o12 + sig * bounds(abs(5 - k))))
            y = -0.5 * ((k + 1.d0) * pi_o12 - bounds(k + 1))
            EA_tab(7) = x
            EA_tab(8) = y * x * x * x
            
            idx = 1.d0 / (bounds(k + 1) - bounds(k))
            
            B0 = idx * (-EA_tab(2) - idx * (EA_tab(1) - idx * pi_o12))
            B1 = idx * (-2 * EA_tab(2) - idx * (EA_tab(1) - EA_tab(7)))
            B2 = idx * (EA_tab(8) - EA_tab(2))
            
            EA_tab(3) = B2 - 4 * B1 + 10 * B0
            EA_tab(4) = (-2 * B2 + 7 * B1 - 15 * B0) * idx
            EA_tab(5) = (B2 - 3 * B1 + 6 * B0) * idx * idx
            
            dx = MA - bounds(k)
            EA = EA_tab(0) &
                 + dx * (EA_tab(1) + dx * (EA_tab(2) + dx * (EA_tab(3) + dx * (EA_tab(4) + dx * EA_tab(5)))))
                 
        end if
                 
        if (EA .lt. pi * 0.25) then
            sE = sine(EA)
            cE = Sqrt(1.d0 - sE * sE)
        else if (EA .gt. 3 * pi * 0.25) then 
            sE = sine(pi - EA)
            cE = - Sqrt(1.d0 - sE * sE)
        else
            cE = sine(pi * 0.5 - EA)
            sE = Sqrt(1.d0 - cE * cE)
        end if
            
        num = (MA - EA) * over_ecc + sE
        denom = over_ecc - cE
        dEA = num * denom / (denom * denom + 0.5 * sE * num)
            
        if ((ecc .lt. 0.78d0) .OR. (MA .gt. 0.4d0)) then 
            sinf(i) = MAsign * (sE * (1.d0 - 0.5 * dEA * dEA) + dEA * cE)
            cosf(i) = cE * (1.d0 - 0.5 * dEA * dEA) - dEA * sE
        else
            dEA = num / (denom + dEA * (0.5 * sE + o6 * cE * dEA))
            sinf(i) = MAsign * (sE * (1.d0 - 0.5 * dEA * dEA) + dEA * cE * (1.d0 - dEA * dEA * o6))
            cosf(i) = cE * (1.d0 - 0.5 * dEA * dEA) - dEA * sE * (1.d0 - dEA * dEA * o6)
        end if
    end do
    
    x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
    tanfhalf = x * sinf / (1.d0 + cosf)
    tanfhalf2 = tanfhalf * tanfhalf
    d = 1.d0 / (tanfhalf2 + 1.d0)
    cosf = (1.d0 - tanfhalf2) * d
    sinf = 2 * tanfhalf * d
    
    x = 1.d0 - ecc**2.d0
    f_M = (1.d0 + ecc * cosf)**2.d0 / x**1.5d0
    f_e = (2.d0 + ecc * cosf) * sinf / x
end

subroutine coords_sat(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, mm, j, &
    &xp, yp, zp, xm, ym, zm) bind(C, name="coords_sat")
    
    integer (c_int), bind(C) :: j
    real*8 :: np, nm
    real*8 :: mrp, mrm, comegam, somegam
    real*8 :: cwp, swp, cip, cwm, swm, cim
    real*8, dimension(j) :: x, y, z, xbc, ybc, zbc
    real*8, dimension(j) :: r, cosf, sinf, cosfw, sinfw
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
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

    mrp = - 1.d0 / (1.d0 + mm)
    mrm = - mm * mrp
    
    call kepler_solve_RPP(np * (t - t0p), ep, cosf, sinf, j)
    
    r = ap * (1.d0 - ep * ep) / (1.d0 + ep * cosf)
    cosfw = cwp * cosf - swp * sinf
    sinfw = swp * cosf + sinf * cwp
    xbc = -r * cosfw
    ybc = -r * sinfw * cip
    zbc = r * sinfw * Sin(ip)
    
    call kepler_solve_RPP(nm * (t - t0m), em, cosf, sinf, j)
    
    r = am * (1.d0 - em * em) / (1.d0 + em * cosf)
    cosfw = cwm * cosf - swm * sinf
    sinfw = swm * cosf + sinf * cwm
    x = -r * (comegam * cosfw - somegam * sinfw * cim)
    y = -r * (somegam * cosfw + comegam * sinfw * cim)
    z = r * sinfw * Sin(im)
        
    xp = xbc + x * mrm
    yp = ybc + y * mrm
    zp = zbc + z * mrm

    xm = xbc + x * mrp
    ym = ybc + y * mrp
    zm = zbc + z * mrp
    
end

subroutine coords_conf(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, j, &
    &xp, yp, zp, xm, ym, zm) bind(C, name="coords_conf")
    
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

subroutine grad_coords_sat(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, mm, j, &
    &xp, yp, xm, ym, dxp, dyp, dxm, dym) bind(C, name="grad_coords_sat")

    integer (c_int), bind(C) :: j
    real*8 :: np, nm
    real*8 :: mrp, mrm, comegam, somegam
    real*8 :: cwp, swp, cip, cwm, swm, cim, sip, sim
    real*8, dimension(j) :: x, y, xbc, ybc
    real*8, dimension(j) :: r, cosf, sinf, cosfw, sinfw, denom
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), dimension(j), intent(out) :: xp, yp, xm, ym
    
    real (c_double), bind(C), dimension(j, 14), intent(out) :: dxp, dyp, dxm, dym
    real*8, dimension(j) :: xbc_a, ybc_a, x_a, y_a
    
    real*8 :: np_Pp, nm_Pm
    real*8 :: np_ms, np_mp, np_mm, nm_mp, nm_mm
    real*8 :: mrp_mm, mrm_mm, mp
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

    mrp = - 1.d0 / (1.d0 + mm)
    mrm = - mm * mrp
    
    mrp_mm = - mrp / (1.d0 + mm)
    mrm_mm = mrp_mm
    
    mp = 1.d0 / mm
    
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
    xbc = -r * ccmss
    xbc_a = - r * ccmss / ap
    ybc = -r * scpcs
    ybc_a = - r * scpcs / ap
    
    dxp(:, 1) = xbc_a
    dxm(:, 1) = xbc_a
    dxp(:, 2) = - r_t0 * ccmss - r * f_M * np * cspsc
    dxm(:, 2) = dxp(:, 2)
    dxp(:, 3) = - r_e * ccmss + r * f_e * cspsc
    dxm(:, 3) = dxp(:, 3)
    dxp(:, 4) = - r_Pp * ccmss + r * f_M * (t - t0p) * np_Pp * cspsc
    dxm(:, 4) = dxp(:, 4)
    dxp(:, 5) = r * cspsc
    dxm(:, 5) = dxp(:, 5)
    dxp(:, 6) = 0.d0
    dxm(:, 6) = 0.d0
        
    ! ms
    dyp(:, 1) = ybc_a
    dym(:, 1) = ybc_a
    ! t0p
    dyp(:, 2) = - r_t0 * scpcs - r * f_M * np * ssmcc
    dym(:, 2) = dyp(:, 2)
    ! ep
    dyp(:, 3) = - r_e * scpcs + r * f_e * ssmcc
    dym(:, 3) = dyp(:, 3)
    ! Pp
    dyp(:, 4) = - r_Pp * scpcs + r * f_M * (t - t0p) * np_Pp * ssmcc
    dym(:, 4) = dyp(:, 4)
    ! wp
    dyp(:, 5) = r * ssmcc
    dym(:, 5) = dyp(:, 5)
    ! ip
    dyp(:, 6) = r * sinfw * sip
    dym(:, 6) = dyp(:, 6)
    
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
    x = -r * ccmss
    x_a = -r * ccmss / am
    y = -r * scpcs
    y_a = -r * scpcs / am
        
    ! ms, t0p, ep, Pp, Op, wp, ip, mp, t0m, em, Pm, wm, Om, im, mm
    xp = xbc + x * mrm
    
    dxp(:, 7) = x_a * mrm
    dxp(:, 8) = (- r_t0 * ccmss - r * f_M * nm * cspsc) * mrm
    dxp(:, 9) = (- r_e * ccmss + r * f_e * cspsc) * mrm
    dxp(:, 10) = (- r_Pm * ccmss + r * f_M * (t - t0m) * nm_Pm * cspsc) * mrm
    dxp(:, 11) = r * cspsc * mrm
    dxp(:, 12) = r * scpcs * mrm
    dxp(:, 13) = -r * somegam * sinfw * sim * mrm
    dxp(:, 14) = x * mrm_mm
    
    yp = ybc + y * mrm

    dyp(:, 7) = y_a * mrm
    dyp(:, 8) = (- r_t0 * scpcs - r * f_M * nm * ssmcc) * mrm
    dyp(:, 9) = (- r_e * scpcs + r * f_e * ssmcc) * mrm
    dyp(:, 10) = (- r_Pm * scpcs + r * f_M * (t - t0m) * nm_Pm * ssmcc) * mrm
    dyp(:, 11) = r * ssmcc * mrm
    dyp(:, 12) = - r * ccmss * mrm
    dyp(:, 13) = r * comegam * sinfw * sim * mrm
    dyp(:, 14) = y * mrm_mm

    xm = xbc + x * mrp

    dxm(:, 7) = x_a * mrm
    dxm(:, 8) = -dxp(:, 8) * mp
    dxm(:, 9) = -dxp(:, 9) * mp
    dxm(:, 10) = -dxp(:, 10) * mp
    dxm(:, 11) = -dxp(:, 11) * mp
    dxm(:, 12) = -dxp(:, 12) * mp
    dxm(:, 13) = -dxp(:, 13) * mp
    dxm(:, 14) = x * mrp_mm
    
    ym = ybc + y * mrp

    dym(:, 7) = y_a * mrp
    dym(:, 8) = -dyp(:, 8) * mp
    dym(:, 9) = -dyp(:, 9) * mp
    dym(:, 10) = -dyp(:, 10) * mp
    dym(:, 11) = -dyp(:, 11) * mp
    dym(:, 12) = -dyp(:, 12) * mp
    dym(:, 13) = -dyp(:, 13) * mp
    dym(:, 14) = y * mrp_mm
    
end

subroutine grad_coords_conf(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, j, &
    &xp, yp, xm, ym, dxp, dyp, dxm, dym) bind(C, name="grad_coords_conf")

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
    dxm(:, 11) = r * cspsc
    dxm(:, 12) = r * scpcs
    dxm(:, 13) = -r * somegam * sinfw * sim
    
    dym(:, 7) = - r * scpcs / am
    dym(:, 8) = - r_t0 * scpcs - r * f_M * nm * ssmcc
    dym(:, 9) = - r_e * scpcs + r * f_e * ssmcc
    dym(:, 10) = - r_Pm * scpcs + r * f_M * (t - t0m) * nm_Pm * ssmcc
    dym(:, 11) = r * ssmcc
    dym(:, 12) = -r * ccmss
    dym(:, 13) = r * comegam * sinfw * sim
    
end

subroutine grad_impacts_sat(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, mm, j, &
    &bp, bpm, theta, dbp, dbpm, dtheta) bind(C, name="grad_impacts_sat")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm, bm2, bp2, bpm2, sth
    real*8, dimension(j) :: xp, yp, xm, ym
    real*8, dimension(j) :: bm_xm, bm_ym, bp_xp, bp_yp, bpm_xm, bpm_ym, bpm_xp, bpm_yp
    real*8, dimension(j) :: theta_bp, theta_bpm, theta_bm, denom
    real*8, dimension(j, 14) :: dxp, dyp, dxm, dym, dbm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    real (c_double), bind(C), intent(out), dimension(j, 14) :: dbp, dbpm, dtheta
    
    call grad_coords_sat(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, &
                Pm, Om, wm, im, mm, j, xp, yp, xm, ym, &
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
            
        theta(i) = 2 * Atan(Sqrt(((a - b) + c) * mu / ((a + (b + c)) * ((a - c) + b))))
    end do
    
    sth = Sin(theta)
    denom = 1.d0 / (bp2 * bpm2 * sth)
    theta_bm = bm * bp * bpm * denom
    theta_bp = ((bpm - bm) * (bpm + bm) - bp2) * 0.5 * bpm * denom
    theta_bpm = ((bp - bm) * (bp + bm) - bpm2) * 0.5 * bp * denom
    
    do i=1,14,1
        dbm(:, i) = bm_xm * dxm(:, i) + bm_ym * dym(:, i)
        dbp(:, i) = bp_xp * dxp(:, i) + bp_yp * dyp(:, i)
        dbpm(:, i) = bpm_xm * dxm(:, i) + bpm_ym * dym(:, i) &
                   + bpm_xp * dxp(:, i) + bpm_yp * dyp(:, i)
        dtheta(:, i) = theta_bm * dbm(:, i) &
                     + theta_bp * dbp(:, i) &
                     + theta_bpm * dbpm(:, i)
    end do

    bp = bp * au_rsun
    bpm = bpm * au_rsun
    dbp = dbp * au_rsun
    dbpm = dbpm * au_rsun
            
end

subroutine grad_impacts_conf(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, j, &
    &bp, bpm, theta, dbp, dbpm, dtheta) bind(C, name="grad_impacts_conf")

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
    
    call grad_coords_conf(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, &
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
            
        theta(i) = 2 * Atan(Sqrt(((a - b) + c) * mu / ((a + (b + c)) * ((a - c) + b))))
    end do
    
    sth = Sin(theta)
    denom = 1.d0 / (bp2 * bpm2 * sth)
    theta_bm = bm * bp * bpm * denom
    theta_bp = ((bpm - bm) * (bpm + bm) - bp2) * 0.5 * bpm * denom
    theta_bpm = ((bp - bm) * (bp + bm) - bpm2) * 0.5 * bp * denom
    
    do i=1,13,1
        dbm(:, i) = bm_xm * dxm(:, i) + bm_ym * dym(:, i)
        dbp(:, i) = bp_xp * dxp(:, i) + bp_yp * dyp(:, i)
        dbpm(:, i) = bpm_xm * dxm(:, i) + bpm_ym * dym(:, i) &
                   + bpm_xp * dxp(:, i) + bpm_yp * dyp(:, i)
        dtheta(:, i) = theta_bm * dbm(:, i) &
                     + theta_bp * dbp(:, i) &
                     + theta_bpm * dbpm(:, i)
    end do

    bp = bp * au_rsun
    bpm = bpm * au_rsun
    dbp = dbp * au_rsun
    dbpm = dbpm * au_rsun
            
end

subroutine impacts_sat(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, mm, j, bp, bpm, theta) bind(C, name="impacts_sat")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm, bm2, bp2, bpm2
    real*8, dimension(j) :: xp, yp, xm, ym, zp, zm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im, mm
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    
    call coords_sat(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, &
                Pm, Om, wm, im, mm, j, xp, yp, zp, xm, ym, zm)
    
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
        
        if ((bpm(i) .eq. 0.d0) .or. ((a - b + c) .lt. 1.d-15)) then
            theta(i) = 0.d0
        else if ((a - c + b) .lt. 1.d-15) then 
            theta(i) = pi
        else
            theta(i) = 2 * Atan(Sqrt(((a - b) + c) * mu / ((a + (b + c)) * ((a - c) + b))))
        end if
    end do

    bp = bp * au_rsun
    bpm = bpm * au_rsun
            
end

subroutine impacts_conf(t, ap, t0p, ep, Pp, wp, ip, am, &
    &t0m, em, Pm, Om, wm, im, j, bp, bpm, theta) bind(C, name="impacts_conf")

    integer :: i
    real*8 :: a, b, c, tmp, mu
    real*8, dimension(j) :: bm, bm2, bp2, bpm2
    real*8, dimension(j) :: xp, yp, xm, ym, zp, zm
    integer (c_int), bind(C) :: j
    real (c_double), bind(C) :: ap, t0p, ep, Pp, wp, ip, am 
    real (c_double), bind(C) :: t0m, em, Pm, Om, wm, im
    real (c_double), bind(C), dimension(j) :: t
    real (c_double), bind(C), intent(out), dimension(j) :: bp, bpm, theta
    
    call coords_conf(t, ap, t0p, ep, Pp, wp, ip, am, t0m, em, &
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
        
        if ((bpm(i) .eq. 0.d0) .or. ((a - b + c) .lt. 1.d-15)) then
            theta(i) = 0.d0
        else if ((a - c + b) .lt. 1.d-15) then 
            theta(i) = pi
        else
            theta(i) = 2 * Atan(Sqrt(((a - b) + c) * mu / ((a + (b + c)) * ((a - c) + b))))
        end if
    end do

    bp = bp * au_rsun
    bpm = bpm * au_rsun
            
end

end