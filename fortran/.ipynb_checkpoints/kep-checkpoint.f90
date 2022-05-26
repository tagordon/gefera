module kepler
use iso_c_binding
implicit none

real*8, parameter :: pi = 4.d0 * Atan(1.d0), G = 8.887692445125634d-10
real*8, parameter :: o3 = 1.d0 / 3.d0, o6 = 1.d0 / 6.d0, pi_o12 = Atan(1.d0) / 3.d0
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
         tol = 1.d-10

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

    integer :: j
    integer :: i, k
    real*8, bind(C) :: ecc
    real*8, dimension(j), bind(C) :: M, cosf, sinf 
    real*8 :: tanfhalf, tanfhalf2, d
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
        
        denom = 1.d0 + cosf(i)
        if (denom .gt. 1.d-10) then 
            x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
            tanfhalf = x * sinf(i) / (1.d0 + cosf(i))
            tanfhalf2 = tanfhalf * tanfhalf
            d = 1.d0 / (tanfhalf2 + 1.d0)
            cosf(i) = (1.d0 - tanfhalf2) * d
            sinf(i) = 2 * tanfhalf * d
        else
            cosf(i) = -1.d0
            sinf(i) = 0.d0
        end if
        
    end do
    
end

! pretty much verbatim from DFM's exoplanet code 
subroutine grad_kepler_solve_RPP(M, ecc, cosf, sinf, f_e, f_M, j)

    integer :: j, i, k
    real*8 :: ecc
    real*8, dimension(j) :: M, cosf, sinf
    real*8 :: tanfhalf, tanfhalf2, d
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
        
        denom = 1.d0 + cosf(i)
        if (denom .gt. 1.d-10) then 
            x = Sqrt((1.d0 + ecc) / (1.d0 - ecc))
            tanfhalf = x * sinf(i) / (1.d0 + cosf(i))
            tanfhalf2 = tanfhalf * tanfhalf
            d = 1.d0 / (tanfhalf2 + 1.d0)
            cosf(i) = (1.d0 - tanfhalf2) * d
            sinf(i) = 2 * tanfhalf * d
        else
            cosf(i) = -1.d0
            sinf(i) = 0.d0
        end if
        
    end do
    
    x = 1.d0 - ecc**2.d0
    f_M = (1.d0 + ecc * cosf)**2.d0 / x**1.5d0
    f_e = (2.d0 + ecc * cosf) * sinf / x
end

end