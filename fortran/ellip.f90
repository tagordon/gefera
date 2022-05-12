module ellip
implicit none

!real*8 :: el1, el2, el3
real*8, parameter :: pi = 4.d0 * Atan(1.d0)
real*8, parameter :: ln2 = log(2.d0)
real*8, parameter :: pihalf = 2.d0 * Atan(1.d0)

contains

real*8 function el1(x, kc)

    integer :: j
    real*8 :: x, kc, e
    !real, intent(out) :: e
    
    real*8 :: g, m, y
    real*8 :: D = 15
    real*8 :: ca, cb
    !real*8 :: pi = 4.d0 * Atan(1.d0)
    integer :: l, i
    
    ca = 10**(-D/2.d0)
    cb = 10**(-D+2.d0)
    
    y = abs(1.d0/x)
    kc = abs(kc)
    m = 1.d0
    
    if (kc == 0.d0) then
        e = log(x + 1.d0 / cos(atan(x)))
        goto 2
    end if
    
    l = 0
        
1   e = m * kc
    g = m
    m = kc + m
    y = - (e/y) + y
    if (y == 0) then
        y = sqrt(e) * cb
    end if 
    if (abs(g - kc) .gt. ca * g) then
        kc = 2 * sqrt(e)
        l = 2 * l
        if (y .lt. 0) then
            l = l + 1
        end if 
        goto 1
    end if
    if (y .lt. 0) then 
        l = l + 1
    end if 
    e = (atan(m/y) + pi * l) / m
    if (x .lt. 0) then
        e = -e
    end if
    
2 el1 = e
return
end

real*8 function el2(x, kc, a, b)

    real*8 :: x, kc, a, b, e
    !real, intent(out) :: e
    
    real*8 :: c, dd, f, g, ik, m, p, y, z
    real*8 :: D = 15
    real*8 :: ca, cb
    !real*8 :: pi = 4.d0 * Atan(1.d0)
    integer :: l, i
    
    ca = 10**(-D/2.d0)
    cb = 10**(-D+2.d0)
    
    c = x**2
    dd = c + 1.d0
    p = sqrt((1.d0 + kc**2 * c)/dd)
    dd = x / dd
    c = dd / (2.d0 * p)
    z = a - b
    ik = a
    a = (b + a) * 0.5
    y = abs(1.d0/x)
    f = 0.d0
    kc = abs(kc)
    m = 1.d0
        
    if (kc == 0.d0) then
        e = sin(atan(x))
        goto 2
    end if

    l = 0
1   b = ik * kc + b
    e = m * kc
    g = e / p
    dd = f * g + dd
    f = c
    ik = a
    p = g + p
    c = (dd / p + c) / 2.d0
    g = m
    m = kc + m
    a = (b / m + a) * 0.5
    y = - (e/y) + y
    if (y == 0) then
        y = sqrt(e) * cb
    end if 
    if (abs(g - kc) .gt. ca * g) then
        kc = 2 * sqrt(e)
        l = 2 * l
        if (y .lt. 0) then
            l = l + 1
        end if 
        goto 1
    end if
    if (y .lt. 0) then 
        l = l + 1
    end if 
    e = (atan(m/y) + pi * l) * (a / m)
    if (x .lt. 0) then
        e = -e
    else
        e = e + c * z
    end if
    
2 el2 = e
return
end

real*8 function el3(x, kc, p)

    real*8 :: x, kc, p, e
    !real, intent(out) :: e
    
    real*8 :: am, ap, c, d, de, f, fa, g
    real*8 :: h, hh, p1, pm, pz, q, r, s 
    real*8 :: t, u, v, w, y, ye, z, zd
    real*8 :: ca, cb
    real*8 :: cD = 15
    !real*8 :: pi = 4.d0 * Atan(1.d0)
    !real*8 :: ln2 = 0.6931471805599453
    integer :: l, m, n, ND, i, k
    logical :: bo, bk
    real*8 :: ra(10), rb(10), rr(10)
    
    !ca = 1E-6
    !cb = 1E-14
    ca = 10**(-cD/2.d0)
    cb = 10**(-cD+2.d0)
    ND = 10
    
    hh = x * x
    f = p * hh
        
    if (kc == 0.d0) then
        s = ca / (1.d0 + abs(x))
    else
        s = kc
    end if
    t = s * s
    pm = t * 0.5
    e = hh * t
    z = abs(f)
    r = abs(p)
    h = 1.d0 + hh
    if ((e .lt. 1.d0) .AND. (z .lt. 0.1) .AND. (t .lt. 1.d0) .AND. (r .lt. 1.d0)) then
        goto 1
    end if 
    w = 1.d0 + f
    if (w == 0.d0) then
        goto 4
    end if
    if (p == 0.d0) then
        p1 = cb / hh
    else
        p1 = p
    end if
    s = abs(s)
    y = abs(x)
    g = p1 - 1.d0
    if (g == 0.d0) then
        g = cb
    end if
    f = p1 - t
    if (f == 0.d0) then
        f = cb * t
    end if 
    am = 1.d0 - t
    ap = 1.d0 + e
    r = p1 * h
    fa = g / (f * p1)
    bo = fa .gt. 0.d0
    fa = abs(fa)
    pz = abs(g * f)
    de = sqrt(pz)
    q = sqrt(abs(p1))
    if (pm .gt. 0.5) then 
        pm = 0.5
    end if
    pm = p1 - pm
    if (pm .ge. 0.d0) then
        u = sqrt(r * ap)
        v = y * de
        if (g .lt. 0.d0) then
            v = -v
        end if
        d = 1.d0 / q
        c = 1.d0
    else
        u = sqrt(h * ap * pz)
        ye = y * q
        v = am * ye
        q = -de / g
        d = -am / de
        c = 0.d0
        pz = ap - r
    end if
    if (bo) then
        r = v / u
        z = 1.d0
        k = 1
        if (pm .lt. 0.d0) then
            h = y * sqrt(h / (ap * fa))
            h = 1.d0 / h - h
            z = h - r - r
            r = 2.d0 + r * h
            if (r == 0.d0) then
                r = cb
            end if
            if (z == 0.d0) then 
                z = h * cb
            end if
            r = r / z
            z = r
            w = pz
        end if
        u = u / w
        v = v / w
    else
        t = u + abs(v)
        bk = .TRUE.
        if (p1 .lt. 0.d0) then
            de = v / pz
            ye = u * ye
            ye = ye + ye
            u = t / pz
            v = (-f - g * e) / t
            t = pz * abs(w)
            z = (hh * r * f - g * ap + ye) / t
            ye = ye / t
        else
            de = v / w 
            ye = 0.d0
            u = (e + p1) / t
            v = t / w
            z = 1.d0
        end if
        if (s .gt. 1.d0) then
            h = u
            u = v
            v = h
        end if
    end if
    y = 1.d0 / y
    e = s
    n = 1
    t = 1.d0
    m = 0
    l = 0
3   y = y - e / y
    if (y == 0.d0) then
        y = sqrt(e) * cb
    end if
    f = c
    c = d / q + c
    g = e / q
    d = f * g + d
    d = d + d
    q = g + q
    g = t
    t = s + t
    n = n + n
    m = m + m
    if (bo) then 
        if (z .lt. 0.d0) then 
            m = k + m
        end if 
        k = int(sign(1.d0, r))
        h = e / (u * u + v * v)
        u = u * (1.d0 + h)
        v = v * (1.d0 - h)
    else
        r = u / v 
        h = z * r
        z = h * z
        hh = e / v
        if (bk) then
            de = de / u
            ye = ye * (h + 1.d0 / h) + de * (1.d0 + r)
            de = de * (u - hh)
            bk = abs(ye) .lt. 1.d0
        else
            k = exponent(z)
            z = fraction(z)
            m = m + k
        end if
    end if
    if (abs(g - s) .gt. ca * g) then 
        if (bo) then
            g = (1.d0 / r - r) * 0.5
            hh = u + v * g
            h = g * u - v
            if (hh == 0.d0) then
                hh = u * cb
            end if
            if (h == 0.d0) then
                h = v * cb
            end if
            z = r * h
            r = hh / h
        else
            u = u + e / u
            v = v + hh
        end if
        s = sqrt(e)
        s = s + s
        e = s * t
        l = l + l
        if (y .lt. 0.d0) then 
            l = 1 + l
        end if 
        goto 3
    end if
    if (y .lt. 0.d0) then
        l = 1 + l
    end if
    e = atan(t / y) + pi * l
    e = e * (c * t + d) / (t * (t + q))
    if (bo) then 
        h = v / (t + u)
        z = 1.d0 - r * h
        h = r + h
        if (z == 0.d0) then
            z = cb
        end if
        if (z .lt. 0.d0) then
            m = m + int(sign(1.d0, h))
        end if
        s = atan(h / z) + m * pi
    else
        if (bk) then
            s = asinh(ye)
        else
            s = log(z) + m * ln2
        end if
        s = s * 0.5
    end if
    e = (e + sqrt(fa) * s) / n
    if (x .le. 0.d0) then
        e = -e
    end if
    goto 4

1   do k=2,ND,1
        rb(k) = 0.5 / k
        ra(k) = 1.d0 - rb(k)
    end do
    zd = 0.5 / (ND + 1)
    s = p + pm
    do k=2,ND,1
        rr(k) = s
        pm = pm * t * ra(k)
        s = s * p + pm
    end do
    u = s * zd
    s = u
    bo = .FALSE.
    do k=ND,2,-1
        u = u + (rr(k) - u) * rb(k)
        bo = .NOT. bo
        if (bo) then
            v = -u
        else
            v = u
        end if
        s = s * hh + v
    end do
    if (bo) then
        s = -s
    end if
    u = (u + 1.d0) * 0.5
    e = (u - s * h) * sqrt(h) * x + u * asinh(x)
    
4 el3 = e
return
end

real*8 function cel(kc, p, a, b)

    real*8 :: kc, p, a, b
    real*8 :: CA, e, f, g, m, q
    
    CA = 1.d-7
    
    kc = abs(kc)
    e = kc
    m = 1.d0

    if (p .gt. 0.d0) then
        p = sqrt(p)
        b = b / p
    else
        f = kc * kc
        q = 1.d0 - f
        g = 1.d0 - p
        f = f - p
        q = (b - a * p) * q
        p = sqrt(f / g)
        a = (a - b) / g
        b = -q / (g * g * p) + a * p
    end if 

1   f = a
    a = b / p + a
    g = e / p
    b = f * g + b
    b = b + b
    p = g + p
    g = m
    m = kc + m
    
    if (abs(g - kc) .gt. g * CA) then
        kc = sqrt(e)
        kc = kc + kc
        e = kc * m
        goto 1
    end if
    cel = pihalf * ((a * m + b) / (m * (m + p)))

return
end

end module ellip