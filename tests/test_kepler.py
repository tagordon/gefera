import pytest
import numpy as np
import gefera as gf

import utils as tu

RPP = gf.kep.RPP
newton = gf.kep.newton

def test_solver():

    n = 100
    M = np.random.rand(n) * np.pi * 2
    ecc = np.random.rand(n)
    ecc[-1] = 0.99999999
    ecc[-2] = 0.0
    
    cosf_RPP, sinf_RPP = RPP(M, ecc)
    cosf_newton, sinf_newton = newton(M, ecc)
    assert np.all(np.isclose(np.array(sinf_RPP), np.array(sinf_newton), atol=1e-12))
    assert np.all(np.isclose(np.array(cosf_RPP), np.array(cosf_newton), atol=1e-12))

def test_impact_conf():
    
    d = tu.random_args_conf()
    o1 = gf.PrimaryOrbit(d['a1'], d['t1'], d['e1'], d['p1'], d['w1'], d['i1'])
    o2 = gf.ConfocalOrbit(d['a2'], d['t2'], d['e2'], d['p2'], d['o2'], d['w2'], d['i2'])
    sys = gf.ConfocalSystem(o1, o2)
    
    t = np.linspace(0, np.max([d['p1'], d['p2']]), 1000)
    (x1, y1, z1), (x2, y2, z2) = sys.coords(t)
    bp, bpm, theta = sys.impacts(t)
    
    # check consistency with law of cosines
    assert np.all(
        np.isclose(
            bpm * bpm + bp * bp - 2 * bp * bpm * np.cos(theta),
            x2 * x2 + y2 * y2,
            atol=1e-10
        )
    )
    
def test_impact_hrch():
    
    d = tu.random_args_hrch()
    o1 = gf.PrimaryOrbit(d['a1'], d['t1'], d['e1'], d['p1'], d['w1'], d['i1'])
    o2 = gf.SatelliteOrbit(d['a2'], d['t2'], d['e2'], d['p2'], d['o2'], d['w2'], d['i2'], d['m2'])
    sys = gf.HierarchicalSystem(o1, o2)
    
    t = np.linspace(0, np.max([d['p1'], d['p2']]), 1000)
    (x1, y1, z1), (x2, y2, z2) = sys.coords(t)
    bp, bpm, theta = sys.impacts(t)
    
    # check consistency with law of cosines
    assert np.all(
        np.isclose(
            bpm * bpm + bp * bp - 2 * bp * bpm * np.cos(theta),
            x2 * x2 + y2 * y2,
            atol=1e-10
        )
    )
    
def test_grad_conf():
    
    args = [
        [8.613763533958783, 7.6911826715470895, 0.810054928381598, 65.0673122426115, 
         6.116839021498093, 0.9090596014452478, 8.522254496506577, 9.253548012338683, 
         0.1353452692423282, 81.69648667461617, 3.527203064379669,1.9466918553511405, 
         0.6229170990214605],
        [6.543808163717847, 1.5307194508070303, 0.5102861338268949, 38.346145696711815, 
         4.752109415235498, 0.3034246098908581, 8.728530164296895, 9.946227259725486, 
         0.5082082737421604, 64.6143918188189, 1.9245179820462204, 0.7250684400125926, 
         1.2429519382855494],
        [4.638367145958105, 9.652138083302301, 0.8258636719163078, 16.65540908205091, 
         1.5355608903323779, 1.5263409783985227, 1.8601526680402882, -9.64638793025562, 
         0.047572401751072455, 21.615350859820236, 1.4760461538247172, 1.7530986787558676, 
         0.26016909353599893],
        [4.69591324353854, 3.0424877618517883, 0.23296401064130862, 86.92699572947448, 
         0.15339956812110442, 0.4439773930744975, 6.5224423731899135, 5.5675225098758485, 
         0.027258381196041775, 71.27304237587627, 2.2453462587899473, 5.423682435643634, 
         0.8202790771897083],
        [5.295774881363932, -1.9240431774039237, 0.3761988177275336, 15.394449929314463, 
         5.15228108580988, 1.3281246873830161, 7.396004566351671, -5.35671764128111, 
         0.33223371889581993, 82.4677110931116, 1.789480468322387, 2.1302746882912555, 
         0.20917047948559375]
    ]
    
    for d in args:
    
        t = np.linspace(0, np.max([d[3], d[9]]), 10)
        dx = 1e-8

        o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
        o2 = gf.ConfocalOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12])
        sys = gf.ConfocalSystem(o1, o2)

        bp, bpm, theta, dbp, dbpm, dtheta = sys.impacts(t, grad=True)
        dbp_fd = np.zeros_like(dbp)
        dbpm_fd = np.zeros_like(dbpm)
        dtheta_fd = np.zeros_like(dtheta)

        for i in range(len(d)):
            d[i] += dx
            o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
            o2 = gf.ConfocalOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12])
            sys_plus = gf.ConfocalSystem(o1, o2)
    
            d[i] -= 2 * dx
            o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
            o2 = gf.ConfocalOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12])
            sys_minus = gf.ConfocalSystem(o1, o2)
    
            bp_plus, bpm_plus, theta_plus = sys_plus.impacts(t)
            bp_minus, bpm_minus, theta_minus = sys_minus.impacts(t)
    
            dbp_fd[i, :] = (bp_plus - bp_minus) / (2 * dx)
            dbpm_fd[i, :] = (bpm_plus - bpm_minus) / (2 * dx)
            dtheta_fd[i, :] = (theta_plus - theta_minus) / (2 * dx)
    
            assert np.all(np.isclose(dbp[i], dbp_fd[i], atol=1e-3))
            assert np.all(np.isclose(dbpm[i], dbpm_fd[i], atol=1e-3))
            assert np.all(np.isclose(dtheta[i], dtheta_fd[i], atol=1e-3))
    
def test_grad_hrch():
    
    args = [
        [5.070367188646622, -5.496009300446632, 0.5110429620244018, 34.74879048392051, 
         5.20142680351299, 0.4407573129135272, 0.03174437986335461, -8.862796356268053, 
         0.11676654421250278, 1.9138712577560262, 2.8424196000737436, 6.207071566685711, 
         0.632034380383273, 0.5453885771964324],
        [4.146602728821974, -4.086248718528058, 0.5216551449954931, 98.20365851526783, 
         2.395354330991065, 1.4791401778726914, 0.4952969026254794, 1.3975986522127233, 
         0.5111893274046421, 1.8361096171456048, 5.09632580117723, 3.543123142737478, 
         0.7789331943669834, 0.2675331736274108],
        [5.220764058789779, -1.4776992859680167, 0.21193936785186743, 49.37897440646385, 
         5.3862664699374365, 0.9079007716446565, 0.8584981425167375, 5.384450401116753, 
         0.5029800185460523, 3.7214579597955386, 5.421324019914887, 1.5128403718857342, 
         0.9315834433265286, 0.3044291013704816],
        [4.872926224905772, 9.774905710338384, 0.6428190703951411, 32.32357812537864, 
         6.190777142137789, 1.180774605518913, 0.5555608738183109, 8.844084790078796, 
         0.16832012349512962, 5.202411454019981, 2.2542656502017766, 3.1407283988638266, 
         1.2849372193668953, 0.10771452678990201],
        [2.9078110378050805, -8.160712480255077, 0.7593436485333901, 20.774078580569, 
         0.0042449085400546906, 0.4494363877542704, 0.2729653902580974, -2.027511209323234, 
         0.20314301082170994, 2.28025396976346, 3.496920437928921, 1.0168878224857407, 
         1.2214339490597386, 0.6117761894719778]
    ]
    
    for d in args:
    
        t = np.linspace(0, np.max([d[3], d[9]]), 10)
        dx = 1e-8

        o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
        o2 = gf.SatelliteOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13])
        sys = gf.HierarchicalSystem(o1, o2)

        bp, bpm, theta, dbp, dbpm, dtheta = sys.impacts(t, grad=True)
        dbp_fd = np.zeros_like(dbp)
        dbpm_fd = np.zeros_like(dbpm)
        dtheta_fd = np.zeros_like(dtheta)

        for i in range(len(d)):
            d[i] += dx
            o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
            o2 = gf.SatelliteOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13])
            sys_plus = gf.HierarchicalSystem(o1, o2)
    
            d[i] -= 2 * dx
            o1 = gf.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
            o2 = gf.SatelliteOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13])
            sys_minus = gf.HierarchicalSystem(o1, o2)
        
            bp_plus, bpm_plus, theta_plus = sys_plus.impacts(t)
            bp_minus, bpm_minus, theta_minus = sys_minus.impacts(t)
    
            dbp_fd[i, :] = (bp_plus - bp_minus) / (2 * dx)
            dbpm_fd[i, :] = (bpm_plus - bpm_minus) / (2 * dx)
            dtheta_fd[i, :] = (theta_plus - theta_minus) / (2 * dx)
    
        assert np.all(np.isclose(dbp[i], dbp_fd[i], atol=1e-4))
        assert np.all(np.isclose(dbpm[i], dbpm_fd[i], atol=1e-4))
        assert np.all(np.isclose(dtheta[i], dtheta_fd[i], atol=1e-4))