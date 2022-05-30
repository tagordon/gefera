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
    o1 = gf.orbits.PrimaryOrbit(d['a1'], d['t1'], d['e1'], d['p1'], d['w1'], d['i1'])
    o2 = gf.orbits.ConfocalOrbit(d['a2'], d['t2'], d['e2'], d['p2'], d['o2'], d['w2'], d['i2'])
    sys = gf.systems.ConfocalSystem(o1, o2)
    
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
    o1 = gf.orbits.PrimaryOrbit(d['a1'], d['t1'], d['e1'], d['p1'], d['w1'], d['i1'])
    o2 = gf.orbits.SatelliteOrbit(d['a2'], d['t2'], d['e2'], d['p2'], d['o2'], d['w2'], d['i2'], d['m2'])
    sys = gf.systems.HierarchicalSystem(o1, o2)
    
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
        [348.2170441566445, 7.307628034226763, 0.8228888917492587, 44.88210852044057, 
         5.835757651640773, 1.3604059438056473, 79.93434712267178, 8.546206641861492, 
         0.5517615558457643, 87.62949992459167, 0.14653478108919407, 3.3625532506160436, 
         0.9337756404381136],
        [488.1493132357637, -3.0802965763826347, 0.6085074457928773, 46.55996385963626, 
         1.9905924701131195, 0.7663172484294362, 494.36179021626776, 3.1604782371959104, 
         0.17032527638132208, 18.533464379848063, 4.684024361474528, 0.03722382515289265, 
         0.5627893681157864],
        [271.74383871372356, 9.065441614111315, 0.035549456295412805, 12.34715387019951, 
         0.6921819980790107, 1.2470043497235688, 269.38064060956384, -1.2755597561147685, 
         0.4064394144194549, 4.844329641707633, 1.6139028876341346, 4.807773800053764, 
         0.38568116805602864],
        [306.2383793905593, 6.787311603481054, 0.8375364075208153, 51.738580027961156, 
         0.9332360851851078, 1.1301872420133285, 24.911748976774103, -5.82934339171941, 
         0.4303049597431203, 52.57265444712054, 0.1481945023967446, 2.72256992017862, 
         1.0439735478650642],
        [126.49577568012296, -1.6887551108507317, 0.7887703304444428, 27.324435848357574, 
         2.300698662352017, 0.6358422399132542, 373.1833389804012, -2.722821719341704, 
         0.5009328156813628, 1.5256607723395088, 0.7942670794023093, 3.903946851711146, 
         0.5755550294555625]
    ]
    
    for d in args:
    
        t = np.linspace(0, np.max([d[3], d[9]]), 10)
        dx = 1e-8

        o1 = gf.orbits.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
        o2 = gf.orbits.ConfocalOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12])
        sys = gf.systems.ConfocalSystem(o1, o2)

        bp, bpm, theta, dbp, dbpm, dtheta = sys.impacts(t, grad=True)
        dbp_fd = np.zeros_like(dbp)
        dbpm_fd = np.zeros_like(dbpm)
        dtheta_fd = np.zeros_like(dtheta)

        for i in range(len(d)):
            d[i] += dx
            o1 = gf.orbits.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
            o2 = gf.orbits.ConfocalOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12])
            sys_plus = gf.systems.ConfocalSystem(o1, o2)
    
            d[i] -= 2 * dx
            o1 = gf.orbits.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
            o2 = gf.orbits.ConfocalOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12])
            sys_minus = gf.systems.ConfocalSystem(o1, o2)
    
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
        [230.4361514495303, -6.284444848735305, 0.8573851643114856, 68.00723789962558, 
         0.9162470417186847, 0.9178826926279238, 2.77799446175148, -6.172556210928488, 
         0.1349734654411347, 3.329156499672567, 5.298506546855079, 5.0835284537239, 
         0.45203857280226817, 0.8224527219851816],
        [138.51719650310912, -5.235504192860693, 0.37459635046317974, 90.96225459983489, 
         4.7445238730507615, 0.2502317144568418, 18.81770309161984, 7.2210919571100725, 
         0.8154787720817648, 1.0620475597178445, 0.572819535052827, 1.5491110992431238, 
         0.9374472582027031, 0.42764275977666233],
        [281.740033959577, 6.791226631754277, 0.5426955367361197, 82.74887086913452, 
         3.5209799910964343, 0.9693576873940464, 1.861550968547006, 2.9086339675184583, 
         0.16563389670329337, 8.289225238532481, 2.051018848374333, 4.751376938946587, 
         1.2259328144119372, 0.10505253193378872],
        [256.7457907711693, -3.78314822663315, 0.6556335056791278, 62.62352939216658, 
         4.239748874408342, 1.196723912140346, 15.321391331025593, -3.748802161120195, 
         0.12818322331487872, 3.5355244410352817, 3.3411939939674724, 1.93635780128358, 
         0.24604918439230836, 0.3992063390153622],
        [42.97037482843434, -3.570552097422272, 0.26299679980815677, 25.614888087169746, 
         0.6766356646206445, 0.4396619723203977, 5.5822612129743625, -0.19661182243682518, 
         0.9905306264801919, 8.239183971371848, 5.304137747469186, 1.8214306644469704, 
         0.4019925429852727, 0.7749014230664877]
    ]
    
    for d in args:
    
        t = np.linspace(0, np.max([d[3], d[9]]), 10)
        dx = 1e-8

        o1 = gf.orbits.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
        o2 = gf.orbits.SatelliteOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13])
        sys = gf.systems.HierarchicalSystem(o1, o2)

        bp, bpm, theta, dbp, dbpm, dtheta = sys.impacts(t, grad=True)
        dbp_fd = np.zeros_like(dbp)
        dbpm_fd = np.zeros_like(dbpm)
        dtheta_fd = np.zeros_like(dtheta)

        for i in range(len(d)):
            d[i] += dx
            o1 = gf.orbits.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
            o2 = gf.orbits.SatelliteOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13])
            sys_plus = gf.systems.HierarchicalSystem(o1, o2)
    
            d[i] -= 2 * dx
            o1 = gf.orbits.PrimaryOrbit(d[0], d[1], d[2], d[3], d[4], d[5])
            o2 = gf.orbits.SatelliteOrbit(d[6], d[7], d[8], d[9], d[10], d[11], d[12], d[13])
            sys_minus = gf.systems.HierarchicalSystem(o1, o2)
        
            bp_plus, bpm_plus, theta_plus = sys_plus.impacts(t)
            bp_minus, bpm_minus, theta_minus = sys_minus.impacts(t)
    
            dbp_fd[i, :] = (bp_plus - bp_minus) / (2 * dx)
            dbpm_fd[i, :] = (bpm_plus - bpm_minus) / (2 * dx)
            dtheta_fd[i, :] = (theta_plus - theta_minus) / (2 * dx)
    
        assert np.all(np.isclose(dbp[i], dbp_fd[i], atol=1e-4))
        assert np.all(np.isclose(dbpm[i], dbpm_fd[i], atol=1e-4))
        assert np.all(np.isclose(dtheta[i], dtheta_fd[i], atol=1e-4))