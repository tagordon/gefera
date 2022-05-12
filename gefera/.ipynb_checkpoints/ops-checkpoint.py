import aesara
from aesara.graph.op import Op
from aesara.graph.basic import Apply
from aesara.tensor import as_tensor_variable
from aesara import tensor as tt
import numpy as np

import ctypes
from ctypes import byref
clib = ctypes.CDLL("../fortran/wrapper.so")

clib.grad_impacts.restype = None
clib.flux.restype = None

class FluxOp(Op):
    
    """
    An op to compute the lightcurve of a transiting planet/moon system.
    
    Args:
        u1 (scalar): The first limb-darkening parameter
        u2 (scalar): The second limb-darkening parameter
        rp (scalar): The radius of the planet
        rm (scalar): The radius of the moon 
        bp (tensor): The impact parameters for the planet
        bpm (tensor): The planet-moon separation in stellar radii 
        theta (tensor): The angle between the vectors connecting
            the center of the star to the center of the planet and the center 
            of the planet to the center of the moon.
    """
    
    def make_node(self, *args):
        inputs = list(map(as_tensor_variable, args))
        outputs = [
            tt.tensor(
                broadcastable=tuple(inputs[0].broadcastable) + (False,),
                dtype=inputs[0].dtype,
            )
            for _ in range(8)
        ]
        return Apply(self, inputs, outputs)
    
    def perform(self, node, inputs, outputs):
        u1, u2, rp, rm, bp, bpm, theta = inputs
        j = len(bp)
        
        u1, u2, rp, rm = list(
            map(
                lambda a: byref(ctypes.c_double(a)), 
                (u1, u2, rp, rm)
            )
        )
        cth = np.cos(theta)
        sth = np.sin(theta)
        bp, bpm, cth, sth = list(
            map(
                lambda a: byref((ctypes.c_double * j).from_buffer(a)), 
                (bp, bpm, cth, sth)
            )
        )
        
        out = ((ctypes.c_double * 8) * j).from_buffer(np.zeros((8, j)))
        
        clib.flux(u1, u2, rp, rm, bp, bpm, cth, sth, out, byref(ctypes.c_int(j)))
        out = np.array(out).T
        
        # could change parameter order output from clib.flux so that this doesn't look so weird 
        outputs[0][0] = out[0]
        outputs[1][0] = out[6]
        outputs[2][0] = out[7]
        outputs[3][0] = out[1]
        outputs[4][0] = out[2]
        outputs[5][0] = out[3]
        outputs[6][0] = out[4]
        outputs[7][0] = out[5]
        
    def grad(self, inputs, gradients):
        outs = self(*inputs)
        dcdf = gradients[0]
                
        if not isinstance(dcdf.type, aesara.gradient.DisconnectedType):
            return [
                tt.dot(dcdf, outs[1]),
                tt.dot(dcdf, outs[2]),
                tt.dot(dcdf, outs[3]),
                tt.dot(dcdf, outs[4]),
                dcdf * outs[5],
                dcdf * outs[6],
                dcdf * outs[7]
            ]
        else:
            return [
                tt.as_tensor_variable(0.0),
                tt.as_tensor_variable(0.0),
                tt.as_tensor_variable(0.0),
                tt.as_tensor_variable(0.0),
                tt.zeros_like(dcdf),
                tt.zeros_like(dcdf),
                tt.zeros_like(dcdf)
            ]
        
    def R_op(self, inputs, eval_points):
        if eval_points[0] is None:
            return eval_points
        return self.grad(inputs, eval_points)

class ImpactsOp(Op):
    
    """
    An op to compute the inputs required to find the lightcurve for the system. 
    To instantiate the op, an array of times at which to compute the lightcurve 
    should be passed. 
    
    Args:
        ab (scalar): Semimajor axis
        tb (scalar): Time of periastron passage
        eb (scalar): Eccentricity
        pb (scalar): Period
        wb (scalar): Argument of periastron (in radians)
        ib (scalar): Inclination (in radians)
        am (scalar): Semimajor axis
        tm (scalar): Time of periastron passage
        em (scalar): Eccentricity
        pm (scalar): Period
        om (scalar): Longitude of ascending node (in radians)
        wm (scalar): Argument of periastron (in radians)
        im (scalar): Inclination (in radians)
        mm (scalar): Moon/planet mass ratio
        
    Returns:
        bp (tensor): The impact parameters for the planet
        bpm (tensor): The planet-moon separation in stellar radii 
        theta (tensor): The angle between the vectors connecting
            the center of the star to the center of the planet and the center 
            of the planet to the center of the moon.
    """
        
    def __init__(self, t):
        self.t = t
    
    def make_node(self, *args):
        inputs = list(map(as_tensor_variable, args))
        outputs = [
            tt.tensor(
                broadcastable=tuple(inputs[0].broadcastable) + (False,),
                dtype=inputs[0].dtype,
            )
            for _ in range(45)
        ]
        return Apply(self, inputs, outputs)
        
    def perform(self, node, inputs, outputs):        
        j = len(self.t)
        
        bp = (ctypes.c_double * j).from_buffer(np.zeros(j))
        bpm = (ctypes.c_double * j).from_buffer(np.zeros(j))
        theta = (ctypes.c_double * j).from_buffer(np.zeros(j))
        
        dbp = ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14)))
        dbpm = ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14)))
        dtheta = ((ctypes.c_double * j) * 14).from_buffer(np.zeros((j, 14)))

        args = list(map(lambda a: byref(ctypes.c_double(a)), inputs))
        t = byref((ctypes.c_double * j).from_buffer(self.t))
        clib.grad_impacts(t, *args, byref(ctypes.c_int(j)), 
                        bp, bpm, theta, dbp, dbpm, dtheta)
        
        outputs[0][0] = np.array(bp)        
        outputs[1][0] = np.array(bpm)
        outputs[2][0] = np.array(theta)
        
        for i in range(14):
            outputs[3 + i][0] = np.array(dbp[i])
            outputs[17 + i][0] = np.array(dbpm[i])
            outputs[31 + i][0] = np.array(dtheta[i])
                    
    def grad(self, inputs, gradients):
        outs = self(*inputs)
        
        bp = outs[0]
        bpm = outs[1]
        theta = outs[2]
        
        dbp = [outs[3 + i] for i in range(14)]
        dbpm = [outs[17 + i] for i in range(14)]
        dtheta = [outs[31 + i] for i in range(14)]
        
        dcdf = gradients
        
        # might wanna change this to a tensor? 
        g = [0] * 14
        
        for i in range(14):
            if not isinstance(dcdf[0].type, aesara.gradient.DisconnectedType):
                g[i] += tt.dot(dcdf[0], dbp[i]) 
            if not isinstance(dcdf[1].type, aesara.gradient.DisconnectedType):
                g[i] += tt.dot(dcdf[1], dbpm[i])
            if not isinstance(dcdf[2].type, aesara.gradient.DisconnectedType):
                g[i] += tt.dot(dcdf[2], dtheta[i])
                    
        return g
        
    def R_op(self, inputs, eval_points):
        if eval_points[0] is None:
            return eval_points
        return self.grad(inputs, eval_points)