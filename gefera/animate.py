import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import gefera as gf

def draw(sys, ax, t, r1, r2, ld_params=None, cmap=plt.cm.copper, fill=True):
    
    """
    Draw the state of the system.
    
    Args:
        sys (System): The system
        ax (matplotlib.pyplot.axis): Axis on which to draw
        t (float): Time at which to draw snapshot of system
        r1 (float): Radius of the first body
        r2 (float): Radius of the second body
        ld_params (optional, float): Limb-darkening parameters for the star
        cmap (optional, matplotlib.colors.Colormap): Colormap for coloring the star according to ld_params
        fill (optional, bool): If true, fill transiting bodies, if false draw outlines only
    """
        
    if isinstance(t, np.ndarray):
        raise Exception("Argument t should be a scalar, not an array.")
            
    x1, y1, _, x2, y2, _ = sys.kep.coords(np.array([t]), sys.pdict) 
        
    b1 = plt.Circle(
        (-x1, y1),
        radius=r1,
        color='k',
        fill=fill
    )
    b2 = plt.Circle(
        (-x2, y2),
        radius=r2,
        color='k',
        fill=fill
    )
    if ld_params is None:
        star = plt.Circle((0, 0), radius=1, color=cmap(1.0), fill=True)
        ax.add_patch(star)
        
    else:
        u1, u2 = ld_params
            
        d = 0.01
        x = np.arange(-1.0, 1.0, d)
        y = np.arange(-1.0, 1.0, d)
        x, y = np.meshgrid(x, y)
        shape = np.shape(x)
        x = x.flatten()
        y = y.flatten()
        z = np.zeros(len(x))
        with np.errstate(invalid='ignore'):
            for i, (x, y) in enumerate(zip(x, y)):
                z[i] = 1 - u1 * (1 - np.sqrt(1 - x**2 - y**2)) - u2 * (1 - np.sqrt(1 - x**2 - y**2))**2
        z = z.reshape(shape)
            
        cmap = cmap
        cmap.set_bad('white')

        edge = plt.Circle(
            (0.005, -0.0), 
            radius=1.0, 
            color='w',
            fill=False, 
            linewidth=4
        )
            
        ax.imshow(z, interpolation='bilinear', extent=(-1, 1, -1, 1), cmap=cmap)
        ax.add_patch(edge)
        
    ax.add_patch(b1)
    ax.add_patch(b2)
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_axis_off()
    
def draw_series(sys, ax, t, r1, r2, ld_params=None, cmap=plt.cm.copper, fill=True, alpha=1.0):
            
    if ld_params is None:
        star = plt.Circle((0, 0), radius=1, color=cmap(1.0), fill=True)
        ax.add_patch(star)
        
    else:
        u1, u2 = ld_params
            
        d = 0.01
        x = np.arange(-1.0, 1.0, d)
        y = np.arange(-1.0, 1.0, d)
        x, y = np.meshgrid(x, y)
        shape = np.shape(x)
        x = x.flatten()
        y = y.flatten()
        z = np.zeros(len(x))
        with np.errstate(invalid='ignore'):
            for i, (x, y) in enumerate(zip(x, y)):
                z[i] = 1 - u1 * (1 - np.sqrt(1 - x**2 - y**2)) - u2 * (1 - np.sqrt(1 - x**2 - y**2))**2
        z = z.reshape(shape)
            
        cmap = cmap
        cmap.set_bad('white')

        edge = plt.Circle(
            (0.005, -0.0), 
            radius=1.0, 
            color='w',
            fill=False, 
            linewidth=4
        )
            
        ax.imshow(z, interpolation='bilinear', extent=(-1, 1, -1, 1), cmap=cmap)
        ax.add_patch(edge)
    
    for ti in t:
        x1, y1, _, x2, y2, _ = sys.kep.coords(np.array([ti]), sys.pdict) 
        
        b1 = plt.Circle(
            (-x1, y1),
            radius=r1,
            color='k',
            fill=fill,
            alpha=alpha
        )
        b2 = plt.Circle(
            (-x2, y2),
            radius=r2,
            color='k',
            fill=fill,
            alpha=alpha
        )
        
        ax.add_patch(b1)
        ax.add_patch(b2)
        
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
    ax.set_axis_off()
    
def snapshots(sys, axs, t, r1, r2, ld_params=None, cmap=plt.cm.copper, fill=True):
    
    """
    Draw a series of snapshots representing the system and multiple times.
    
    Args:
        sys (System): The system
        axs (list): A list of axes of length len(t) on which to plot the snapshots
        t (array): Time at which to draw snapshot of system
        r1 (float): Radius of the first body
        r2 (float): Radius of the second body
        ld_params (optional, float): Limb-darkening parameters for the star
        cmap (optional, matplotlib.colors.Colormap): Colormap for coloring the star according to ld_params
        fill (optional, bool): If true, fill transiting bodies, if false draw outlines only
    """
    
    if len(axs) != len(t):
        raise ValueError('One axis object is required for each time point.')
    
    for ax, ti in zip(axs, t):
        draw(sys, ax, ti, r1, r2, ld_params=ld_params, cmap=cmap, fill=fill)
        
def animate(sys, fig, t, r1, r2, duration=5, ld_params=None, cmap=plt.cm.copper):
    
    """
    Make an animation of the system.
    
    Args:
        sys (System): The system
        fig (matplotlib.pyplot.figure): Figure object
        t (array): Array of times representing each frame of the animation
        r1 (float): Radius of the first body
        r2 (float): Radius of the second body
        ld_params (optional, float): Limb-darkening parameters for the star
        cmap (optional, matplotlib.colors.Colormap): Colormap for coloring the star according to ld_params
        
    Returns:
        matplotlib.animation.FuncAnimation: The animation.
    """
        
    t = np.ascontiguousarray(t)
    ax = fig.gca()
            
    x1, y1, z1, x2, y2, z2 = sys.kep.coords(t, sys.pdict)
        
    if ld_params is None:
        star = plt.Circle((0, 0), radius=1, color=cmap(1.0), fill=True)
        ax.add_patch(star)
        
    else:
        u1, u2 = ld_params
            
        d = 0.01
        x = np.arange(-1.0, 1.0, d)
        y = np.arange(-1.0, 1.0, d)
        x, y = np.meshgrid(x, y)
        shape = np.shape(x)
        x = x.flatten()
        y = y.flatten()
        z = np.zeros(len(x))
        with np.errstate(invalid='ignore'):
            for i, (x, y) in enumerate(zip(x, y)):
                z[i] = 1 - u1 * (1 - np.sqrt(1 - x**2 - y**2)) - u2 * (1 - np.sqrt(1 - x**2 - y**2))**2
        z = z.reshape(shape)
            
        cmap = cmap
        cmap.set_bad('white')

        edge = plt.Circle(
            (0.005, -0.0), 
            radius=1.0, 
            color='w',
            fill=False, 
            linewidth=4
        )
        
        ax.imshow(z, interpolation='bilinear', extent=(-1, 1, -1, 1), cmap=cmap)
        ax.add_patch(edge)
        
    interval = int(duration * 1000 / len(t))
        
    p_patch = ax.add_patch(plt.Circle((0, 0), 0, color='k'))
    m_patch = ax.add_patch(plt.Circle((0, 0), 0, color='k'))
        
    ax.set_xlim(-2, 2)
    ax.set_ylim(-2, 2)
        
    ax.set_axis_off()
        
    def init():

        p_patch.set_center((0, 0))
        p_patch.set_radius(r1)
        m_patch.set_center((0, 0))
        m_patch.set_radius(r2)
        return p_patch, m_patch
        
    def update(i):

        p_patch.set_center((-x1[i], y1[i]))
        m_patch.set_center((-x2[i], y2[i]))
        return p_patch, m_patch
        
    return animation.FuncAnimation(
        fig, 
        update, 
        frames=np.arange(len(t)), 
        init_func=init, 
        blit=True, 
        interval=interval
    )