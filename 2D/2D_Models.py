import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

'''
Models for testing two population scenarios.
'''

def no_divergence(notused, ns, pts):
    """
    Standard neutral model, populations never diverge.
    """
    
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def no_mig(params, ns, pts):
    """
    Split into two populations, no migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    """
    nu1, nu2, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=0, m21=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def sym_mig(params, ns, pts):
    """
    Split into two populations, with symmetric migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    """
    nu1, nu2, m, T = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    return fs


def asym_mig(params, ns, pts):
    """
    Split into two populations, with different migration rates.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
	"""
    nu1, nu2, m12, m21, T = params
    xx = Numerics.default_grid(pts)
    
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    
    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)
    fs = Spectrum.from_phi(phi, ns, (xx,xx))
    
    return fs    

