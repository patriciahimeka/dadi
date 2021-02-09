import numpy
from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum

'''
Models for testing three population scenarios in Wepfer et al., The oceanographic isolation of the Ogasawara Islands and genetic divergence in a reef-building coral.
Adapted from Dan Portik's pipeline (version 2018, python 2.7): https://github.com/dportik/dadi_pipeline/blob/master/README.md#V
Cite: https://onlinelibrary.wiley.com/doi/abs/10.1111/mec.14266

Patricia Wepfer
August, 2018
'''

##########################################################################################
#Basic models of (no gene flow / gene flow) between (all / some) population pairs
##########################################################################################

def split_nomig(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration does not occur between any population pair.
    """
    #6 parameters	
    nu1, nuA, nu2, nu3, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=0, m21=0)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=0, m21=0, m23=0, m32=0, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def split_symmig_all(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is symmetrical between all population pairs (ie 1<->2, 2<->3, and 1<->3).
    """
    #10 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, m3, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=m3, m31=m3)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs


def split_asymmig_all(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3.
    Migration is asymmetrical between population pairs (ie 1<->2 and 1<->3), symmetrical between (2<->3).
    """
    #13 parameters
    nu1, nuA, nu2, nu3, mA1, mA2, m12, m21, m13, m31, m3, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA1, m21=mA2)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m12, m21=m21, m23=m3, m32=m3, m13=m13, m31=m31)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs



def split_symmig_adjacent(params, ns, pts):
    """
    Model with split between pop 1 and (2,3), then split between 2 and 3. Assume 2 occurs
    in between populations 1 and 3, which do not come in to contact with one another.
    Migration is symmetrical between 'adjacent' population pairs (ie 1<->2, 2<->3, but not 1<->3).
    """
    #9 parameters
    nu1, nuA, nu2, nu3, mA, m1, m2, T1, T2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA, m21=mA)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)
    
    phi = Integration.three_pops(phi, xx, T2, nu1=nu1, nu2=nu2, nu3=nu3, m12=m1, m21=m1, m23=m2, m32=m2, m13=0, m31=0)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

def starsplit(params, ns, pts):
    """
    Model with split between pop 1, 2 and 3 at the same time.
    Migration is symmetrical between population pairs (ie 2<->3 (D, O), asymmetrical between (1<->2, 1<->3).
    """
    #9 parameters
    nu1, nu2, nu3, m12, m21, m13, m31, m3, T = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    #phi = Integration.two_pops(phi, xx, T1, nu1=nu1, nu2=nuA, m12=mA1, m21=mA2)

    phi = PhiManip.phi_2D_to_3D_split_2(xx, phi)

    phi = Integration.three_pops(phi, xx, T, nu1=nu1, nu2=nu2, nu3=nu3, m12=m12, m21=m21, m23=m3, m32=m3, m13=m13, m31=m31)

    fs = Spectrum.from_phi(phi, ns, (xx,xx,xx))
    return fs

