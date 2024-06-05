import numpy as np
import astropy.units as u
import astropy.constants as c
from scipy.optimize import fsolve
from scipy.interpolate import UnivariateSpline, interp1d
import matplotlib.pyplot as plt
import lal
import lalsimulation as lalsim

#from .data import PACKAGE_FILENAMES


def compute_isco(chi_bh):
    '''
    This function takes as input the aligned spin component of the BH and
    returns the innermost stable circular orbit radius normalized by the mass
    of the BH.
    '''
    if not np.all([np.abs(chi_bh) <= 1.0]):
        raise ValueError('|chi1| must be less than 1.0')
    Z1 = 1.0 + ((1.0 - chi_bh**2)**(1. / 3.)) * ((1 + chi_bh)**(1. / 3.) +
                                             (1 - chi_bh)**(1. / 3.))
    Z2 = np.sqrt(3 * chi_bh**2 + Z1**2)
    r_isco = 3 + Z2 - np.sign(chi_bh) * np.sqrt((3 - Z1) * (3 + Z1 + 2 * Z2))
    return r_isco

def spins_to_chi_eff(m1, m2, chi1, chi2):
    '''
    This function returns the effective spin of the system
    '''
    chi_eff = (chi1*m1 + chi2*m2)/(m1 + m2)
    return chi_eff


def _lalsim_neutron_star_radius(m, max_mass, fam):
    '''Return neutron star radius in meters'''
    if m > max_mass:
        return m * 1500  # 1 solar mass = 1500m
    else:
        try:
            return lalsim.SimNeutronStarRadius(m * lal.MSUN_SI, fam)
        except RuntimeError:
            # FIXME handle RuntimeError for edge cases
            # Raised if m is close to max_mass
            # GSL function failed: interpolation error (errnum=1)
            # XLAL Error - <GSL function> (interp.c:150): Generic failure
            return m * 1500


def _r_ns_from_lal_simulation(m_ns, eosname):
    eos = lalsim.SimNeutronStarEOSByName(eosname)
    fam = lalsim.CreateSimNeutronStarFamily(eos)
    max_mass = lalsim.SimNeutronStarMaximumMass(fam) / c.M_sun.value
    try:
        iter(m_ns)
        R_ns = np.array([_lalsim_neutron_star_radius(m, max_mass, fam)
                        for m in m_ns])
    except TypeError:
        R_ns = _lalsim_neutron_star_radius(m_ns, max_mass, fam)
    return R_ns, max_mass


def _compactness_baryon_mass(m_ns, r_ns):
    C_ns = c.G * m_ns * c.M_sun / (r_ns * u.m * c.c**2)
    C_ns = C_ns.value

    d1 = 0.619
    d2 = 0.1359
    BE = (d1 * C_ns + d2 * C_ns * C_ns) * m_ns  # arXiv:1601.06083 (Eq: 21)
    m2_b = m_ns + BE  # Baryonic mass - Gravitational mass = Binding Energy
    return [C_ns, m2_b]


def max_mass_from_eosname(eosname):
    if eosname == "2H":
        max_mass = 2.834648092299807
    else:
        eos = lalsim.SimNeutronStarEOSByName(eosname)
        fam = lalsim.CreateSimNeutronStarFamily(eos)
        max_mass = lalsim.SimNeutronStarMaximumMass(fam) / c.M_sun.value
    return max_mass


def computeCompactness(M_ns, eosname='2H', max_mass=None):
    '''
    Return the neutron star compactness as a function of mass
    and equation of state or radius

    Parameters
    ----------
    M_ns : array_like
        Neutron star mass in solar masses
    eosname : str or interp1d
        Neutron star equation of state to be used
    max_mass : float
        Maximum mass of neutron star.

    Returns
    -------
    [C_ns, m2_b, max_mass]
        Compactness, baryon mass and maximum neutron star mass
        in solar masses.

    Notes
    -----
    The radius and maximum mass of the neutron star is
    inferred based on the equation of state supplied.
    Max mass only needs to be supplied for EoS marginalization.

    Examples
    --------
    >>> computeCompactness(2.8)
    [array(0.298), array(3.354), 2.834]
    >>> computeDiskMass.computeCompactness(2.9, eosname='AP4')
    [0.5, 0.0, 2.212]
    >>> m_ns = np.array([1.1, 1.2, 1.3])
    >>> computeDiskMass.computeCompactness(m_ns, eosname='AP4')
    [array([0.141, 0.154, 0.167]), array([1.199, 1.318, 1.439]), 2.212]
    '''
    if isinstance(eosname, interp1d):
        # find R as a function of M
        R_ns = eosname(M_ns)
        C_ns, m2_b = _compactness_baryon_mass(M_ns, R_ns)
    elif eosname != '2H':
        # infer radius and maximum mass based on lalsimulation EoS
        R_ns, max_mass = _r_ns_from_lal_simulation(M_ns, eosname)
        C_ns, m2_b = _compactness_baryon_mass(M_ns, R_ns)
    else:
        with open(PACKAGE_FILENAMES['equil_2H.dat'], 'rb') as f:
            M_g, M_b, Compactness = np.loadtxt(f, unpack=True)
        max_mass = max_mass_from_eosname("2H")
        s = UnivariateSpline(M_g, Compactness, k=5)
        s_b = UnivariateSpline(M_g, M_b, k=5)
        C_ns = s(M_ns)
        m2_b = s_b(M_ns)
    try:
        C_ns[M_ns > max_mass] = 0.5  # BH compactness set to 0.5
        m2_b[M_ns > max_mass] = 0.0  # BH baryon mass set to 0.0
    except TypeError:  # if C_ns is not an array
        if M_ns > max_mass:
            C_ns = 0.5
            m2_b = 0.0
    return [C_ns, m2_b, max_mass]

def zeta_diskejecta(Q):
	"""
	Compute viscous outflowas as a function of the mass ratio
	Equation (12) from arxiv 2102.11569
	"""
	zeta1_lower=0.04
	zeta2_lower=0.14
	zeta1_upper=0.32
	zeta2_upper=0.44
	zeta_lower=zeta1_lower + (zeta2_lower - zeta1_lower) / (1 + np.exp(1.5*(Q-3)))
	zeta_upper=zeta1_upper + (zeta2_upper - zeta1_upper) / (1 + np.exp(1.5*(Q-3)))
	return zeta_lower, zeta_upper
	
'''
For Newtonian: Obtained from arXiv 1807.00011
alpha = 0.406
beta = 0.139
gamma = 0.255
delta = 1.761

For Kerr: Obtained from private message from F. Foucart
alpha = 0.24037857
beta = 0.15184471
gamma = 0.35607169
delta = 1.81491784
'''


def xi_implicit(xi, kappa, chi1, eta):
    xi_implicit = 3 * (xi**2 - 2 * kappa * xi + kappa**2 * chi1**2)
    xi_implicit /= (xi**2 - 3 * kappa * xi + 2 * chi1 * np.sqrt(kappa**3 * xi))
    xi_implicit -= eta * xi**3
    return xi_implicit


def computeEjectae(m1, m2, chi1, chi2, eosname='2H', mode="computeDiskMass",
                   kerr=False, R_ns=None, max_mass=None):
    """
    This function computes the remnant disk mass after the coalescence using
    the equation (4) arXiv 1807.00011.

    and

    this function computes the dynamical ejecta after the coalescence of NS-BH
    using the equation (9) arXiv 2002.07728

    Parameters
    ----------
    m1, m2 : array_like
        primary and secondary mass(es)
    chi1, chi2 : array_like
        primary and secondary spin(s)
    eosname : str
        Name of the equation of state to be used. `AP4`
        `None` when supplying `R_ns`.
    mode : str
        Selection of the type of ejectae computed: computeDiskMass (O1-O4a
        calculation for HasRemnant), computeDynaEjecta (only for NSBH), both
    kerr : bool
        Supply to use the relativistic tidal parameter. See Fishbone (1971).
    R_ns : float
        Radius of the secondary in m, assuming it is a neutron star.
    max_mass : float
        Maximum mass of a neutron star. To be supplied if not supplying EoS.

    Example
    -------
    >>> computeEjectae(5.0, 2.0, 0.99, 0.)
    0.6321412881595185
    >>> m1 = np.array([5., 6., 7.])
    >>> m2 = np.array([1.0, 1.2, 1.6])
    >>> chi1 = np.zeros(3)
    >>> chi2 = np.zeros(3)
    >>> computeEjectae(m1, m2, chi1, chi2)
    array([0.12833991, 0.05054819, 0.])
    >>> computeEjectae(m1, m2, chi1, chi2, eosname='AP4')
    array([0.00851525, 0., 0.])
    >>> max_mass = 3.0  # in m_sun
    >>> r_ns = 15000.  # in meters
    >>> computeEjectae(m1, m2, chi1, chi2,
    ...                 R_ns=r_ns, max_mass=max_mass)
    array([0.12833991 0.05054819 0.])
    >>> # m1=2.0, m2=2.0 is a BNS event assuming 2H EOS
    >>> # DiskMass == 1.0 is a ad hoc value assigned for BNS events
    >>> masses_spins = np.array([2.0, 2.0, 0., 0.])
    >>> computeEjectae(masses_spins[0], masses_spins[1],
    ...                 masses_spins[2], masses_spins[3], eosname="2H")
    1.0
    >>> masses_spins = np.array([5.0, 2.0, 0.99, 0.])
    >>> computeEjectae(masses_spins[0], masses_spins[1],
    ...                 masses_spins[2], masses_spins[3], eosname="2H")
    0.6321412881595185
    >>> computeEjectae(masses_spins[0], masses_spins[1],masses_spins[2],
    ...                 masses_spins[3], mode="computeDiskMass",eosname="2H")
    0.6321412881595185
    >>> computeEjectae(masses_spins[0], masses_spins[1],masses_spins[2],
    ...                 masses_spins[3], mode="computeDynaEjecta",eosname="2H")
    0.05321543216488239
    >>> computeEjectae(masses_spins[0], masses_spins[1],masses_spins[2],
    ...                 masses_spins[3], mode="both",eosname="2H")
    0.14803662538881018

    Notes
    -----
    The primary mass, by convention, is larger one. If arrays are supplied,
    `m1`, `m2`, `chi1`, `chi1` should be of the same size.
    """
    try:
        _ = iter(m1)
        assert m1.size == m2.size == chi1.size == chi2.size, "masses and spins must have the same dimensions"
        # Making sure that the supplied m1 >= m2 ##
        largerMass1 = m1 >= m2
        largerMass2 = m2 > m1
        m_primary = np.zeros_like(m1)
        m_secondary = np.zeros_like(m2)
        chi_primary = np.zeros_like(chi1)
        chi_secondary = np.zeros_like(chi2)

        m_primary[largerMass1] = m1[largerMass1]
        m_primary[~largerMass1] = m2[largerMass2]
        m_secondary[~largerMass2] = m2[~largerMass2]
        m_secondary[largerMass2] = m1[~largerMass1]

        chi_primary[largerMass1] = chi1[largerMass1]
        chi_primary[~largerMass1] = chi2[~largerMass1]
        chi_secondary[~largerMass2] = chi2[~largerMass2]
        chi_secondary[largerMass2] = chi1[~largerMass1]
    except TypeError:
        if m1 >= m2:
            m_primary = m1
            m_secondary = m2
            chi_primary = chi1
            chi_secondary = chi2
        else:
            m_primary = m2
            m_secondary = m1
            chi_primary = chi2
            chi_secondary = chi1

    [C_ns, m2_b, max_mass] = computeCompactness(m_secondary, eosname=eosname,
                                                max_mass=max_mass)
    #print("C_ns",np.median(C_ns), "m2_b",np.median(m2_b), "max_mass",max_mass)

    eta = m_primary * m_secondary / (m_primary + m_secondary)**2
    BBH = m_secondary > max_mass
    BNS = m_primary < max_mass
    if not isinstance(BNS, np.ndarray):
        if BNS or BBH:
            return float(not BBH or BNS)

    if mode == "computeDiskMass" or mode == "both":
        if not kerr:
            alpha = 0.406
            beta = 0.139
            gamma = 0.255
            delta = 1.761
            firstTerm = alpha * (1 - 2 * C_ns) / (eta**(1 / 3))
            secondTerm = beta * compute_isco(chi1) * C_ns / eta
            thirdTerm = gamma
        else:
            alpha = 0.24037857
            beta = 0.15184471
            gamma = 0.35607169
            delta = 1.81491784
            kappa = (m_primary / m_secondary) * C_ns
            try:
                Num = len(kappa)
                xi = fsolve(xi_implicit, 100 * np.ones(Num),
                            args=(kappa, chi1, eta))
            except TypeError:
                xi = fsolve(xi_implicit, 100, args=(kappa, chi1, eta))

            firstTerm = alpha * xi * (1 - 2 * C_ns)
            secondTerm = beta * compute_isco(chi1) * C_ns / eta
            thirdTerm = gamma

        combinedTerms = firstTerm - secondTerm + thirdTerm

        if hasattr(combinedTerms, 'ndim') and combinedTerms.ndim >= 1:
            lessthanzero = combinedTerms < 0.0
            combinedTerms[lessthanzero] = 0.0
        else:
            combinedTerms = 0 if combinedTerms < 0 else combinedTerms
        M_rem_disk = (combinedTerms**delta) * m2_b

    if mode == "computeDynaEjecta" or mode == "both":
        a1 = 7.11595154E-03
        a2 = 1.43636803E-03
        a4 = -2.76202990E-02
        n1 = 8.63604211E-01
        n2= 1.68399507


        Q = m1 / m2
        #chi_eff=spins_to_chi_eff(m1, m2, chi1, chi2)
        firstTerm = a1 * np.power(Q, n1) * (1.0 - 2.0 * C_ns) / C_ns
        secondTerm = -a2*np.power(Q,n2)*compute_isco(chi1) #normalized by the m1 BH
        thirdTerm = a4
        #print(firstTerm* m2_b,secondTerm* m2_b,thirdTerm * m2_b)
        combinedTerms = firstTerm + secondTerm + thirdTerm

        if hasattr(combinedTerms, 'ndim') and combinedTerms.ndim >= 1:
            lessthanzero = combinedTerms < 0.0
            combinedTerms[lessthanzero] = 0.0
        else:
            combinedTerms = 0 if combinedTerms < 0 else combinedTerms

        M_rem_dyn = combinedTerms * m2_b

    # Computation of the total ejecta based on the sum of contributions
    # from the disk wind and dynamical ejecta component with an
    # arbitrary zeta = 0.15 (Fernandez (2014), Siegel and Metzger
    # 2018

    if mode == "both":
    	zeta_upper, zeta_upper=zeta_diskejecta(Q)
    	#print(np.round(zeta_lower,2))
    	M_rem = M_rem_dyn + zeta_upper * M_rem_disk

    if mode == "computeDynaEjecta":
        M_rem = M_rem_dyn

    if mode == "computeDiskMass":
    	Q = m1 / m2
    	zeta_lower, zeta_upper=zeta_diskejecta(Q)
    	M_rem = zeta_lower*M_rem_disk

    """
    if isinstance(BNS, np.ndarray):
        M_rem[BBH] = 0.0  # BBH is always EM-Dark
        M_rem[BNS] = 1.0  # Ad hoc value, assuming BNS is always EM-Bright
    """

    return M_rem
