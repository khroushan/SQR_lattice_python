# Program to calculate the dc conductance across a square lattice.
# Authur: Amin Ahmdi
# Date: Nov 25, 2016
# ################################################################
import numpy as np
import numpy.linalg as lg

N_energy = 200
eta = 1.e-7j
mlat = 10
I = np.eye(mlat,dtype=complex)

def make_h_tau(mlat):
    """ Constructs the Hamiltonian and the connection 
    matrix of a square lattice
    """
    tau = -np.eye(mlat,dtype=complex)
    h = np.zeros((mlat,mlat), dtype=complex)
    
    for i in range(mlat-1):
        h[i,i+1] = -1.0
        h[i+1,i] = np.conjugate(h[i,i+1])

    return h, tau



def g_lead(E, tau, h):
    """ Compute the lead's Green's function using decimation
    method. Inputs:
    E: energy
    tau: connection matrix between slices
    h: The Hamiltonian of the slice
    Return:
    Lead's Green's function
    """
    eta = 1.e-7j               # tiny imaginary for retarded G.Fs.
    I = np.eye(mlat,dtype=complex)
    
    ee = E + eta
    # initialize alpha, beta, eps, eps_s
    alpha = tau
    beta = np.conjugate(tau.T)
    eps = h
    eps_s = h

    for i_dec in range(30):
        aux = lg.inv(ee*I - eps)
        aux1 = np.dot(alpha,np.dot(aux,beta))
        eps_s = eps_s + aux1
        eps = eps + aux1
        aux1 = np.dot(beta, np.dot(aux, alpha))
        eps = eps + aux1
        alpha = np.dot(alpha,np.dot(aux, alpha))
        beta = np.dot(beta, np.dot(aux, beta))


    g_l = lg.inv(ee*I - eps_s)
    return g_l

# construct the Hamiltonian
h, tau = make_h_tau(mlat)


# loop over energy 
for ie in range(N_energy):

    energy = -4.0 +  ie*(8.0/N_energy)

    # The lead's Green's function
    g_l = g_lead(energy, tau, h)
    
    # The self-energy due to the leads
    sigma_l = np.dot(tau, np.dot(g_l,tau))
    sigma_r = sigma_l                       #in case of SQR lattice
    sigma_l_dg = np.conjugate(sigma_l.T)

    # the Green's function of device in presence of leads' self-energies
    G_d = lg.inv(energy*I - h - 2*sigma_l)
    G_d_dg = np.conjugate(G_d.T)
    
    gamma_l = -1j * ( sigma_l - sigma_l_dg)
    gamma_r = gamma_l

    # G * gamma_l * G_dg * gamma_r
    auxg = np.dot(G_d,np.dot(gamma_l,np.dot(G_d_dg,gamma_r)))

    gg = np.trace(auxg)

    print energy, gg.real, gg.imag
