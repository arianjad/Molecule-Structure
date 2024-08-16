import numpy as np
from sympy.physics.wigner import wigner_3j,wigner_6j,wigner_9j

# Kronecker delta function
def kronecker(a,b):
    if a==b:
        return 1
    else:
        return 0

# Rotation
def Rot_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(K0,K1)*kronecker(N0,N1)*kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1):
        return 0
    else:
        return N0*(N0+1)-K0**2

# Spin Rotation
def SR_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(K0,K1)*kronecker(J0,J1)*kronecker(N0,N1)*kronecker(F0,F1)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(N0+J0+S)*np.sqrt(S*(S+1)*(2*S+1)*N0*(N0+1)*(2*N0+1))*\
            wigner_6j(S,N0,J0,N0,S,1)

# bF Fermi Contact
def IS_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(K0,K1)*kronecker(N0,N1)*kronecker(F0,F1)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(J1 + F0 + I + J0 + N0 + S + 1)*\
            np.sqrt((2*J1+1)*(2*J0+1)*S*(S+1)*(2*S+1)*I*(I+1)*(2*I+1))*\
            wigner_6j(I,J1,F0,J0,I,1)*wigner_6j(J0,S,N0,S,J1,1)

# c Dipole dipole
def T2IS_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(K0,K1):
        return 0
    else:
        return -np.sqrt(5/3)*(-1)**(3*F0-2*M0+I+J0+N0-K0)*np.sqrt((2*I+1)*(I+1)*I)*\
            np.sqrt((2*J0+1)*(2*J1+1)*3*(2*S+1)*(S+1)*S)*wigner_6j(I,J1,F0,J0,I,1)*\
            wigner_9j(S,N1,J1,1,2,1,S,N0,J0)*np.sqrt((2*N0+1)*(2*N1+1))*wigner_3j(N0,2,N1,-K0,0,K1)

# From when I used to calculate IzSz
def Iz_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1):
        return 0
    else:
        return (-1)**(I+F0+S-K0+1)*np.sqrt((2*J1+1)*(2*J0+1)*(2*N0+1)*(2*N1+1)*I*(I+1)*(2*I+1))*\
            wigner_6j(J1,I,F0,I,J0,1)*wigner_6j(N1,J1,S,J0,N0,1)*wigner_3j(N0,1,N1,-K0,0,K1)

# For EDM
def Sz_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1):
        return 0
    else:
        return (-1)**(N1+N0+S+J0-K0)*np.sqrt((2*N0+1)*(2*N1+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(N1,S,J0,S,N0,1)*wigner_3j(N0,1,N1,-K0,0,K1)

# For Spin Rotation in bending mode
def NzSz_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(F0,F1)*kronecker(M0,M1)*kronecker(J0,J1)*kronecker(K0,K1):
        return 0
    else:
        return K0*(-1)**(N1+N0+S+J0-K0)*np.sqrt((2*N0+1)*(2*N1+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(N1,S,J0,S,N0,1)*wigner_3j(N0,1,N1,-K0,0,K1)

# Zeeman Bz interaction, just from electron spin
def ZeemanZ_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(K0,K1)*kronecker(M0,M1)*kronecker(N0,N1):
        return 0
    else:
        return (-1)**(F0-M0+F1+2*J0+I+N0+S)*np.sqrt((2*F0+1)*(2*F1+1)*(2*J0+1)*(2*J1+1)*S*(S+1)*(2*S+1))*\
            wigner_6j(F0,J0,I,J1,F1,1)*wigner_6j(J0,S,N0,S,J1,1)*wigner_3j(F0,1,F1,-M0,0,M1)

# Stark Ez interaction with molecule frame dipole
def StarkZ_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(M0,M1):
        return 0
    else:
        return (-1)**(F0-M0+J0+I+2+F1+S+J1-K0)*np.sqrt((2*F0+1)*(2*F1+1)*(2*J0+1)*(2*J1+1)*(2*N0+1)*(2*N1+1))*\
            wigner_3j(F0,1,F1,-M0,0,M1)*wigner_6j(J1,F1,I,F0,J0,1)*wigner_6j(N1,J1,S,J0,N0,1)*wigner_3j(N0,1,N1,-K0,0,K1)

# l-doubling p term
def p_lD_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(M0,M1)*kronecker(F0,F1)*kronecker(J0,J1)*(not kronecker(K0,K1)):
        return 0
    else:
        return sum([(-1)**(-2*M0+2*I+3*J0+S+N0)*-np.sqrt(10/6)*wigner_6j(S,N0,J0,N1,S,1)*np.sqrt(S*(S+1)*(2*S+1))*\
            np.sqrt(3)*(-1)**(1+N0+N1)*wigner_6j(2,1,1,N1,N0,N0)*np.sqrt(N0*(N0+1)*(2*N0+1))*\
            (-1)**(N0-K0)*wigner_3j(N0,2,N1,-K0,2*q,K1)*np.sqrt((2*N0+1)*(2*N1+1))*kronecker(K0,K1+2*q) for q in [-1,1]])

# l-doubling q term
def q_lD_bBJ(K0,N0,J0,F0,M0,K1,N1,J1,F1,M1,S=1/2,I=1/2):
    if not kronecker(M0,M1)*kronecker(F0,F1)*kronecker(J0,J1)*(not kronecker(K0,K1)):
        return 0
    else:
        return sum([(-1)**(-2*M0+2*I+2*J0)*np.sqrt((2*J1+1))*(-1)**(S+J1+N1)*wigner_6j(N1,J1,S,J0,N0,0)*\
            np.sqrt(5)*(-1)**(N0+N1)*wigner_6j(2,2,0,N1,N0,N0)*(-1)**(N0+K0)*wigner_3j(N0,2,N1,-K0,2*q,K1)*\
            1/(2*np.sqrt(6))*np.sqrt((2*N0-1)*(2*N0)*(2*N0+1)*(2*N0+2)*(2*N0+3)*(2*N0+1)*(2*N1+1))*kronecker(K0,K1+2*q) for q in [-1,1]])
