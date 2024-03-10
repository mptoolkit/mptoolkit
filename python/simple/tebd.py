#!/usr/bin/python3
# Matrix Product Toolkit http://mptoolkit.qusim.net/
#
# python/simple/tebd.py
#
# Copyright (C) 2020 Ian McCulloch <ian@qusim.net>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Research publications making use of this software should include
# appropriate citations and acknowledgements as described in
# the file CITATIONS in the main source directory.
#----------------------------------------------------------------------------
# ENDHEADER

#####################################################
# Simple TEBD code using MPS/MPO representations    #
# Ian McCulloch 25/5/2020                           #
#####################################################

import numpy as np
import scipy
import scipy.sparse.linalg
import scipy.sparse as sparse
import math
import cmath
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import sys
import copy

np.set_printoptions(threshold=sys.maxsize)

# pseudo-inverse of a diagonal matrix.
# We want to cut off values smaller than the threshold. We use
# Tikhonov regularization, s -> s / (s^2 + a^2), where a is the cutoff
def pseudo_inverse(S):
    a2 = 1e-7 ** 2
    return np.array([x / (x**2 + a2) for x in S])

def norm_frob(A):
    return np.sum(np.multiply(np.conj(A),A))

## take an MPS in right orthogonal form and construct the canonical form
## Lambda_0 Gamma_0 Lambda_1 Gamma_1 ... Gamma_N Lambda_N+1
def canonicalize(MPS):
    Lambda = [np.array([1])]
    Gamma = []
    V = np.array([[1]])
    lastS = np.array([1])
    for i in range(len(MPS)):
        # reshape into m x (dm)
        Temp = np.einsum("ij,sjk->sik", V, MPS[i])
        s = np.shape(Temp)
        A = np.reshape(Temp, (s[0]*s[1], s[2]))
        U, S, V = np.linalg.svd(A, full_matrices=0)
        Gamma = Gamma + [np.einsum("i,sij->sij", pseudo_inverse(lastS), np.reshape(U, (s[0], s[1], np.shape(U)[1])))]
        Lambda = Lambda + [S]
        V = np.einsum("i,ij->ij", S, V)
        lastS = S
    # assert V is 1x1 identity matrix
    return Lambda, Gamma

def orthonormalize(MPS):
    # sweep left-to-right with SVD's
    V = np.array([[1]])
    for i in range(len(MPS)):
        s = np.shape(MPS[i])
        A = np.reshape(np.einsum("ij,sjk->sik",V,MPS[i]), (s[0]*s[1], s[2]))
        U,S,V = np.linalg.svd(A, full_matrices=0)
        MPS[i] = np.reshape(U, (s[0], s[1], np.shape(U)[1]))
        V = np.einsum("i,ij->ij", S, V)
    # sweep right to left
    U = np.array([[1]])
    for i in range(len(MPS)-1, -1, -1):
        s = np.shape(MPS[i])
        A = np.reshape(np.einsum("sij,jk->isk", MPS[i], U), (s[1], s[0]*np.shape(U)[1]))
        U,S,V = np.linalg.svd(A, full_matrices=0)
        MPS[i] = np.einsum("isk->sik", np.reshape(V, (np.shape(V)[0], s[0], s[2])))
        U = np.einsum("ij,j->ij", U, S)
    # the final U is the norm, which we discard
    return MPS    

def TEBD_Gate(Lambda0, Gamma1, Lambda1, Gamma2, Lambda2, Gate, m):
    # contract the outside lambda's into Gamma
    Gamma1 = np.einsum("i,sij->sij", Lambda0, Gamma1)
    Gamma2 = np.einsum("sij,j->sij", Gamma2, Lambda2)
    # apply the centre lambda
    Gamma2 = np.einsum("i,sij->sij", Lambda1, Gamma2)
    # make 2-body tensor Gamma1 Lambda1 Gamma2
    Temp = np.einsum("sij,tjk->stik", Gamma1, Gamma2)
    # apply the 2-body gate and reorder to apply the svd
    TShape = np.shape(Temp)
    GShape = (TShape[0],TShape[1],TShape[0],TShape[1])
    Temp = np.einsum("uvst,stik->uivk", np.reshape(Gate,GShape), Temp)
    # SVD
    s = np.shape(Temp)
    A = np.reshape(Temp, (s[0]*s[1], s[2]*s[3]))
    U, S, V = np.linalg.svd(A, full_matrices=0)
    # truncate
    m = min(len(S), m)
    S = S[0:m]
    U = U[:,0:m]
    V = V[0:m,:]
    # reshape back to Gamma
    Gamma1 = np.reshape(U, (s[0], s[1], np.shape(U)[1]))
    Gamma2 = np.einsum("isj->sij", np.reshape(V, (np.shape(V)[0], s[2], s[3])))
    # invert the outer Lambda matrices
    Gamma1 = np.einsum("i,sij->sij", pseudo_inverse(Lambda0), Gamma1)
    Gamma2 = np.einsum("sij,j->sij", Gamma2, pseudo_inverse(Lambda2))
    return Gamma1, S, Gamma2

def EvolveEvenBonds(Lambda, Gamma, Gate, m):
    for i in range(0, len(Gamma)-1, 2):
        Gamma[i], Lambda[i+1], Gamma[i+1] = TEBD_Gate(Lambda[i], Gamma[i], Lambda[i+1], Gamma[i+1], Lambda[i+2], Gate, m)
    return Lambda, Gamma

def EvolveOddBonds(Lambda, Gamma, Gate, m):
    for i in range(1, len(Gamma)-1, 2):
        Gamma[i], Lambda[i+1], Gamma[i+1] = TEBD_Gate(Lambda[i], Gamma[i], Lambda[i+1], Gamma[i+1], Lambda[i+2], Gate, m)
    return Lambda, Gamma

def MPS_from_canonical(Lambda, Gamma):
    MPS = []
    for i in range(len(Gamma)):
        MPS = MPS + [np.einsum("i,sij->sij", Lambda[i], Gamma[i])]
    return MPS

# evolve a timestep, using simple 2nd order Lie-Suzuki-Trotter decomposition
def EvolveTimestep(Lambda, Gamma, H, m, dt):
    Lambda,Gamma = EvolveEvenBonds(Lambda, Gamma, scipy.linalg.expm(complex(0,-dt/2)*H), m)
    Lambda,Gamma = EvolveOddBonds(Lambda, Gamma, scipy.linalg.expm(complex(0,-dt)*H), m)
    Lambda,Gamma = EvolveEvenBonds(Lambda, Gamma, scipy.linalg.expm(complex(0,-dt/2)*H), m)
    return Lambda,Gamma
    
def LocalExpectationValue(Lambda, Gamma, Site, Op):
    G = np.einsum("i,sij->sij", Lambda[Site], Gamma[Site])
    G = np.einsum("sij,j->sij", G, Lambda[Site+1])
    H = np.einsum("st,tij->sij", Op, G)
    return np.einsum("sij,sij", np.conj(G), H)

# return an array of the local expectation value at every site
def Density(Lambda, Gamma, Op):
    return np.array([LocalExpectationValue(Lambda, Gamma, n, Op) for n in range(0,len(Gamma))])

## initial E and F matrices for the left and right vacuum states
def initial_E(W):
    E = np.zeros((W.shape[0],1,1))
    E[0] = 1
    return E

def initial_F(W):
    F = np.zeros((W.shape[1],1,1))
    F[-1] = 1
    return F

def contract_from_right(W, A, F, B):
    # the einsum function doesn't appear to optimize the contractions properly,
    # so we split it into individual summations in the optimal order
    #return np.einsum("abst,sij,bjl,tkl->aik",W,A,F,B, optimize=True)
    Temp = np.einsum("sij,bjl->sbil", np.conj(A), F)
    Temp = np.einsum("sbil,abst->tail", Temp, W)
    return np.einsum("tail,tkl->aik", Temp, B)

def contract_from_left(W, A, E, B):
    # the einsum function doesn't appear to optimize the contractions properly,
    # so we split it into individual summations in the optimal order
    # return np.einsum("abst,sij,aik,tkl->bjl",W,A,E,B, optimize=True)
    Temp = np.einsum("sij,aik->sajk", np.conj(A), E)
    Temp = np.einsum("sajk,abst->tbjk", Temp, W)
    return np.einsum("tbjk,tkl->bjl", Temp, B)

# construct the F-matrices for all sites except the first
def construct_F(Alist, MPO, Blist):
    F = [initial_F(MPO[-1])]

    for i in range(len(MPO)-1, 0, -1):
        F.append(contract_from_right(MPO[i], Alist[i], F[-1], Blist[i]))
    return F

def construct_E(Alist, MPO, Blist):
    return [initial_E(MPO[0])]



# 2-1 coarse-graining of two site MPO into one site
def coarse_grain_MPO(W, X):
    return np.reshape(np.einsum("abst,bcuv->acsutv",W,X),
                      [W.shape[0], X.shape[1],
                       W.shape[2]*X.shape[2],
                       W.shape[3]*X.shape[3]])

def product_W(W1, W2):
    return np.reshape(np.einsum("abst,cdtu->acbdsu", W1, W2), [W1.shape[0]*W2.shape[0],
                                                               W1.shape[1]*W2.shape[1],
                                                               W1.shape[2],W2.shape[3]])

def product_MPO(M1, M2):
    assert len(M1) == len(M2)
    Result = []
    for i in range(0, len(M1)):
        Result.append(product_W(M1[i], M2[i]))
    return Result


# 2-1 coarse-graining of two-site MPS into one site
def coarse_grain_MPS(A,B):
    return np.reshape(np.einsum("sij,tjk->stik",A,B),
                      [A.shape[0]*B.shape[0], A.shape[1], B.shape[2]])

def fine_grain_MPS(A, dims):
    assert A.shape[0] == dims[0] * dims[1]
    Theta = np.transpose(np.reshape(A, dims + [A.shape[1], A.shape[2]]),
                         (0,2,1,3))
    M = np.reshape(Theta, (dims[0]*A.shape[1], dims[1]*A.shape[2]))
    U, S, V = np.linalg.svd(M, full_matrices=0)
    U = np.reshape(U, (dims[0], A.shape[1], -1))
    V = np.transpose(np.reshape(V, (-1, dims[1], A.shape[2])), (1,0,2))
    # assert U is left-orthogonal
    # assert V is right-orthogonal
    #print(np.dot(V[0],np.transpose(V[0])) + np.dot(V[1],np.transpose(V[1])))
    return U, S, V

def truncate_MPS(U, S, V, m):
    m = min(len(S), m)
    trunc = np.sum(S[m:])
    S = S[0:m]
    U = U[:,:,0:m]
    V = V[:,0:m,:]
    return U,S,V,trunc,m

def Expectation(AList, MPO, BList):
    E = [[[1]]]
    for i in range(0,len(MPO)):
        E = contract_from_left(MPO[i], AList[i], E, BList[i])
    return E[0][0][0]

class HamiltonianMultiply(sparse.linalg.LinearOperator):
    def __init__(self, E, W, F):
        self.E = E
        self.W = W
        self.F = F
        self.dtype = np.dtype('d')
        self.req_shape = [W.shape[2], E.shape[1], F.shape[2]]
        self.size = self.req_shape[0]*self.req_shape[1]*self.req_shape[2]
        self.shape = [self.size, self.size]

    def _matvec(self, A):
        # the einsum function doesn't appear to optimize the contractions properly,
        # so we split it into individual summations in the optimal order
        #R = np.einsum("aij,sik,abst,bkl->tjl",self.E,np.reshape(A, self.req_shape),
        #              self.W,self.F, optimize=True)
        R = np.einsum("aij,sik->ajsk", self.E, np.reshape(A, self.req_shape))
        R = np.einsum("ajsk,abst->bjtk", R, self.W)
        R = np.einsum("bjtk,bkl->tjl", R, self.F)
        return np.reshape(R, -1)

## optimize a single site given the MPO matrix W, and tensors E,F
def optimize_site(A, W, E, F):
    H = HamiltonianMultiply(E,W,F)
    # we choose tol=1E-8 here, which is OK for small calculations.
    # to bemore robust, we should take the tol -> 0 towards the end
    # of the calculation.
    E,V = sparse.linalg.eigsh(H,1,v0=A,which='SA', tol=1E-8)
    return (E[0],np.reshape(V[:,0], H.req_shape))

def optimize_two_sites(A, B, W1, W2, E, F, m, dir):
    W = coarse_grain_MPO(W1,W2)
    AA = coarse_grain_MPS(A,B)
    H = HamiltonianMultiply(E,W,F)
    E,V = sparse.linalg.eigsh(H,1,v0=AA,which='SA')
    AA = np.reshape(V[:,0], H.req_shape)
    A,S,B = fine_grain_MPS(AA, [A.shape[0], B.shape[0]])
    A,S,B,trunc,m = truncate_MPS(A,S,B,m)
    if (dir == 'right'):
        B = np.einsum("ij,sjk->sik", np.diag(S), B)
    else:
        assert dir == 'left'
        A = np.einsum("sij,jk->sik", A, np.diag(S))
    return E[0], A, B, trunc, m

def two_site_dmrg(MPS, MPO, m, sweeps):
    E = construct_E(MPS, MPO, MPS)
    F = construct_F(MPS, MPO, MPS)
    F.pop()
    for sweep in range(0,int(sweeps/2)):
        for i in range(0, len(MPS)-2):
            Energy,MPS[i],MPS[i+1],trunc,states = optimize_two_sites(MPS[i],MPS[i+1],
                                                                     MPO[i],MPO[i+1],
                                                                     E[-1], F[-1], m, 'right')
            print("Sweep {:} Sites {:},{:}    Energy {:16.12f}    States {:4} Truncation {:16.12f}"
                     .format(sweep*2,i,i+1, Energy, states, trunc))
            E.append(contract_from_left(MPO[i], MPS[i], E[-1], MPS[i]))
            F.pop();
        for i in range(len(MPS)-2, 0, -1):
            Energy,MPS[i],MPS[i+1],trunc,states = optimize_two_sites(MPS[i],MPS[i+1],
                                                                     MPO[i],MPO[i+1],
                                                                     E[-1], F[-1], m, 'left')
            print("Sweep {} Sites {},{}    Energy {:16.12f}    States {:4} Truncation {:16.12f}"
                     .format(sweep*2+1,i,i+1, Energy, states, trunc))
            F.append(contract_from_right(MPO[i+1], MPS[i+1], F[-1], MPS[i+1]))
            E.pop();
    return MPS


            

d=2   # local bond dimension
N=40 # number of sites

InitialA1 = np.zeros((d,1,1))
InitialA1[0,0,0] = 1
InitialA2 = np.zeros((d,1,1))
InitialA2[1,0,0] = 1

## initial state |01010101>
MPS = [InitialA1, InitialA2] * int(N/2)

## Local operators
I = np.identity(2)
Z = np.zeros((2,2))
Sz = np.array([[0.5,  0  ],
             [0  , -0.5]])
Sp = np.array([[0, 1],
             [0, 0]])
Sm = np.array([[0, 0],
             [1, 0]])

## Hamiltonian MPO
W = np.array([[I, Sz, 0.5*Sp,  0.5*Sm,  Z],
              [Z,  Z,      Z,       Z, Sz], 
              [Z,  Z,      Z,       Z, Sm],
              [Z,  Z,      Z,       Z, Sp],
              [Z,  Z,      Z,       Z,  I]])

Wfirst = np.array([[I, Sz, 0.5*Sp, 0.5*Sm,   Z]])

Wlast = np.array([[Z], [Sz], [Sm], [Sp], [I]])

# the complete MPO
MPO = [Wfirst] + ([W] * (N-2)) + [Wlast]

HamSquared = product_MPO(MPO, MPO)

m = 10
sweeps = 8
MPS = two_site_dmrg(MPS, MPO, m, sweeps)

Energy = Expectation(MPS, MPO, MPS)
print("Final energy expectation value {}".format(Energy))

H2 = Expectation(MPS, HamSquared, MPS)
print("variance = {:16.12f}".format(H2 - Energy*Energy))

# sage the groundstate, we need it later
GS = copy.deepcopy(MPS)

# Act on the MPS with a local operator at some site
MPS[N//2] = np.einsum("st,tij->sij", Sp, MPS[N//2])

# orthonormalize the MPS
MPS = orthonormalize(MPS)

# canonicalize in preparation for TEND
Lambda,Gamma = canonicalize(MPS)

# 2 sites of the Hamiltonian.  Full Hamiltonian is sum over translations of H2site
H2site = 0.5 * (np.kron(Sp, Sm) + np.kron(Sm, Sp)) + np.kron(Sz, Sz)

# matrix of Sz values as a function of time and position.  Rows are timesteps, columns are position
SzMatrix = [Density(Lambda, Gamma, Sz)]

# for the correlation function, we want f(x,t) = <0|S-(t) S+(0)|0> = <0|S- e^{iHt} S+(0)|0> (ignoring phase factor)
IdentityMPO = [np.array([[I]]) for i in range(N)]
Correlator = []
print(I)
for i in range(N):
    c = copy.deepcopy(IdentityMPO)
    c[i][0,0] = Sm
    Correlator = Correlator + [c]

# initial correlation
Correlation = [ [Expectation(GS, Correlator[n], MPS) for n in range(N)] ]

# timestep
dt = 0.02
NT = 600 # number of timesteps
Skip = 1 # only calculate correlations after this many timesteps
for i in range(NT):
    print("timestep {} / {}".format(i,NT))
    Lambda,Gamma = EvolveTimestep(Lambda, Gamma, H2site, 20, dt)
    SzMatrix = SzMatrix + [Density(Lambda, Gamma, Sz)]
    # correlation function
    if (i % Skip == 0):
        MPS = MPS_from_canonical(Lambda, Gamma)
        Correlation = Correlation + [ [Expectation(GS, Correlator[n], MPS) for n in range(N)] ]
    

X = np.outer(np.linspace(0,N-1,N), np.ones(NT+1))
T = np.outer(np.ones(N), np.linspace(0,dt,NT+1))

SzMatrix = np.transpose(np.array(SzMatrix).real)

#np.set_printoptions(threshold=sys.maxsize)
#print(SzMatrix)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.plot_surface(X, T, SzMatrix, cmap='viridis', edgecolor='none')
ax.set_xlabel("X")
ax.set_ylabel("Time")
ax.set_zlabel("Sz")
ax.set_title('Evolution of the magnetization')
fig.show()

# correlation function
Correlation = np.transpose(np.array(Correlation))

figc = plt.figure()
axc = plt.axes(projection='3d')
axc.plot_surface(X, T, np.absolute(Correlation), cmap='viridis', edgecolor='none')
axc.set_xlabel("X")
axc.set_ylabel("Time")
axc.set_zlabel("Sz")
axc.set_title('Time-dependent correlation')
figc.show()

# Calculate the spectral function
# time evolution is symmetric, so we can add the -ve time component
Correlation = np.concatenate((np.flip(Correlation[:,1:], axis=1), Correlation), axis=1)

# FFT
Spec = np.absolute(np.fft.fft2(Correlation))

# we didn't correct for the groundstate energy, so just plot by hand the relevant part
Spec = Spec[:,68:52:-1]
NumOmega=np.shape(Spec)[1]

K = np.outer(np.linspace(0,np.pi*2,N), np.ones(NumOmega)) 
W = np.outer(np.ones(N), np.linspace(0,NumOmega*np.pi/(dt*NT), NumOmega))

print(np.shape(Spec))
print(np.shape(K))
print(np.shape(W))

fig2 = plt.figure()
ax2 = plt.axes(projection='3d')
ax2.plot_surface(K, W, np.absolute(Spec), cmap='viridis', edgecolor='none')
ax2.set_xlabel("k")
ax2.set_ylabel("Energy")
ax2.set_zlabel("Spectral density")
ax2.set_title('Spectral function')
fig2.show()
input()
