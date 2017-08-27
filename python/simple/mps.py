#####################################################
# Simple DMRG code using MPS/MPO representations    #
# Ian McCulloch August 2017                         #
#####################################################

import numpy as np
import scipy
import scipy.sparse.linalg
import scipy.sparse as sparse
import math

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
    return np.einsum("abst,sij,bjl,tkl->aik",W,A,F,B)

def contract_from_left(W, A, E, B):
    return np.einsum("abst,sij,aik,tkl->bjl",W,A,E,B)

# construct the E-matrices for all sites except the first
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
        R = np.einsum("aij,abst,bkl,sik->tjl",self.E,self.W,self.F,np.reshape(A, self.req_shape))
        return np.reshape(R, -1)

## optimize a single site given the MPO matrix W, and tensors E,F
def optimize_site(A, W, E, F):
    H = HamiltonianMultiply(E,W,F)
    E,V = sparse.linalg.eigsh(H,1,v0=A,which='SA')
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
    for sweep in range(0,sweeps/2):
        for i in range(0, len(MPS)-2):
            Energy,MPS[i],MPS[i+1],trunc,states = optimize_two_sites(MPS[i],MPS[i+1],
                                                                     MPO[i],MPO[i+1],
                                                                     E[-1], F[-1], m, 'right')
            print("Sweep {} Sites {},{}    Energy {:16.12f}    States {:4} Truncation {:16.12f}") \
                     .format(sweep*2,i,i+1, Energy, states, trunc)
            E.append(contract_from_left(MPO[i], MPS[i], E[-1], MPS[i]))
            F.pop();
        for i in range(len(MPS)-2, 0, -1):
            Energy,MPS[i],MPS[i+1],trunc,states = optimize_two_sites(MPS[i],MPS[i+1],
                                                                     MPO[i],MPO[i+1],
                                                                     E[-1], F[-1], m, 'left')
            print("Sweep {} Sites {},{}    Energy {:16.12f}    States {:4} Truncation {:16.12f}") \
                     .format(sweep*2+1,i,i+1, Energy, states, trunc)
            F.append(contract_from_right(MPO[i+1], MPS[i+1], F[-1], MPS[i+1]))
            E.pop();
    return MPS
            

d=2   # local bond dimension
N=100 # number of sites

InitialA1 = np.zeros((d,1,1))
InitialA1[0,0,0] = 1
InitialA2 = np.zeros((d,1,1))
InitialA2[1,0,0] = 1

## initial state |01010101>
MPS = [InitialA1, InitialA2] * (N/2)

## Local operators
I = np.identity(2)
Z = np.zeros((2,2))
Sz = np.array([[0.5,  0  ],
             [0  , -0.5]])
Sp = np.array([[0, 0],
             [1, 0]])
Sm = np.array([[0, 1],
             [0, 0]])

## Hamiltonian MPO
W = np.array([[I, Sz, 0.5*Sp, 0.5*Sm,   Z],
              [Z,  Z,      Z,       Z, Sz], 
              [Z,  Z,      Z,       Z, Sm],
              [Z,  Z,      Z,       Z, Sp],
              [Z,  Z,      Z,       Z,  I]])

Wfirst = np.array([[I, Sz, 0.5*Sp, 0.5*Sm,   Z]])

Wlast = np.array([[Z], [Sz], [Sm], [Sp], [I]])

# the complete MPO
MPO = [Wfirst] + ([W] * (N-2)) + [Wlast]

HamSquared = product_MPO(MPO, MPO)

MPS = two_site_dmrg(MPS, MPO, 10, 6)

Energy = Expectation(MPS, MPO, MPS)
print("Final energy expectation value {}").format(Energy)

H2 = Expectation(MPS, HamSquared, MPS)
print("variance = {:16.12f}").format(H2 - Energy*Energy)
