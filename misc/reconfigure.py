#!/usr/bin/env python
# $Id$

import sys,os,string

# read the number of sites from mp-info on the wavefunction
def ReadNumSites(File):
    Command = "mp-info "+File+" | grep 'Number' | sed 's/Number of sites = //'"
    return int(os.popen(Command, 'r').readline()[:-1])

def ReadAttribute(File, Attr):
    return os.popen("mp-attr "+File+" "+Attr, 'r').readline()[:-1]

def SetAttribute(File, Attr, Value):
    os.system("mp-attr "+File+" "+Attr+"=\""+str(Value)+"\"")

def Energy(File, Lattice, Ham, LocalOp, i):
    return float(os.popen("mp-expectation -r "+File+" "+Lattice+":\"adjoint("
                          +LocalOp+"("+str(i)+"))*("+Ham
                          +")*"+LocalOp+"("+str(i)+")\"").readline()[:-1])

def Norm(File, Lattice, LocalOp, i):
    return float(os.popen("mp-expectation -r "+File+" "+Lattice+":\"adjoint("+LocalOp+"("+str(i)+"))*"
                          +LocalOp+"("+str(i)+")\"").readline()[:-1])

def main():
    argv = sys.argv[1:]

    if len(argv) != 3:
        sys.exit("reconfigure.py: error: expected parameters <in-psi> <local-operator> <out-psi>")

    Wavefunc = argv[0]
    LocalOp = argv[1]
    NewWavefunc = argv[2]

    NumSites = ReadNumSites(Wavefunc)/2
    Ham = ReadAttribute(Wavefunc, "Hamiltonian")
    
    [Lattice,Op] = string.split(Ham, ':', 2)

    MinSite = 0
    MinE = 1000
    for i in range(1,NumSites+1):
        N = Norm(Wavefunc, Lattice, LocalOp, i)
        if (N > 1e-16):
            E = Energy(Wavefunc, Lattice, Op, LocalOp, i) / N
            print "site="+str(i)+" energy="+str(E)
            if (E < MinE):
                MinSite = i
                MinE = E

    if MinSite == 0:
       sys.exit("reconfigure.py: unexpected: could not find a site to modify!")

    print "Minimum energy is "+str(MinE)+" at site "+str(MinSite)

    #NewWavefunc = Wavefunc+".r"
    os.system("mp-apply "+Lattice+":\""+LocalOp+"("+str(MinSite)+")\" "+Wavefunc+" "+NewWavefunc)
    os.system("mp-normalize "+NewWavefunc)
    SetAttribute(NewWavefunc, "Hamiltonian", Lattice+":"+Op)
    SetAttribute(NewWavefunc, "Energy", MinE)
    print "New wavefunction saved as "+NewWavefunc
        
if __name__ == '__main__': main()
