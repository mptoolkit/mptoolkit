#!/usr/bin/env python
# $Id$

import sys,os,string

NumStates = 200

NumStatesStr = str(NumStates)

def Normalize(x):
    os.system("mp-normalize "+x)

def Apply(Op, x, res):
    os.system("cp "+x+" "+res)
    os.system("mp-apply-opt -m "+NumStatesStr+" --no-resid "+Op+" "+x+" "+res)
    os.system("mp-apply-opt -m "+NumStatesStr+" -f 0 -2 --no-resid "+Op+" "+x+" "+res)

def Apply2(Op1, x1, Op2, x2, res):
    os.system("cp "+x1+" "+res)
    os.system("mp-apply-multiple -m "+NumStatesStr+" -w "+res
              +" --operator "+Op1+" --rhs "+x1
              +" --operator "+Op2+" --rhs "+x2)
    os.system("mp-apply-multiple -m "+NumStatesStr+" -f 0 -2 -w "+res
              +" --operator "+Op1+" --rhs "+x1
              +" --operator "+Op2+" --rhs "+x2)

def Apply3(Op1, x1, Op2, x2, Op3, x3, res):
    os.system("cp "+x1+" "+res)
    os.system("mp-apply-multiple -m "+NumStatesStr+" -w "+res
              +" --operator "+Op1+" --rhs "+x1
              +" --operator "+Op2+" --rhs "+x2
              +" --operator "+Op3+" --rhs "+x3)
    os.system("mp-apply-multiple -m "+NumStatesStr+" -f 0 -2 -w "+res
              +" --operator "+Op1+" --rhs "+x1
              +" --operator "+Op2+" --rhs "+x2
              +" --operator "+Op3+" --rhs "+x3)

def Norm(x, y):
    return float(os.popen("mp-overlap -r "+x+" "+y).readline()[:-1])

def Expectation(File, Op, i):
    return float(os.popen("mp-expectation -r "+File+" "+Op+" "+i).readline()[:-1])

def main():
    argv = sys.argv[1:]

    if len(argv) != 4:
        sys.exit("reconfigure.py: error: expected parameters <in-psi> <H> <length> <out-prefix>")

    Wavefunc = argv[0]
    Ham = argv[1]
    Length = int(argv[2])
    OutPrefix = argv[3]

    [Lattice,Op] = string.split(Ham, ':', 2)

    # copy the initial wavefunction as krylov vector 0
    os.system("cp "+Wavefunc+" "+OutPrefix+".0")
    Normalize(OutPrefix+".0")

    # first step in the Krylov sequence
    a = Expectation(OutPrefix+".0", Ham, OutPrefix+".0")
    print "first tridiagonal coefficient a=",a
    Apply2(Ham, OutPrefix+".0",  Lattice+":"+str(-a), OutPrefix+".0", OutPrefix+".1")
    Normalize(OutPrefix+".1")

    for i in range(1,Length-1):
        a = Expectation(OutPrefix+"."+str(i), Ham, OutPrefix+"."+str(i))
        b2 = Expectation(OutPrefix+"."+str(i-1), Ham, OutPrefix+"."+str(i))
        print "tridiagonal coefficients a=",a," b^2=",b2
        Apply3(Ham, OutPrefix+"."+str(i), Lattice+":"+str(-a), OutPrefix+"."+str(i),
               Lattice+":"+str(-b2), OutPrefix+"."+str(i-1), OutPrefix+"."+str(i+1))
        Normalize(OutPrefix+"."+str(i+1))
    
if __name__ == '__main__': main()
