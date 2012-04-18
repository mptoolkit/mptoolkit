// -*- C++ -*- $Id$
//
// functions for operators acting on wavefunctions

#if !defined(OPERATOR_ACTIONS_H_SDH47589FOIJHO9JEW)
#define OPERATOR_ACTIONS_H_SDH47589FOIJHO9JEW

#include "mps/linearwavefunction.h"
#include "mpo/linear_operator.h"


MatrixOperator transfer_from_left(MatrixOperator const& m, 
                                  LinearOperator const& Op, 
                                  LinearWavefunction const& Psi);

MatrixOperator transfer_from_right(MatrixOperator const& m, 
                                   LinearOperator const& Op, 
                                   LinearWavefunction const& Psi);

#endif
