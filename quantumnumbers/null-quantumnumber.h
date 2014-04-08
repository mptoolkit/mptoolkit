/* -*- C++ -*- $Id$

  The trivial null quantum number.
*/

#if !defined(NULL_QUANTUMNUMBER_H_HG768HJ)
#define NULL_QUANTUMNUMBER_H_HG768HJ

#include "common/niftycounter.h"
#include "common/halfint.h"
#include "quantumnumber.h"

namespace QuantumNumbers
{

class NullQN
{
   public:
      typedef NullQN ProjectionType;
      typedef NullQN QuantumNumberType;
      typedef StaticQuantumNumberFactory<NullQN> FactoryType;

      NullQN() {}
      explicit NullQN(int const* InIter);
      explicit NullQN(std::string const& s);

      std::string ToString() const;

      int* Convert(int* OutIter) const
      { *OutIter = 0; return OutIter+1; }

      static char const* Type() { return "Null"; }

      static int Size() { return 1; }

      static int num_casimir() { return 0; }
      static std::string casimir_name(std::string const& QName, int) { return QName; }

      // Registation is automatic via a nifty counter
      static void Register();

      static char const* Suffix() { return ""; }

      bool operator==(NullQN const&) const { return true; }
      bool operator!=(NullQN const&) const { return false; }
};

// the empty symmetry list Null:Null
extern SymmetryList NullSymmetryList;

// the null quantum number
extern QuantumNumber Null;

inline
std::ostream& operator<<(std::ostream& out, NullQN const& s)
{
   return out;
}

typedef NullQN NullQNProjection;

//
// inlines
//

inline
NullQN::NullQN(int const* InIter)
{
}

inline
int degree(NullQN const&)
{
   return 1;
}

inline
double trace(NullQN const&)
{
   return 1;
}

inline
int multiplicity(NullQN const&, NullQN const&, NullQN const&)
{
   return 1;
}

inline
double clebsch_gordan(NullQN const& q1, NullQN const& q2, NullQN const& q,
		      NullQNProjection const& m1,  NullQNProjection const& m2,  NullQNProjection const& m)
{
   return 1;
}

inline
double product_coefficient(NullQN const& k1, NullQN const& k2, NullQN const& k,
			   NullQN const& qp, NullQN const& q, NullQN const& qpp)
{
   return 1;
}

inline
double inverse_product_coefficient(NullQN const& k1, NullQN const& k2, NullQN const& k,
				   NullQN const& qp, NullQN const& q, NullQN const& qpp)
{
   return 1;
}

inline
double tensor_coefficient(NullQN const& q1,  NullQN const& q2,  NullQN const& q,
			  NullQN const& k1,  NullQN const& k2,  NullQN const& k,
			  NullQN const& q1p, NullQN const& q2p, NullQN const& qp)
{
   return 1;
}

inline
double inverse_tensor_coefficient(NullQN const& q1,  NullQN const& q2,  NullQN const& q,
				  NullQN const& k1,  NullQN const& k2,  NullQN const& k,
				  NullQN const& q1p, NullQN const& q2p, NullQN const& qp)
{
   return 1;
}

inline
double recoupling(NullQN const& q1, NullQN const& q2, NullQN const& q12,
		  NullQN const& q3, NullQN const& q, NullQN const& q23)
{
   return 1;
}

inline
double recoupling_12_3__13_2(NullQN const& q1, NullQN const& q2, NullQN const& q12,
                             NullQN const& q3, NullQN const& q, NullQN const& q13)
{
   return 1;
}

inline
NullQN adjoint(NullQN const& q)
{
   return NullQN();
}

inline
double adjoint_coefficient(NullQN const& qp, NullQN const& k, NullQN const& q)
{
   return 1;
}

inline
double conj_phase(NullQN const& qp, NullQN const& k, NullQN const& q)
{
   return 1;
}

inline
bool is_transform_target(NullQN const& q1, NullQN const& q2, NullQN const& q)
{
   return true;
}

inline
int num_transform_targets(NullQN const& q1, NullQN const& q2)
{
   return 1;
}

template <typename OutIter>
inline
void transform_targets(NullQN const& q1, NullQN const& q2, OutIter Out)
{
   *Out++ = NullQN();
}

inline
int num_inverse_transform_targets(NullQN const& q1, NullQN const& q)
{
   return 1;
}

template <typename OutIter>
inline
void inverse_transform_targets(NullQN const& q1, NullQN const& q, OutIter Out)
{
   *Out++ = NullQN();
}

template <typename OutIter>
inline
void enumerate_projections(NullQN const& q, OutIter Out)
{
   *Out++ = NullQNProjection();
}

inline
bool is_delta(NullQN const& q1, NullQN const& Q, NullQNProjection const& P, NullQN const& q2)
{
   return true;
}

inline
NullQNProjection difference(NullQN const& q1, NullQN const& q2)
{
   return NullQNProjection();
}

inline
NullQNProjection negate(NullQNProjection const& q1)
{
   return NullQNProjection();
}

inline
NullQNProjection sum(NullQNProjection const& q1, NullQNProjection const& q2)
{
   return NullQNProjection();
}

inline
bool is_possible(NullQN const& q, NullQNProjection const& p)
{
   return true;
}

inline
bool is_projection(NullQN const& q, NullQNProjection const& p)
{
   return true;
}

inline
NullQN change(NullQN const& q, NullQNProjection const& p)
{
   return NullQN();
}

inline
NullQN heighest_weight(NullQNProjection const& p)
{
   return NullQN();
}

inline
double weight(NullQNProjection const& p)
{
   return 0;
}

inline
double delta_shift_coefficient(NullQN const& qp, NullQN const& k, NullQN const& q, NullQN const& Delta)
{
   return 1;
}

inline
double casimir(NullQN const&, int)
{
   return 0;
}

//
// NullQNProjection
//

namespace
{
   NiftyCounter::nifty_counter<NullQN::Register> NullQNHalfIntCounter;
}

} // namespace QuantumNumbers

#endif

