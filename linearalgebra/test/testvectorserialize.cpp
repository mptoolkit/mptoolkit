// -*- C++ -*- $Id$

#include "linearalgebra/vector.h"
//#include "linearalgebra/sparsematrix.h"
#include "linearalgebra/pstreamio.h"
#include "pstream/pfilestream.h"
#include "linearalgebra/vector_utility.h"

using namespace LinearAlgebra;

int main()
{
   Vector<double> M = random_vector<double>(10);

   {
      PStream::opfilestream Out;
      Out.open("hello.test.xdr", O_WRONLY | O_TRUNC | O_CREAT);
      Out.set_format(PStream::format::XDR);
      Out << M;
      Out.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      Vector<double> MCheck;
      In >> MCheck;
      CHECK_CLOSE(M, MCheck);
      In.close();
   }

   {
      PStream::opfilestream Out;
      Out.open("hello.test.xdr", O_WRONLY | O_TRUNC | O_CREAT);
      Out.set_format(PStream::format::XDR);
      Out << M[range(2,8)];
      Out.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      Vector<double> MCheck;
      In >> MCheck;
      CHECK_CLOSE(M[range(2,8)], MCheck);
      In.close();
   }

   // test complex
   Vector<std::complex<double> > V = random_vector<std::complex<double> >(10);

   {
      PStream::opfilestream Out;
      Out.open("hello.test.xdr", O_WRONLY | O_TRUNC | O_CREAT);
      Out.set_format(PStream::format::XDR);
      Out << V;
      Out.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      Vector<std::complex<double> > VCheck;
      In >> VCheck;
      CHECK_CLOSE(V, VCheck);
      In.close();
   }

   // test real/imag
   {
      PStream::opfilestream Out;
      Out.open("hello.test.xdr", O_WRONLY | O_TRUNC | O_CREAT);
      Out.set_format(PStream::format::XDR);
      Out << real(V);
      Out << imag(V);
      Out.close();
   }

   {
      PStream::ipfilestream In;
      In.open("hello.test.xdr", O_RDONLY);
      In.set_format(PStream::format::XDR);
      Vector<double> VCheck;
      In >> VCheck;
      CHECK_CLOSE(real(V), VCheck);
      In >> VCheck;
      CHECK_CLOSE(imag(V), VCheck);
      In.close();
   }

}
