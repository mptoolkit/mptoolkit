// -*- C++ -*-
//----------------------------------------------------------------------------
// Matrix Product Toolkit http://mptoolkit.qusim.net/
//
// utils/SSF-FourierTransform-YC_OldNumbering.cpp
//
// Copyright (C) 2016 Seyed Saadatmand <s.saadatmand@uq.edu.au>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Research publications making use of this software should include
// appropriate citations and acknowledgements as described in
// the file CITATIONS in the main source directory.
//----------------------------------------------------------------------------
// ENDHEADER

// -*- C++ -*- $Id: SSF-DiscreteFourierTransform.cpp 1490 2015-05-19 09:15:06Z seyed $
// Authors: Ian P. McCulloch and Seyed N. Saadatmand
// Contact: s.saadatmand@uq.edu.au

#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>
#include "common/terminal.h"
#include "common/math_const.h"
#include <boost/program_options.hpp>
#include "mp/copyright.h"

using namespace std;

int main(int argc, char** argv)
{

   // Default values of variables
   //int GridSize = 100;
   //int w = 0;
   //string DataFile;
    
   if (argc != 4)
     {
         cerr << "usage: " << basename(argv[0]) << " <width> <grid size> <real-space correlation file> [all required]\n\n"
              << "description: Performs a planar discrete Fourier transform of real-space correlations on a\n"
              << "           \"YC_OldNumbering\" structure into a 2D momentum space (inverse lattice space).\n\n"
              << "expected file format for \"real-space correlation file\":\n"
              << "    <n2> <correlation>\n\n"
              << "outputs:\n"
              << "[default output]	<k_x> <K_y> <SSF>\n";
         return 1;
     }

   cout.precision(16);

   int w = boost::lexical_cast<int>(argv[1]);
   int GridSize = boost::lexical_cast<int>(argv[2]);
   string DataFile = argv[3];

   vector<int> n1, n2;
   vector<double> corr;

   ifstream file( DataFile.c_str() );

   int my_n1, my_n2;
   double my_corr;
   while (file >> my_n1 >> my_n2 >> my_corr) 
   {
      n1.push_back(my_n1);
      n2.push_back(my_n2);
      corr.push_back(my_corr);
   }

   double N = 0.5 * ( -1 + sqrt( 1.0 + 8*corr.size() ) );
   //cout << "N=" << N << endl; //just a test 

   for (int i_kx = 0; i_kx < GridSize; ++i_kx)
   {
      double kx = 2*math_const::pi*( i_kx*(2.0/(GridSize-1)) - 1.0 );

      for (int i_ky = 0; i_ky < GridSize; ++i_ky)
      {
	 double ky = 2*math_const::pi*( i_ky*(2.0/(GridSize-1)) - 1.0 );

	 double sum = 0;

         // performing summation,
	 for (unsigned i = 0; i < corr.size(); ++i)
	 {
	    // get real-space coordinates of the points n1[i] and n2[i]

	    int cell1 = n1[i]/w;
	    int site1 = n1[i]%w;

	    int cell2 = n2[i]/w;
	    int site2 = n2[i]%w;

	    double x1 = cell1 * sqrt(3.0) / 2.0;
	    double y1 = site1 + (cell1/2.0);

	    double x2 = cell2 * sqrt(3.0) / 2.0;
	    double y2 = site2 + (cell2/2.0);
	    
	    sum += 2 * cos( kx*(x1-x2) + ky*(y1-y2) ) * corr[i];
            //cout << "for n1=" << n1[i] << ", n2=" << n2[i] << ", and corr=" << corr[i] << " ---> sum=" << sum << endl; //just a test     
	 }
          
	 double SSF = sum/N - 0.75;

	 cout << kx << '\t' << ky << '\t' << abs(SSF) << endl;

      }
   
    cout << '\n';

   }

   return 0;
}
