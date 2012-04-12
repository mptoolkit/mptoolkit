// -*- C++ -*- $Id$


      double Broadening = 0;
      double GroundstateEnergy;


         ("GroundstateEnergy,G", prog_opt::value<double>(&GroundstateEnergy),
          "Groundstate energy of the Hamiltonian (wavefunction attribute \"GroundstateEnergy\")")

         ("Frequency,F", prog_opt::value<std::vector<double> >(),
          "Frequency of the corresponding correction vector (not required if the \"Frequency\" attribute is set")
         ("Broadening,B", prog_opt::value<double>(&Broadening),
          "Magnitude of the Lorentzian broadening to use (wavefunction attribute \"Broadening\")")


      // Determine the groundstate energy
      if (vm.count("GroundstateEnergy") == 0)
      {
         if (Psi[0].Attributes().count("GroundstateEnergy") == 0)
         {
            std::cerr << "fatal: No groundstate energy specified.\n"
               "Specity the energy either with --GroundstateEnergy parameter, "
               "or as an attribute of the first wavefunction.\n";
            return 1;
         }
         GroundstateEnergy = Psi[0].Attributes()["GroundstateEnergy"].as<double>();
      }
