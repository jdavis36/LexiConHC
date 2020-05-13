#ifndef LEXICONHC_COUPLINGS_H
#define LEXICONHC_COUPLINGS_H

#include <math.h>
#include <string>

// Default values of parameters as in mod_Parameters.F90
#define DEFVAL_MZ 91.1876
#define DEFVAL_MW 80.399
#define DEFVAL_SW 0.23119
#define DEFVAL_CW sqrt (1- pow(DEFVAL_SW,2))
#define DEFVAL_LAMBDA_VI 10000.
#define DEFVAL_e .3028619
#define DEFVAL_gs 1

// Coupling definitions
// Added gluon couplings that appear in H->VV
#define AMPLITUDE_JHUGEN_COUPLING_COMMANDS \
  COUPLING_COMMAND(ghz1, ampjhu, 0.) \
  COUPLING_COMMAND(ghz1_prime2, ampjhu, 0.) \
  COUPLING_COMMAND(ghz2, ampjhu, 0.) \
  COUPLING_COMMAND(ghz4, ampjhu, 0.) \
  COUPLING_COMMAND(ghw1, ampjhu, 0.) \
  COUPLING_COMMAND(ghw1_prime2, ampjhu, 0.) \
  COUPLING_COMMAND(ghw2, ampjhu, 0.) \
  COUPLING_COMMAND(ghw4, ampjhu, 0.) \
  COUPLING_COMMAND(ghzgs1_prime2, ampjhu, 0.) \
  COUPLING_COMMAND(ghzgs2, ampjhu, 0.) \
  COUPLING_COMMAND(ghzgs4, ampjhu, 0.) \
  COUPLING_COMMAND(ghgsgs2, ampjhu, 0.) \
  COUPLING_COMMAND(ghgsgs4, ampjhu, 0.) \
  COUPLING_COMMAND(ghg2, ampjhu, 0.) \
  COUPLING_COMMAND(ghg4, ampjhu, 0.) \
  COUPLING_COMMAND(dV_Z, ampjhu, 0.) \
  COUPLING_COMMAND(dV_A, ampjhu, 0.) \
  COUPLING_COMMAND(dP_Z, ampjhu, 0.) \
  COUPLING_COMMAND(dP_A, ampjhu, 0.) \
  COUPLING_COMMAND(dM_Z, ampjhu, 0.) \
  COUPLING_COMMAND(dM_A, ampjhu, 0.) \
  COUPLING_COMMAND(dFour_Z, ampjhu, 0.) \
  COUPLING_COMMAND(dFour_A, ampjhu, 0.) \
  COUPLING_COMMAND(dZZWpWm, ampjhu, 0.) \
  COUPLING_COMMAND(dZAWpWm, ampjhu, 0.) \
  COUPLING_COMMAND(dAAWpWm, ampjhu, 0.) \

#define EFT_JHUGEN_COUPLING_COMMANDS \
  COUPLING_COMMAND(ghz1, eftjhu, 0.) \
  COUPLING_COMMAND(ghz1_prime2, eftjhu, 0.) \
  COUPLING_COMMAND(ghz2, eftjhu, 0.) \
  COUPLING_COMMAND(ghz4, eftjhu, 0.) \
  COUPLING_COMMAND(ghzgs2, eftjhu, 0.) \
  COUPLING_COMMAND(ghzgs4, eftjhu, 0.) \
  COUPLING_COMMAND(ghgsgs2, eftjhu, 0.) \
  COUPLING_COMMAND(ghgsgs4, eftjhu, 0.) \
  COUPLING_COMMAND(ghg2, eftjhu, 0.) \
  COUPLING_COMMAND(ghg4, eftjhu, 0.) \

// Renamed the EFT_Higgs basis couplings to match rosetta naming convention
#define HIGGSBASIS_COUPLING_COMMANDS \
  COUPLING_COMMAND(dCz, hbasis, 0.) \
  COUPLING_COMMAND(Czz, hbasis, 0.) \
  COUPLING_COMMAND(Czbx, hbasis, 0.) \
  COUPLING_COMMAND(tCzz, hbasis, 0.) \
  COUPLING_COMMAND(dCw, hbasis, 0.) \
  COUPLING_COMMAND(Cww, hbasis, 0.) \
  COUPLING_COMMAND(Cwbx, hbasis, 0.) \
  COUPLING_COMMAND(tCww, hbasis, 0.) \
  COUPLING_COMMAND(Cza, hbasis, 0.) \
  COUPLING_COMMAND(tCza, hbasis, 0.) \
  COUPLING_COMMAND(Cabx, hbasis, 0.) \
  COUPLING_COMMAND(Caa, hbasis, 0.) \
  COUPLING_COMMAND(tCaa, hbasis, 0.) \
  COUPLING_COMMAND(Cgg, hbasis, 0.) \
  COUPLING_COMMAND(tCgg, hbasis, 0.)

#define EFT_HIGGSBASIS_COUPLING_COMMANDS \
  COUPLING_COMMAND(dCz, efthbasis, 0.) \
  COUPLING_COMMAND(Czbx, efthbasis, 0.) \
  COUPLING_COMMAND(Czz, efthbasis, 0.) \
  COUPLING_COMMAND(tCzz, efthbasis, 0.) \
  COUPLING_COMMAND(Cza, efthbasis, 0.) \
  COUPLING_COMMAND(tCza, efthbasis, 0.) \
  COUPLING_COMMAND(Caa, efthbasis, 0.) \
  COUPLING_COMMAND(tCaa, efthbasis, 0.) \
  COUPLING_COMMAND(Cgg, efthbasis, 0.) \
  COUPLING_COMMAND(tCgg, efthbasis, 0.)

namespace LexiConHCCouplings{
#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) coupl_##PREFIX##_##NAME,

  enum Amplitude_JHUGen_CouplingType{
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS
    nAmplitude_JHUGen_CouplingTypes
  };
  enum EFT_JHUGen_CouplingType{
    EFT_JHUGEN_COUPLING_COMMANDS
    nEFT_JHUGen_CouplingTypes
  };
  enum HiggsBasis_CouplingType{
    HIGGSBASIS_COUPLING_COMMANDS
    nHiggsBasis_CouplingTypes
  };
  enum EFT_HiggsBasis_CouplingType{
    EFT_HIGGSBASIS_COUPLING_COMMANDS
    nEFT_HiggsBasis_CouplingTypes
  };
  

#undef COUPLING_COMMAND

  // Functions to get the coupling names from the indices
  std::string getCouplingName(LexiConHCCouplings::Amplitude_JHUGen_CouplingType type);
  std::string getCouplingName(LexiConHCCouplings::EFT_JHUGen_CouplingType type);
  std::string getCouplingName(LexiConHCCouplings::HiggsBasis_CouplingType type); 
  std::string getCouplingName(LexiConHCCouplings::EFT_HiggsBasis_CouplingType type);
}


#endif
