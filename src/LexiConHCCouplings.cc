#include <cassert>
#include <iostream>
#include <string>
#include "LexiConHCCouplings.h"


using namespace std;


#define COUPLING_COMMAND(NAME, PREFIX, DEFVAL) \
case coupl_##PREFIX##_##NAME: \
  return #NAME;

std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::Amplitude_JHUGen_CouplingType type){
  switch (type){
    AMPLITUDE_JHUGEN_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::EFT_JHUGen_CouplingType type){
  switch (type){
    EFT_JHUGEN_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::HiggsBasis_CouplingType type){
  switch (type){
    HIGGSBASIS_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::EFT_HiggsBasis_CouplingType type){
  switch (type){
    EFT_HIGGSBASIS_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::Amplitude_JHUGen_Include_Triple_CouplingType type){
  switch (type){
    AMPLITUDE_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::EFT_JHUGen_Include_Triple_CouplingType type){
  switch (type){
    EFT_JHUGEN_INCLUDE_TRIPLE_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::HiggsBasis_Include_Triple_CouplingType type){
  switch (type){
    HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::EFT_HiggsBasis_Include_Triple_CouplingType type){
  switch (type){
    EFT_HIGGSBASIS_INCLUDE_TRIPLE_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
// Charged and Neutral Current options (Will only be seen on output)
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::Amplitude_JHUGen_Neutral_Current_CouplingType type){
  switch (type){
    AMPLITUDE_JHUGEN_NEUTRAL_CURRENT_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::Amplitude_JHUGen_Charged_Current_CouplingType type){
  switch (type){
    AMPLITUDE_JHUGEN_CHARGED_CURRENT_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::Amplitude_JHUGen_Include_Triple_Charged_Current_CouplingType type){
  switch (type){
    AMPLITUDE_JHUGEN_INCLUDE_TRIPLE_CHARGED_CURRENT_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}
std::string LexiConHCCouplings::getCouplingName(LexiConHCCouplings::Amplitude_JHUGen_Include_Triple_Neutral_Current_CouplingType type){
  switch (type){
    AMPLITUDE_JHUGEN_INCLUDE_TRIPLE_NEUTRAL_CURRENT_COUPLING_COMMANDS
  default:
    cerr << "LexiConHCCouplings::getCouplingName: Type " << type << " is undefined." << endl;
    assert(0);
  }
}





#undef COUPLING_COMMAND
