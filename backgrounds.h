auto primaryMuon = [](int pdg, std::string process, bool matched, int origin ){
  return ( abs(pdg) == 13 && process == "primary" && matched && origin == 4 );
};

auto isDecay = [](std::string process, int origin){
  return( origin == 4 && process == "Decay" );
};

auto isCosmic = [](int origin){
  return( origin == 2 );
};

auto isExtraBeam = [](std::string process, bool matched, int origin){
  return( process == "primary" && !matched && origin == 4 );
};

auto upstreamInt = [](std::string process, int origin){
  return( ( process.find("Inelastic") != std::string::npos ) && origin == 4 );
};
