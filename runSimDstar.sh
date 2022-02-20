#tion!/bin/bash
parallel -j50 "root -b -q -l 'SimulateDstarPolarization.cc(1000000, kEvtGen, kMonash, kHardQCD, 13000, "{}", \"AnalysisResults_DstarPolarization_Monash_HardQCD_EvtGen_"{}".root\")'" ::: {1..50}
