/**
 * @file GenTrack.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Generated track description
 * @version 0.1
 * @date 2023-10-19
 *
 * @copyright Copyright (c) 2023
 *
 */

// Jet analysis headers
#include "GenTrack.h"
#include "TString.h"
#include <iostream>
ClassImp(GenTrack)

    //___________________
    GenTrack::GenTrack() : BaseTrack{}, fPDGID{0}, fSube{0}
{
}

void GenTrack::print()
{
    std::cout << "-- Generated Tracks Info --\n"
              << Form("TrkPt: %5.2f TrkEta: %3.2f TrkPhi: %3.2f PDGID: %i Sube: %c", this->TrkPt(), this->TrkEta(), this->TrkPhi(), this->TrkPDGID(), this->TrkSube());
}
