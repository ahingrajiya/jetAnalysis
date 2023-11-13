/**
 * @file EventCut.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief Example of the event cut
 * @version 0.1
 * @date 2023-10-19
 * 
 * @copyright Copyright (c) 2023
 * 
 */

// Jet analysis headers
#include "EventCut.h"

// ROOT headers
#include "TMath.h"
#include "TString.h"

// C++ headers
#include <iostream>

ClassImp(EventCut)

//________________
EventCut::EventCut() : fVx{-1e9, 1e9}, fVy{-1e9, 1e9}, fVz{-1e9, 1e9},
    fShiftVx{0}, fShiftVy{0}, fVR{1e9}, 
    fHiBin{-1000, 1000}, fCentVal{-1000., 1000.},
    fPtHat{-1e9, 1e9}, fPtHatWeight{-1e9, 1e9}, fVerbose{kFALSE},
    fPPrimaryVertexFilter{kFALSE},
    fHBHENoiseFilterResultRun2Loose{kFALSE},
    fCollisionEventSelectionAODc2{kFALSE},
    fPhfCoincFilter2Th4{kFALSE},
    fPPAprimaryVertexFilter{kFALSE},
    fPBeamScrapingFilter{kFALSE},
    fPClusterCompatibilityFilter{kFALSE},
    fEventsPassed{0}, fEventsFailed{0} {
    fLumi[0] = 0;
    fLumi[1] = std::numeric_limits<unsigned int>::max();
}

//________________
EventCut::~EventCut() {
    /* Empty */
}

//________________
void EventCut::report() {
    TString report = "\nReporting from EventCut";
    report += TString::Format( "Vx              :\t %f - %f\n", fVx[0], fVx[1] );
    report += TString::Format( "Vy              :\t %f - %f\n", fVy[0], fVy[1] );
    report += TString::Format( "Vz              :\t %f - %f\n", fVz[0], fVz[1] );
    report += TString::Format( "HiBin           :\t %d - %d\n", fHiBin[0], fHiBin[1] );
    report += TString::Format( "Centrality      :\t %f - %f\n", fCentVal[0], fCentVal[1] );
    report += TString::Format( "pThat           :\t %f - %f\n", fPtHat[0], fPtHat[1] );
    report += TString::Format( "ptHatWeight     :\t %f - %f\n", fPtHatWeight[0], fPtHatWeight[1] );
    report += TString::Format( "Events passed   :\t %lld\n", fEventsPassed );
    report += TString::Format( "Events failed   :\t %lld\n", fEventsFailed );
    std::cout << report.Data() << std::endl;
}

//________________
Bool_t EventCut::pass(const Event* ev) {
    
    if (fVerbose) {
        std::cout << "\n----- Event cut -----\n";
    }

    const Bool_t goodVx = kTRUE;
    const Bool_t goodVy = kTRUE;

    const Bool_t goodVz = ( fVz[0] <= ev->vz() ) &&
                          ( ev->vz() < fVz[1] );
    if (fVerbose) {
        std::cout << Form("vz        : %5.2f <= %5.2f < %5.2f \t %s \n",
                          fVz[0], ev->vz(), fVz[1], ( goodVz ) ? "true" : "false" );
    }

    const Bool_t goodHiBin = ( fHiBin[0] <= ev->hiBin() ) &&
                             ( ev->hiBin() < fHiBin[1] );
    if (fVerbose) {
        std::cout << Form("hiBin        : %d <= %d < %d \t %s \n",
                          fHiBin[0], ev->hiBin(), fHiBin[1], ( goodHiBin ) ? "true" : "false" );
    }
    const Bool_t goodCent = ( fCentVal[0] <=  ev->centrality() ) &&
                            ( ev->centrality() < fCentVal[1] );

    if (fVerbose) {
        std::cout << Form("centrality   : %5.2f <= %5.2f < %5.2f \t %s \n",
                          fCentVal[0], ev->centrality(), fCentVal[1], ( goodCent ) ? "true" : "false" );
    }

    const Bool_t goodPtHat = ( fPtHat[0] <= ev->ptHat() ) &&
                             ( ev->ptHat() < fPtHat[1] );

    if (fVerbose) {
        std::cout << Form("ptHat        : %9.2f <= %9.2f < %9.2f \t %s \n",
                          fPtHat[0], ev->ptHat(), fPtHat[1], ( goodPtHat ) ? "true" : "false" );
    }

    const Bool_t goodPtHatWeight = ( fPtHatWeight[0] <= ev->ptHatWeight() ) &&
                                   ( ev->ptHatWeight() < fPtHatWeight[1] );
                
    if (fVerbose) {
        std::cout << Form("eventWeight  : %7.2f <= %7.2f < %7.2f \t %s \n",
                          fPtHatWeight[0], ev->ptHatWeight(), fPtHatWeight[1], ( goodPtHatWeight ) ? "true" : "false" );
    }

    Bool_t goodFilters = kTRUE;
    if ( fPPrimaryVertexFilter ) {
        if ( ev->trigAndSkim()->pPAprimaryVertexFilter() == 0 ) goodFilters = kFALSE;
    }
    if ( fHBHENoiseFilterResultRun2Loose ) {
        if ( ev->trigAndSkim()->HBHENoiseFilterResultRun2Loose() == 0 ) goodFilters = kFALSE;
    }
    if ( fCollisionEventSelectionAODc2 ) {
        if ( ev->trigAndSkim()->collisionEventSelectionAODv2() == 0 ) goodFilters = kFALSE;
    }
    if ( fPhfCoincFilter2Th4 ) {
        if ( ev->trigAndSkim()->phfCoincFilter2Th4() == 0 ) goodFilters = kFALSE;
    }
    if ( fPPAprimaryVertexFilter ) {
        if ( ev->trigAndSkim()->pPAprimaryVertexFilter() == 0 ) goodFilters = kFALSE;
    }
    if ( fPBeamScrapingFilter ) {
        if ( ev->trigAndSkim()->pBeamScrapingFilter() == 0 ) goodFilters = kFALSE;
    }
    if ( fPClusterCompatibilityFilter ) {
        if ( ev->trigAndSkim()->pClusterCompatibilityFilter() == 0 ) goodFilters = kFALSE;
    }
    if (fVerbose) {
        std::cout << Form("Event filters passed: %s\n", (goodFilters) ? "true" : "false");
    }

    Bool_t passEvent = goodVx && goodVy && goodVz && goodHiBin &&
                       goodCent && goodPtHat && goodPtHatWeight;
    ( passEvent ) ? fEventsPassed++ : fEventsFailed++;
    
    return passEvent;
}