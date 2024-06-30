/**
 * @file ForestAODReader.cc
 * @author Grigory Nigmatkulov (gnigmat@uic.edu)
 * @brief CMS ForestAOD reader
 * @version 0.1
 * @date 2023-10-19
 *
 * @copyright Copyright (c) 2023
 *
 */

// Jet analysis headers
#include "ForestAODReader.h"

// ROOT headers
#include "TFile.h"

// C++ headers
#include <cstddef>
#include <cstring>
#include <fstream>

ClassImp(ForestAODReader)

    //_________________
    ForestAODReader::ForestAODReader()
    : fEvent{nullptr}, fInFileName{nullptr}, fEvents2Read{0}, fEventsProcessed{0}, fIsMc{kFALSE}, fUseHltBranch{kTRUE},
      fUseSkimmingBranch{kTRUE}, fJetCollection{"ak4PFJetAnalyzer"}, fUseJets{kTRUE}, fUseTrackBranch{kFALSE},
      fUseGenTrackBranch{kFALSE}, fHltTree{nullptr}, fSkimTree{nullptr}, fEventTree{nullptr},
      fTrkTree{nullptr}, fGenTrkTree{nullptr},
      fJEC{nullptr}, fJECFiles{}, fJEU{nullptr}, fJEUFiles{},
      fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020}, fYearOfDataTaking{2018}, fDoJetPtSmearing{kFALSE},
      fFixJetArrays{kFALSE}, fEventCut{nullptr}, fJetCut{nullptr}, fRecoJet2GenJetId{}, fGenJet2RecoJet{}
{
    // Initialize many variables
    clearVariables();
}

//_________________
ForestAODReader::ForestAODReader(const Char_t *inputStream, const Bool_t &useHltBranch, const Bool_t &useSkimmingBranch,
                                 const Char_t *jetCollection, const Bool_t &useJets,
                                 const Bool_t &useTrackBranch, const Bool_t &useGenTrackBranch, const Bool_t &isMc, const Bool_t &setStoreLocation)
    : fEvent{nullptr}, fInFileName{inputStream}, fEvents2Read{0}, fEventsProcessed{0}, fIsMc{isMc},
      fUseHltBranch{useHltBranch}, fUseSkimmingBranch{useSkimmingBranch}, fJetCollection{"ak4PFJetAnalyzer"},
      fUseJets{useJets}, fUseTrackBranch{useTrackBranch}, fUseGenTrackBranch{useGenTrackBranch},
      fJEC{nullptr}, fJECFiles{}, fJEU{nullptr},
      fJEUFiles{}, fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020}, fYearOfDataTaking{2018},
      fDoJetPtSmearing{kFALSE}, fFixJetArrays{kFALSE}, fEventCut{nullptr}, fJetCut{nullptr}, fIsInStore{setStoreLocation}
{
    // Initialize many variables
    clearVariables();
}

//_________________
ForestAODReader::~ForestAODReader()
{
    if (fEvent)
        delete fEvent;
    if (fHltTree)
        delete fHltTree;
    if (fSkimTree)
        delete fSkimTree;
    if (fEventTree)
        delete fEventTree;
    if (fUseJets)
        delete fJetTree;
    if (fTrkTree)
        delete fTrkTree;
    if (fGenTrkTree)
        delete fGenTrkTree;
    if (fJEC)
        delete fJEC;
    if (fJEU)
        delete fJEU;
    if (fEventCut)
        delete fEventCut;
    if (fJetCut)
        delete fJetCut;
}

//_________________
void ForestAODReader::clearVariables()
{
    fRunId = {0};
    fEventId = {0};
    fLumi = {0};
    fVertexZ = {-999.f};
    fHiBin = {-1};
    fPtHatWeight = {-1.f};
    fPtHat = {-1.f};

    // bad jets and multiplicity to be added

    fNRecoJets = {0};
    fNGenJets = {0};
    fNTracks = {0};

    fHLT_HIAK4PFJet15_v1 = {0};
    fHLT_HIAK4PFJet15_v1_Prescl = {0};
    fHLT_HIAK4PFJet30_v1 = {0};
    fHLT_HIAK4PFJet30_v1_Prescl = {0};
    fHLT_HIAK4PFJet40_v1 = {0};
    fHLT_HIAK4PFJet40_v1_Prescl = {0};
    fHLT_HIAK4PFJet60_v1 = {0};
    fHLT_HIAK4PFJet60_v1_Prescl = {0};
    fHLT_HIAK4PFJet80_v1 = {0};
    fHLT_HIAK4PFJet80_v1_Prescl = {0};
    fHLT_HIAK4PFJet120_v1 = {0};
    fHLT_HIAK4PFJet120_v1_Prescl = {0};

    fHLT_HIAK8PFJet15_v1 = {0};
    fHLT_HIAK8PFJet15_v1_Prescl = {0};
    fHLT_HIAK8PFJet25_v1 = {0};
    fHLT_HIAK8PFJet25_v1_Prescl = {0};
    fHLT_HIAK8PFJet40_v1 = {0};
    fHLT_HIAK8PFJet40_v1_Prescl = {0};
    fHLT_HIAK8PFJet60_v1 = {0};
    fHLT_HIAK8PFJet60_v1_Prescl = {0};
    fHLT_HIAK8PFJet80_v1 = {0};
    fHLT_HIAK8PFJet80_v1_Prescl = {0};
    fHLT_HIAK8PFJet140_v1 = {0};
    fHLT_HIAK8PFJet140_v1_Prescl = {0};

    fHLT_HIPFJet25_v1 = {0};
    fHLT_HIPFJet25_v1_Prescl = {0};
    fHLT_HIPFJet140_v1 = {0};
    fHLT_HIPFJet140_v1_Prescl = {0};

    fHLT_HIPuAK4CaloJet80Eta5p1_v1 = {0};
    fHLT_HIPuAK4CaloJet100Eta5p1_v1 = {0};

    fHBHENoiseFilterResultRun2Loose = {0};
    fHBHENoiseFilterResultRun2Tight = {0};
    fHBHEIsoNoiseFilterResult = {0};
    fCollisionEventSelectionAODv2 = {0};
    fPhfCoincFilter2Th4 = {0};
    fPPAprimaryVertexFilter = {0};
    fPBeamScrapingFilter = {0};
    fPprimaryVertexFilter = {0};
    fPVertexFilterCutG = {0};
    fPVertexFilterCutGloose = {0};
    fPVertexFilterCutGtight = {0};
    fPVertexFilterCutE = {0};
    fPVertexFilterCutEandG = {0};
    fPClusterCompatibilityFilter = {0};

    // Loop over jets and tracks
    for (Short_t i{0}; i < 9999; i++)
    {

        // Jet variables
        if (i < 100)
        {
            fRawJetPt[i] = {0.f};
            fRecoJetPt[i] = {0.f};
            fRecoJetEta[i] = {0.f};
            fRecoJetPhi[i] = {0.f};
            fRecoJetWTAEta[i] = {0.f};
            fRecoJetWTAPhi[i] = {0.f};
            fRecoJetTrackMax[i] = {0.f};
            if (fIsMc)
            {
                fRefJetPt[i] = {0.f};
                fRefJetEta[i] = {0.f};
                fRefJetPhi[i] = {0.f};
                fRefJetWTAEta[i] = {0.f};
                fRefJetWTAPhi[i] = {0.f};
                fRefJetPartonFlavor[i] = {-999};
                fRefJetPartonFlavorForB[i] = {-99};
                fGenJetPt[i] = {0.f};
                fGenJetEta[i] = {0.f};
                fGenJetPhi[i] = {0.f};
                fGenJetWTAEta[i] = {0.f};
                fGenJetWTAPhi[i] = {0.f};
            }

        } // if (i<100

        // Track variables
        fTrackPt[i] = {0.f};
        fTrackEta[i] = {0.f};
        fTrackPhi[i] = {0.f};
        fTrackPtErr[i] = {0.f};
        fTrackDcaXY[i] = {0.f};
        fTrackDcaZ[i] = {0.f};
        fTrackDcaXYErr[i] = {0.f};
        fTrackDcaZErr[i] = {0.f};
        fTrackChi2[i] = {0.f};
        fTrackNDOF[i] = {0};
        fTrackPartFlowEcal[i] = {0.f};
        fTrackPartFlowHcal[i] = {0.f};
        fTrackMVA[i] = {0.f};
        fTrackAlgo[i] = {0};
        fTrackCharge[i] = {0};
        fTrackNHits[i] = {0};
        fTrackNLayers[i] = {0};
        fTrackHighPurity[i] = {kFALSE};
    } // for (Short_t i{0}; i<9999; i++)

    fGenTrackPt.clear();
    fGenTrackEta.clear();
    fGenTrackPhi.clear();
    fGenTrackCharge.clear();
    fGenTrackPid.clear();
    fGenTrackSube.clear();

    if (fIsMc)
    {
        fRecoJet2GenJetId.clear();
        fGenJet2RecoJet.clear();
    }
}

//_________________
Int_t ForestAODReader::init()
{
    Int_t status = 0;
    // Setup chains to read
    status = setupChains();
    // Setup branches to read
    setupBranches();
    // Setup jet energy correction files and pointer
    setupJEC();
    return status;
}

//________________
void ForestAODReader::setupJEC()
{
    if (fJECFiles.empty())
    {
        if (fJECInputFileName.Length() <= 0)
        {
            setJECFileName();
        }
        fJECFiles.push_back(Form("/Users/gnigmat/work/cms/soft/jetAnalysis/aux_files/%s_%i/JEC/%s",
                                 fCollidingSystem.Data(), fCollidingEnergyGeV, fJECInputFileName.Data()));
        std::cout << Form("Add JEC file: %s\n", fJECFiles.back().c_str());
    }

    fJEC = new JetCorrector(fJECFiles);
}

//________________
Float_t ForestAODReader::jetPtWeight(const Bool_t &isMC, const std::string &system, const Int_t &year,
                                     const Int_t &energy, float jetpt) const
{
    Float_t jetptweight = 1.0;

    // JetPtWeightFunction is derived from MC vs data jet pT spectra.
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *JetPtWeightFunction = new TF1("JetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from all jets above 120
    GeV and JECv6 JetPtWeightFunction->SetParameters(0.79572,0.0021861,-6.35407e-06,6.66435e-09); jetptweight =
    JetPtWeightFunction->Eval(jetpt);
    }
    */
    return jetptweight;
}

//________________
Float_t ForestAODReader::leadJetPtWeight(const Bool_t &isMC, const std::string &system, const Int_t &year,
                                         const Int_t &energy, const Float_t &leadjetpt) const
{
    Float_t leadjetptweight = 1.0;

    // LeadJetPtWeightFunction is derived from MC vs data leading jet pT spectra.
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *LeadJetPtWeightFunction = new TF1("LeadJetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from leading
    jets above 120 GeV and JECv6 LeadJetPtWeightFunction->SetParameters(0.876682,0.00131479,-3.90884e-06,4.40358e-09); ;
        leadjetptweight = LeadJetPtWeightFunction->Eval(leadjetpt);
    }
    */
    return leadjetptweight;
}

//________________
Float_t ForestAODReader::subleadJetPtWeight(const Bool_t &isMC, const std::string &system, const Int_t &year,
                                            const Int_t &energy, const Float_t &subleadjetpt)
{
    Float_t subleadjetptweight = 1.0;

    // SubLeadJetPtWeightFunction is derived from MC vs data subleading jet pT spectra.
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *SubLeadJetPtWeightFunction = new TF1("SubLeadJetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from
    leading jets above 120 GeV and JECv6
        SubLeadJetPtWeightFunction->SetParameters(0.876682,0.00131479,-3.90884e-06,4.40358e-09); ;
        subleadjetptweight = SubLeadJetPtWeightFunction->Eval(subleadjetpt);
    }
    */
    return subleadjetptweight;
}

//________________
Float_t ForestAODReader::jetPtSmeringWeight(const Bool_t &isMC, const std::string &system, const Int_t &year,
                                            const Int_t &energy, const Float_t &jetpt, const Bool_t &dosmearing,
                                            const Float_t resolutionfactor) const
{
    Float_t jetptsmearweight = 1.0;
    if (!dosmearing)
        return jetptsmearweight;

    // JetPtSmearingWeightFunction is derived from MC vs data jet pT spectra.
    if (!isMC && system == "pp" && energy == 5020 && year == 2017)
    {
        TF1 *JetPtSmearingWeightFunction =
            new TF1("JetPtSmearingWeightFunction", "pol3", 0.0, 500.0); // Derived from all jets above 120 GeV and JECv6
        JetPtSmearingWeightFunction->SetParameters(0.174881, -0.00091979, 3.50064e-06, -6.52541e-09, 4.64199e-12);
        jetptsmearweight = JetPtSmearingWeightFunction->Eval(jetpt);
        jetptsmearweight = jetptsmearweight * resolutionfactor;
    }

    return jetptsmearweight;
}

//________________
Float_t ForestAODReader::trkEtaMixWeight(const Bool_t &isMC, const std::string &system, const Int_t &year,
                                         const Int_t &energy, const Float_t &trketa, const Bool_t &reco) const
{
    float trketamixweight = 1.0;

    // TrkEtaMixWeightFunction is derived from trk eta from signal over trk eta from mixing
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017 && !reco){
        TF1 *TrkEtaMixWeightFunction = new TF1("TrkEtaMixWeightFunction", "pol3", 0.0, 500.0);
        TrkEtaMixWeightFunction->SetParameters(0.174881, -0.00091979, 3.50064e-06, -6.52541e-09, 4.64199e-12);
        trketamixweight = TrkEtaMixWeightFunction->Eval(jetpt);
    }
    */
    return trketamixweight;
}

//________________
Float_t ForestAODReader::eventWeight(const Bool_t &isMC, const Bool_t &use_centrality, const std::string &system,
                                     const Int_t &year, const Int_t &energy, const Float_t &vz, const Int_t mult,
                                     const Float_t &weighttree, const Float_t &leadjetpt) const
{

    Float_t vzweight = 1.0;
    Float_t multweight = 1.0;
    Float_t evtweight = 1.0;
    Float_t multefficiency = 1.0;
    Float_t jetefficiency = 1.0;
    Float_t totalweight = 1.0;

    // VzWeightFunction is derived from MC vs data event Vz --> MC only --> vzweight
    // MultCentWeightFunction is derived from MC vs data event multiplicity or centrality --> MC only --> multweight
    // MultTriggerWeightFunction is derived from the turn on plots as function of multiplicity --> RECO only
    // JetTriggerWeightFunction is derived from the turn on plots as function of leading jet pT --> RECO only
    // weighttree is the pthat weight --> MC only

    if (isMC && !use_centrality && system == "pp" && energy == 5020 && year == 2017)
    {

        TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol6", -15.0, 15.0);
        VzWeightFunction->SetParameters(0.973805, 0.00339418, 0.000757544, -1.37331e-06, -2.82953e-07, -3.06778e-10,
                                        3.48615e-09);
        vzweight = VzWeightFunction->Eval(vz);

        TF1 *MultCentWeightFunction = new TF1("MultCentWeightFunction", "pol0", 0.0, 500.0);
        MultCentWeightFunction->SetParameter(0, 1.0);
        multweight = MultCentWeightFunction->Eval(mult);

        TF1 *MultTriggerWeightFunction =
            new TF1("MultTriggerWeightFunction", "pol0", 0.0, 500.0); // fitted from turn on curves
        MultTriggerWeightFunction->SetParameter(0, 1.0);
        Float_t multtrigweight = 1.0;
        multtrigweight = MultTriggerWeightFunction->Eval(mult);
        multefficiency = 1. / multtrigweight;

        TF1 *JetTriggerWeightFunction =
            new TF1("JetTriggerWeightFunction", "pol0", 0.0, 500.0); // fitted from turn on curves
        JetTriggerWeightFunction->SetParameter(0, 1.0);
        Float_t jettrigweight = 1.0;
        jettrigweight = JetTriggerWeightFunction->Eval(leadjetpt);
        jetefficiency = 1. / jettrigweight;

        evtweight = weighttree;
    }

    totalweight = evtweight * multweight * vzweight * multefficiency * jetefficiency;
    return totalweight;
}

//_________________
void ForestAODReader::finish()
{
}

//_________________
Int_t ForestAODReader::setupChains()
{

    // Setup chains (0-good, 1-bad)
    Int_t returnStatus = 1;

    // Setup chains to read

    std::cout << "Setting chains... ";

    // Use event branch
    fEventTree = new TChain("hiEvtAnalyzer/HiTree");

    // Use HLT branch
    if (fUseHltBranch)
    {
        fHltTree = new TChain("hltanalysis/HltTree");
    }
    // Use skimming branch
    if (fUseSkimmingBranch)
    {
        fSkimTree = new TChain("skimanalysis/HltTree");
    }
    // Use jet branch
    if (fUseJets)
    {
        fJetTree = new TChain(Form("%s/t", fJetCollection.Data()));
    }

    // Use reconstructed track branch
    if (fUseTrackBranch)
    {
        fTrkTree = new TChain("ppTrack/trackTree");
    }
    // Use generated track branch
    if (fIsMc && fUseGenTrackBranch)
    {
        fGenTrkTree = new TChain("HiGenParticleAna/hi");
    }
    std::cout << "\t[DONE]\n";

    // Initialize input file name (should switch to const char* processing later)
    TString input(fInFileName);

    // Check input exists
    if (input.Length() <= 0)
    {
        std::cerr << "No normal inputfile. Terminating." << std::endl;
        returnStatus = 1;
        exit(0);
    }
    // Normail input
    else
    {
        // If input is a single ROOT file
        if (input.Index(".root") > 0)
        {
            std::cout << Form("Adding %s file to chains\n", input.Data());
            fEventTree->Add(input.Data());
            if (fUseHltBranch)
                fHltTree->Add(input.Data());
            if (fUseSkimmingBranch)
                fSkimTree->Add(input.Data());
            if (fUseJets)
                fJetTree->Add(input.Data());
            if (fUseTrackBranch)
                fTrkTree->Add(input.Data());
            if (fIsMc && fUseGenTrackBranch)
                fGenTrkTree->Add(input.Data());

            fEvents2Read = fEventTree->GetEntries();
            std::cout << Form("Total number of events to read: %lld\n", fEvents2Read);
            Long64_t fEvents2Read2 = fJetTree->GetEntries();
            std::cout << Form("Total number of events to read2: %lld\n", fEvents2Read2);
        }
        // Assuming that list of files is provided instead of a single file
        else
        {
            std::ifstream inputStream(input.Data());

            if (!inputStream)
                std::cout << Form("ERROR: Cannot open file list: %s\n", input.Data());
            Int_t nFiles = 0;
            std::string file;
            size_t pos;
            while (getline(inputStream, file))
            {
                // NOTE: our external formatters may pass "file NumEvents"
                //       Take only the first part
                // cout << "DEBUG found " <<  file << endl;
                pos = file.find_first_of(" ");
                if (pos != std::string::npos)
                    file.erase(pos, file.length() - pos);
                // cout << "DEBUG found [" <<  file << "]" << endl;

                // Check that file is of a correct name
                if (file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos &&
                    file.find("AOD") != std::string::npos)
                {
                    // Open file
                    if (fIsInStore)
                    {
                        file = "root://cmsxrootd.fnal.gov/" + file;
                    }

                    // Open file
                    TFile *ftmp = TFile::Open(file.c_str());

                    // Check file is not zombie and contains information
                    if (ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())
                    {
                        std::cout << Form("Adding file to chain: %s\n", file.c_str());
                        // Adding file to chains
                        fEventTree->Add(file.c_str());
                        if (fUseHltBranch)
                            fHltTree->Add(file.c_str());
                        if (fUseSkimmingBranch)
                            fSkimTree->Add(file.c_str());
                        if (fUseJets)
                            fJetTree->Add(file.c_str());
                        if (fUseTrackBranch)
                            fTrkTree->Add(file.c_str());
                        if (fIsMc && fUseGenTrackBranch)
                            fGenTrkTree->Add(file.c_str());
                        ++nFiles;
                    } // if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

                    if (ftmp)
                    {
                        ftmp->Close();
                    } // if (ftmp)
                }     // if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos &&
                      // file.find("AOD") != std::string::npos )
            }         // while ( getline( inputStream, file ) )

            std::cout << Form("Total number of files in chain: %d\n", nFiles);
            fEvents2Read = fEventTree->GetEntries();
            std::cout << Form("Total number of events to read: %lld\n", fEvents2Read);
        } // else {   if file list
        returnStatus = 0;
    } // else {   if normal input
    return returnStatus;
}

//_________________
void ForestAODReader::setupBranches()
{

    // Disable all branches - this is important while reading big files
    fEventTree->SetBranchStatus("*", 0);
    if (fUseHltBranch)
        fHltTree->SetBranchStatus("*", 0);
    if (fUseSkimmingBranch)
        fSkimTree->SetBranchStatus("*", 0);
    if (fUseJets)
        fJetTree->SetBranchStatus("*", 0);
    if (fUseTrackBranch)
        fTrkTree->SetBranchStatus("*", 0);
    if (fUseGenTrackBranch && fIsMc)
        fGenTrkTree->SetBranchStatus("*", 0);

    // enable branches of interest -> see definition of each variables above

    // Event quantities
    fEventTree->SetBranchStatus("run", 1);
    fEventTree->SetBranchStatus("evt", 1);
    fEventTree->SetBranchStatus("lumi", 1);
    fEventTree->SetBranchStatus("vz", 1);
    fEventTree->SetBranchStatus("hiBin", 1); // centrality only for PbPb and XeXe
    fEventTree->SetBranchAddress("run", &fRunId);
    fEventTree->SetBranchAddress("evt", &fEventId);
    fEventTree->SetBranchAddress("lumi", &fLumi);
    fEventTree->SetBranchAddress("vz", &fVertexZ);
    fEventTree->SetBranchAddress("hiBin", &fHiBin);

    if (fIsMc)
    {
        fEventTree->SetBranchStatus("weight", 1);
        fEventTree->SetBranchStatus("pthat", 1);
        fEventTree->SetBranchAddress("weight", &fPtHatWeight);
        fEventTree->SetBranchAddress("pthat", &fPtHat);
    }

    if (fUseHltBranch)
    {
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet15_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet15_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet30_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet30_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet40_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet40_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet60_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet60_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet80_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet80_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet120_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK4PFJet120_v1_Prescl", 1);

        fHltTree->SetBranchStatus("HLT_HIAK8PFJet15_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet15_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet25_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet25_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet40_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet40_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet60_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet60_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet80_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet80_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet140_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIAK8PFJet140_v1_Prescl", 1);

        fHltTree->SetBranchStatus("HLT_HIPFJet25_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIPFJet25_v1_Prescl", 1);
        fHltTree->SetBranchStatus("HLT_HIPFJet140_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIPFJet140_v1_Prescl", 1);

        fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet80Eta5p1_v1", 1);
        fHltTree->SetBranchStatus("HLT_HIPuAK4CaloJet100Eta5p1_v1", 1);

        fHltTree->SetBranchAddress("HLT_HIAK4PFJet15_v1", &fHLT_HIAK4PFJet15_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet15_v1_Prescl", &fHLT_HIAK4PFJet15_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet30_v1", &fHLT_HIAK4PFJet30_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet30_v1_Prescl", &fHLT_HIAK4PFJet30_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet40_v1", &fHLT_HIAK4PFJet40_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet40_v1_Prescl", &fHLT_HIAK4PFJet40_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet60_v1", &fHLT_HIAK4PFJet60_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet60_v1_Prescl", &fHLT_HIAK4PFJet60_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet80_v1", &fHLT_HIAK4PFJet80_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet80_v1_Prescl", &fHLT_HIAK4PFJet80_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet120_v1", &fHLT_HIAK4PFJet120_v1);
        fHltTree->SetBranchAddress("HLT_HIAK4PFJet120_v1_Prescl", &fHLT_HIAK4PFJet120_v1_Prescl);

        fHltTree->SetBranchAddress("HLT_HIAK8PFJet15_v1", &fHLT_HIAK8PFJet15_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet15_v1_Prescl", &fHLT_HIAK8PFJet15_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet25_v1", &fHLT_HIAK8PFJet25_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet25_v1_Prescl", &fHLT_HIAK8PFJet25_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet40_v1", &fHLT_HIAK8PFJet40_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet40_v1_Prescl", &fHLT_HIAK8PFJet40_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet60_v1", &fHLT_HIAK8PFJet60_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet60_v1_Prescl", &fHLT_HIAK8PFJet60_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet80_v1", &fHLT_HIAK8PFJet80_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet80_v1_Prescl", &fHLT_HIAK8PFJet80_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet140_v1", &fHLT_HIAK8PFJet140_v1);
        fHltTree->SetBranchAddress("HLT_HIAK8PFJet140_v1_Prescl", &fHLT_HIAK8PFJet140_v1_Prescl);

        fHltTree->SetBranchAddress("HLT_HIPFJet25_v1", &fHLT_HIPFJet25_v1);
        fHltTree->SetBranchAddress("HLT_HIPFJet25_v1_Prescl", &fHLT_HIPFJet25_v1_Prescl);
        fHltTree->SetBranchAddress("HLT_HIPFJet140_v1", &fHLT_HIPFJet140_v1);
        fHltTree->SetBranchAddress("HLT_HIPFJet140_v1_Prescl", &fHLT_HIPFJet140_v1_Prescl);

        fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet80Eta5p1_v1", &fHLT_HIPuAK4CaloJet80Eta5p1_v1);
        fHltTree->SetBranchAddress("HLT_HIPuAK4CaloJet100Eta5p1_v1", &fHLT_HIPuAK4CaloJet100Eta5p1_v1);
    } // if ( fUseHltBranch )

    // Skimming quantities
    if (fUseSkimmingBranch)
    {
        fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
        fSkimTree->SetBranchStatus("HBHENoiseFilterResultRun2Tight", 1);
        fSkimTree->SetBranchStatus("HBHEIsoNoiseFilterResult", 1);
        fSkimTree->SetBranchStatus("collisionEventSelectionAODv2", 1);
        fSkimTree->SetBranchStatus("phfCoincFilter2Th4", 1);
        fSkimTree->SetBranchStatus("pPAprimaryVertexFilter", 1);
        fSkimTree->SetBranchStatus("pBeamScrapingFilter", 1);
        fSkimTree->SetBranchStatus("pprimaryVertexFilter", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutG", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutGloose", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutGtight", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutE", 1);
        fSkimTree->SetBranchStatus("pVertexFilterCutEandG", 1);
        fSkimTree->SetBranchStatus("pclusterCompatibilityFilter", 1);

        fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &fHBHENoiseFilterResultRun2Loose);
        fSkimTree->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &fHBHENoiseFilterResultRun2Tight);
        fSkimTree->SetBranchAddress("HBHEIsoNoiseFilterResult", &fHBHEIsoNoiseFilterResult);
        fSkimTree->SetBranchAddress("collisionEventSelectionAODv2", &fCollisionEventSelectionAODv2);
        fSkimTree->SetBranchAddress("phfCoincFilter2Th4", &fPhfCoincFilter2Th4);
        fSkimTree->SetBranchAddress("pPAprimaryVertexFilter", &fPPAprimaryVertexFilter);
        fSkimTree->SetBranchAddress("pBeamScrapingFilter", &fPBeamScrapingFilter);
        fSkimTree->SetBranchAddress("pprimaryVertexFilter", &fPprimaryVertexFilter);
        fSkimTree->SetBranchAddress("pVertexFilterCutG", &fPVertexFilterCutG);
        fSkimTree->SetBranchAddress("pVertexFilterCutGloose", &fPVertexFilterCutGloose);
        fSkimTree->SetBranchAddress("pVertexFilterCutGtight", &fPVertexFilterCutGtight);
        fSkimTree->SetBranchAddress("pVertexFilterCutE", &fPVertexFilterCutE);
        fSkimTree->SetBranchAddress("pVertexFilterCutEandG", &fPVertexFilterCutEandG);
        fSkimTree->SetBranchAddress("pclusterCompatibilityFilter", &fPClusterCompatibilityFilter);
    } // if ( fUseSkimmingBranch )

    // Jet quantities
    if (fUseJets)
    {

        fJetTree->SetBranchStatus("nref", 1);
        fJetTree->SetBranchStatus("rawpt", 1);
        fJetTree->SetBranchStatus("jtpt", 1);
        fJetTree->SetBranchStatus("trackMax", 1);
        fJetTree->SetBranchStatus("jteta", 1);
        fJetTree->SetBranchStatus("jtphi", 1);
        fJetTree->SetBranchStatus("WTAeta", 1);
        fJetTree->SetBranchStatus("WTAphi", 1);

        fJetTree->SetBranchAddress("nref", &fNRecoJets);
        fJetTree->SetBranchAddress("rawpt", &fRawJetPt);
        fJetTree->SetBranchAddress("jtpt", &fRecoJetPt);
        fJetTree->SetBranchAddress("trackMax", &fRecoJetTrackMax);
        fJetTree->SetBranchAddress("jteta", &fRecoJetEta);
        fJetTree->SetBranchAddress("jtphi", &fRecoJetPhi);
        fJetTree->SetBranchAddress("WTAeta", &fRecoJetWTAEta);
        fJetTree->SetBranchAddress("WTAphi", &fRecoJetWTAPhi);

        // Gen jet quantities
        if (fIsMc)
        {
            fJetTree->SetBranchStatus("ngen", 1);
            fJetTree->SetBranchStatus("genpt", 1);
            fJetTree->SetBranchStatus("geneta", 1);
            fJetTree->SetBranchStatus("genphi", 1);
            fJetTree->SetBranchStatus("WTAgeneta", 1);
            fJetTree->SetBranchStatus("WTAgenphi", 1);
            fJetTree->SetBranchAddress("ngen", &fNGenJets);
            fJetTree->SetBranchAddress("genpt", &fGenJetPt);
            fJetTree->SetBranchAddress("geneta", &fGenJetEta);
            fJetTree->SetBranchAddress("genphi", &fGenJetPhi);
            fJetTree->SetBranchAddress("WTAgeneta", &fGenJetWTAEta);
            fJetTree->SetBranchAddress("WTAgenphi", &fGenJetWTAPhi);

            // Matching Jets
            fJetTree->SetBranchStatus("refpt", 1);
            fJetTree->SetBranchStatus("refeta", 1);
            fJetTree->SetBranchStatus("refphi", 1);
            fJetTree->SetBranchStatus("refWTAeta", 1);
            fJetTree->SetBranchStatus("refWTAphi", 1);
            fJetTree->SetBranchStatus("refparton_flavor", 1);
            fJetTree->SetBranchStatus("refparton_flavorForB", 1);
            fJetTree->SetBranchAddress("refpt", &fRefJetPt);
            fJetTree->SetBranchAddress("refeta", &fRefJetEta);
            fJetTree->SetBranchAddress("refphi", &fRefJetPhi);
            fJetTree->SetBranchAddress("refWTAeta", &fRefJetWTAEta);
            fJetTree->SetBranchAddress("refWTAphi", &fRefJetWTAPhi);
            fJetTree->SetBranchAddress("refparton_flavor", &fRefJetPartonFlavor);
            fJetTree->SetBranchAddress("refparton_flavorForB", &fRefJetPartonFlavorForB);
        }

    } // if ( fUsePartFlowJetBranch )

    // Track quantities
    if (fUseTrackBranch)
    {
        fTrkTree->SetBranchStatus("nTrk", 1);
        fTrkTree->SetBranchStatus("trkPt", 1);
        fTrkTree->SetBranchStatus("trkEta", 1);
        fTrkTree->SetBranchStatus("trkPhi", 1);
        fTrkTree->SetBranchStatus("trkPtError", 1);
        fTrkTree->SetBranchStatus("trkDxy1", 1);
        fTrkTree->SetBranchStatus("trkDxyError1", 1);
        fTrkTree->SetBranchStatus("trkDz1", 1);
        fTrkTree->SetBranchStatus("trkDzError1", 1);
        fTrkTree->SetBranchStatus("trkChi2", 1);
        fTrkTree->SetBranchStatus("trkNdof", 1);
        fTrkTree->SetBranchStatus("trkCharge", 1);
        fTrkTree->SetBranchStatus("trыkNHit", 1);
        fTrkTree->SetBranchStatus("trkNlayer", 1);
        fTrkTree->SetBranchStatus("highPurity", 1);
        fTrkTree->SetBranchStatus("pfEcal", 1);
        fTrkTree->SetBranchStatus("pfHcal", 1);
        fTrkTree->SetBranchStatus("trkMVA", 1);
        fTrkTree->SetBranchStatus("trkAlgo", 1);

        fTrkTree->SetBranchAddress("nTrk", &fNTracks);
        fTrkTree->SetBranchAddress("trkPt", &fTrackPt);
        fTrkTree->SetBranchAddress("trkEta", &fTrackEta);
        fTrkTree->SetBranchAddress("trkPhi", &fTrackPhi);
        fTrkTree->SetBranchAddress("trkPtError", &fTrackPtErr);
        fTrkTree->SetBranchAddress("trkDxy1", &fTrackDcaXY);
        fTrkTree->SetBranchAddress("trkDxyError1", &fTrackDcaXYErr);
        fTrkTree->SetBranchAddress("trkDz1", &fTrackDcaZ);
        fTrkTree->SetBranchAddress("trkDzError1", &fTrackDcaZErr);
        fTrkTree->SetBranchAddress("trkChi2", &fTrackChi2);
        fTrkTree->SetBranchAddress("trkNdof", &fTrackNDOF);
        fTrkTree->SetBranchAddress("trkCharge", &fTrackCharge);
        fTrkTree->SetBranchAddress("trkNHit", &fTrackNHits);
        fTrkTree->SetBranchAddress("trkNlayer", &fTrackNLayers);
        fTrkTree->SetBranchAddress("highPurity", &fTrackHighPurity);
        fTrkTree->SetBranchAddress("pfEcal", &fTrackPartFlowEcal);
        fTrkTree->SetBranchAddress("pfHcal", &fTrackPartFlowHcal);
        fTrkTree->SetBranchAddress("trkMVA", &fTrackMVA);
        fTrkTree->SetBranchAddress("trkAlgo", &fTrackAlgo);
    } // if ( fUseTrackBranch )

    // Gen particle quantities
    if (fIsMc && fUseGenTrackBranch)
    {
        fGenTrkTree->SetBranchStatus("pt", 1);
        fGenTrkTree->SetBranchStatus("eta", 1);
        fGenTrkTree->SetBranchStatus("phi", 1);
        fGenTrkTree->SetBranchStatus("chg", 1);
        fGenTrkTree->SetBranchStatus("pdg", 1);
        fGenTrkTree->SetBranchStatus("sube", 1);

        fGenTrkTree->SetBranchAddress("pt", &fGenTrackPt);
        fGenTrkTree->SetBranchAddress("eta", &fGenTrackEta);
        fGenTrkTree->SetBranchAddress("phi", &fGenTrackPhi);
        fGenTrkTree->SetBranchAddress("chg", &fGenTrackCharge);
        fGenTrkTree->SetBranchAddress("pdg", &fGenTrackPid);
        fGenTrkTree->SetBranchAddress("sube", &fGenTrackSube);
    }
}

//_________________
void ForestAODReader::report()
{
    std::cout << "ForestAODReader::reporting" << std::endl;
}

//_________________
void ForestAODReader::readEvent()
{

    if (fIsMc)
    {
        fRecoJet2GenJetId.clear();
        fGenJet2RecoJet.clear();
    }

    // Or one can call the clearVariables() function (will take more time)
    if (fUseGenTrackBranch && fIsMc)
    {
        fGenTrackPt.clear();
        fGenTrackEta.clear();
        fGenTrackPhi.clear();
        fGenTrackCharge.clear();
        fGenTrackPid.clear();
        fGenTrackSube.clear();
    }

    if (fEventsProcessed >= fEvents2Read)
    {
        std::cerr << "ForestAODReader::readEvent() out of entry numbers\n";
        fReaderStatus = 2; // End of input stream
    }
    fEventTree->GetEntry(fEventsProcessed);
    if (fUseHltBranch)
        fHltTree->GetEntry(fEventsProcessed);
    if (fUseSkimmingBranch)
        fSkimTree->GetEntry(fEventsProcessed);
    if (fUseJets)
        fJetTree->GetEntry(fEventsProcessed);
    if (fUseTrackBranch)
        fTrkTree->GetEntry(fEventsProcessed);
    if (fUseGenTrackBranch)
        fGenTrkTree->GetEntry(fEventsProcessed);
    fEventsProcessed++;

    // std::cout << "Events processed: " << fEventsProcessed << std::endl;
}

//________________
void ForestAODReader::fixIndices()
{

    if (fUseJets)
    {

        // Loop over reconstructed jets
        if (fIsMc)
        {
            for (Int_t iRecoJet{0}; iRecoJet < fNRecoJets; iRecoJet++)
            {

                if (fNGenJets <= 0)
                {
                    fRecoJet2GenJetId.push_back(-1);
                    continue;
                }

                // Must have a gen-matched jet
                if (fRefJetPt[iRecoJet] < 0)
                {
                    fRecoJet2GenJetId.push_back(-1);
                    continue;
                }

                for (Int_t iGenJet{0}; iGenJet < fNGenJets; iGenJet++)
                {
                    // Skip Ref and Gen jets that do not match on pT within computational precision
                    // std::cout << "|Gen pT - Ref pT| = " << TMath::Abs(fPFGenJetPt[iGenJet] - fPFRefJetPt[iRecoJet]) <<
                    // std::endl;

                    // std::cout << "iGen: " << iGenJet << " genPt: " << fPFGenJetPt[iGenJet] << std::endl;
                    if (TMath::Abs(fGenJetPt[iGenJet] - fRefJetPt[iRecoJet]) > 2.f * FLT_EPSILON)
                    {
                        // If it is the last one
                        if (iGenJet == (fNGenJets - 1))
                            fRecoJet2GenJetId.push_back(-1);
                        continue;
                    }
                    else
                    {
                        fRefJetEta[iRecoJet] = fGenJetEta[iGenJet];
                        fRefJetPhi[iRecoJet] = fGenJetPhi[iGenJet];
                        fRefJetWTAEta[iRecoJet] = fGenJetWTAEta[iGenJet];
                        fRefJetWTAPhi[iRecoJet] = fGenJetWTAPhi[iGenJet];
                        fRecoJet2GenJetId.push_back(iGenJet);
                        break;
                    }
                }
            } // for (Int_t iRecoJet=0; iRecoJet<fNPFRecoJets; iRecoJet++)

            // Fill the corresponding index in the reco vector and fill the gen
            for (Int_t iGenJet{0}; iGenJet < fNGenJets; iGenJet++)
            {
                std::vector<Int_t>::iterator it = std::find(fRecoJet2GenJetId.begin(), fRecoJet2GenJetId.end(), iGenJet);
                if (it != fRecoJet2GenJetId.end())
                {
                    fGenJet2RecoJet.push_back(std::distance(fRecoJet2GenJetId.begin(), it));
                }
                else
                {
                    fGenJet2RecoJet.push_back(-1);
                }
            }
        }
    }
}

//_________________
Event *ForestAODReader::returnEvent()
{

    // std::cout << "ForestAODReader::returnEvent" << std::endl;
    readEvent();

    Int_t nBadJets{0};

    if (fFixJetArrays)
    {
        fixIndices();
    }

    fEvent = new Event();

    fEvent->setRunId(fRunId);
    fEvent->setEventId(fEventId);
    fEvent->setLumi(fLumi);
    fEvent->setVz(fVertexZ);
    fEvent->setHiBin(fHiBin);
    if (fIsMc)
    {
        fEvent->setPtHat(fPtHat);
        fEvent->setPtHatWeight(fPtHatWeight);
    }

    // Fill HLT branch
    if (fUseHltBranch)
    {
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet15_v1(fHLT_HIAK4PFJet15_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet15_v1_Prescl(fHLT_HIAK4PFJet15_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet30_v1(fHLT_HIAK4PFJet30_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet30_v1_Prescl(fHLT_HIAK4PFJet30_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet40_v1(fHLT_HIAK4PFJet40_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet40_v1_Prescl(fHLT_HIAK4PFJet40_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet60_v1(fHLT_HIAK4PFJet60_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet60_v1_Prescl(fHLT_HIAK4PFJet60_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet80_v1(fHLT_HIAK4PFJet80_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet80_v1_Prescl(fHLT_HIAK4PFJet80_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet120_v1(fHLT_HIAK4PFJet120_v1);
        fEvent->trigAndSkim()->setHLT_HIAK4PFJet120_v1_Prescl(fHLT_HIAK4PFJet120_v1_Prescl);

        fEvent->trigAndSkim()->setHLT_HIAK8PFJet15_v1(fHLT_HIAK8PFJet15_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet15_v1_Prescl(fHLT_HIAK8PFJet15_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet25_v1(fHLT_HIAK8PFJet25_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet25_v1_Prescl(fHLT_HIAK8PFJet25_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet40_v1(fHLT_HIAK8PFJet40_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet40_v1_Prescl(fHLT_HIAK8PFJet40_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet60_v1(fHLT_HIAK8PFJet60_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet60_v1_Prescl(fHLT_HIAK8PFJet60_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet80_v1(fHLT_HIAK8PFJet80_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet80_v1_Prescl(fHLT_HIAK8PFJet80_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet140_v1(fHLT_HIAK8PFJet140_v1);
        fEvent->trigAndSkim()->setHLT_HIAK8PFJet140_v1_Prescl(fHLT_HIAK8PFJet140_v1_Prescl);

        fEvent->trigAndSkim()->setHLT_HIPFJet25_v1(fHLT_HIPFJet25_v1);
        fEvent->trigAndSkim()->setHLT_HIPFJet25_v1_Prescl(fHLT_HIPFJet25_v1_Prescl);
        fEvent->trigAndSkim()->setHLT_HIPFJet140_v1(fHLT_HIPFJet140_v1);
        fEvent->trigAndSkim()->setHLT_HIPFJet140_v1_Prescl(fHLT_HIPFJet140_v1_Prescl);

        fEvent->trigAndSkim()->setHLT_HIPuAK4CaloJet80Eta5p1_v1(fHLT_HIPuAK4CaloJet80Eta5p1_v1);
        fEvent->trigAndSkim()->setHLT_HIPuAK4CaloJet100Eta5p1_v1(fHLT_HIPuAK4CaloJet100Eta5p1_v1);
    }

    // Fill skimming branch
    if (fUseSkimmingBranch)
    {
        fEvent->trigAndSkim()->setHBHENoiseFilterResultRun2Loose(fHBHENoiseFilterResultRun2Loose);
        fEvent->trigAndSkim()->setHBHENoiseFilterResultRun2Tight(fHBHENoiseFilterResultRun2Tight);
        fEvent->trigAndSkim()->setHBHEIsoNoiseFilterResult(fHBHEIsoNoiseFilterResult);
        fEvent->trigAndSkim()->setCollisionEventSelectionAODv2(fCollisionEventSelectionAODv2);
        fEvent->trigAndSkim()->setPhfCoincFilter2Th4(fPhfCoincFilter2Th4);
        fEvent->trigAndSkim()->setPPAprimaryVertexFilter(fPPAprimaryVertexFilter);
        fEvent->trigAndSkim()->setPBeamScrapingFilter(fPBeamScrapingFilter);
        fEvent->trigAndSkim()->setPprimaryVertexFilter(fPprimaryVertexFilter);
        fEvent->trigAndSkim()->setPVertexFilterCutG(fPVertexFilterCutG);
        fEvent->trigAndSkim()->setPVertexFilterCutGloose(fPVertexFilterCutGloose);
        fEvent->trigAndSkim()->setPVertexFilterCutGtight(fPVertexFilterCutGtight);
        fEvent->trigAndSkim()->setPVertexFilterCutE(fPVertexFilterCutE);
        fEvent->trigAndSkim()->setPVertexFilterCutEandG(fPVertexFilterCutEandG);
        fEvent->trigAndSkim()->setPClusterCompatibilityFilter(fPClusterCompatibilityFilter);
    }

    // fEvent->print();

    // Create particle flow jet instances
    if (fUseJets)
    {

        // Loop over generated jets
        if (fIsMc && !fEvent->isGenJetCollectionFilled())
        {

            // std::cout << "nGenPFJets: " << fNPFGenJets << std::endl;

            for (Int_t iGenJet{0}; iGenJet < fNGenJets; iGenJet++)
            {
                GenJet *jet = new GenJet{};
                jet->setPt(fGenJetPt[iGenJet]);
                jet->setEta(fGenJetEta[iGenJet]);
                jet->setPhi(fGenJetPhi[iGenJet]);
                jet->setWTAEta(fGenJetWTAEta[iGenJet]);
                jet->setWTAPhi(fGenJetWTAPhi[iGenJet]);
                jet->setFlavor(fRefJetPartonFlavor[fGenJet2RecoJet.at(iGenJet)]);
                jet->setFlavorForB(fRefJetPartonFlavorForB[fGenJet2RecoJet.at(iGenJet)]);
                jet->setPtWeight(jetPtWeight(fIsMc, fCollidingSystem.Data(), fYearOfDataTaking, fCollidingEnergyGeV,
                                             fGenJetPt[iGenJet]));
                fEvent->genJetCollection()->push_back(jet);
            } // for (Int_t iGenJet{0}; iGenJet<fNPFGenJets; iGenJet++)

            // Projection from filling the collection several times
            fEvent->setGenJetCollectionIsFilled();
        } // if ( fIsMc )

        // Loop over reconstructed jets

        // std::cout << "nRecoPFJets: " << fNPFRecoJets << std::endl;

        for (Int_t iJet{0}; iJet < fNRecoJets; iJet++)
        {

            // Create a new jet instance
            RecoJet *jet = new RecoJet{};

            if (fIsMc)
            {
                // Count number of reconstructed jets
                // with pT > pThat of the event (wrong )
                if (fRecoJetPt[iJet] > fPtHat)
                {
                    nBadJets++;
                }

                // Add index of the matched GenJet
                jet->setGenJetId(fRecoJet2GenJetId.at(iJet));
            } // if ( fIsMc )

            // Reco
            jet->setPt(fRawJetPt[iJet]);
            jet->setEta(fRecoJetEta[iJet]);
            jet->setPhi(fRecoJetPhi[iJet]);
            jet->setWTAEta(fRecoJetWTAEta[iJet]);
            jet->setWTAPhi(fRecoJetWTAPhi[iJet]);
            if (fJEC)
            {
                fJEC->SetJetPT(fRawJetPt[iJet]);
                fJEC->SetJetEta(fRecoJetEta[iJet]);
                fJEC->SetJetPhi(fRecoJetPhi[iJet]);
                double pTcorr = fJEC->GetCorrectedPT();
                // std::cout << "pTCorr: " << pTcorr << std::endl;
                jet->setPtJECCorr(fJEC->GetCorrectedPT());
                // std::cout << "pTCorr: " << jet->recoJetPtJECCorr() << std::endl;
            }
            else
            { // If no JEC available
                jet->setPtJECCorr(fRecoJetPt[iJet]);
            }

            // jet->print();

            // Check fronе-loaded cut
            if (fJetCut && !fJetCut->pass(jet))
            {
                delete jet;
                continue;
            }

            fEvent->recoJetCollection()->push_back(jet);
        } // for (Int_t iJet{0}; iJet<fNPFRecoJets; iJet++)
    }     // if ( fUsePartFlowJetBranch )

    if (fEventCut && !fEventCut->pass(fEvent))
    {
        delete fEvent;
        fEvent = nullptr;
    }

    return fEvent;
}