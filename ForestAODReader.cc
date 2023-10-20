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
ForestAODReader::ForestAODReader() : fEvent{nullptr}, fInFileName{nullptr}, fEvents2Read{0}, fEventsProcessed{0},
    fIsMc{kFALSE}, fUseHltBranch{kTRUE}, fUseSkimmingBranch{kTRUE}, 
    fUsePartFlowJetBranch{kTRUE}, fUseCaloJetBranch{kFALSE}, 
    fUseTrackBranch{kFALSE}, fUseGenTrackBranch{kFALSE},
    fHltTree{nullptr}, fSkimTree{nullptr}, fEventTree{nullptr}, fCaloJetTree{nullptr},
    fPartFlowJetTree{nullptr}, fTrkTree{nullptr}, fGenTrkTree{nullptr},
    fJEC{nullptr}, fJECFiles{}, fJEU{nullptr}, fJEUFiles{},
    fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020},
    fYearOfDataTaking{2018}, fDoJetPtSmearing{kFALSE}, 
    fFixJetArrays{kFALSE}, fEventCut{nullptr}, fJetCut{nullptr} {
    // Initialize many variables
    clearVariables();
}

//_________________
ForestAODReader::ForestAODReader(const Char_t* inputStream,  
                    const Bool_t& useHltBranch, const Bool_t& useSkimmingBranch, 
                    const Bool_t& usePFJetBranch, const Bool_t& useCaloJetBranch, 
                    const Bool_t& useTrackBranch, const Bool_t& useGenTrackBranch, 
                    const Bool_t& isMc) : 
    fEvent{nullptr}, fInFileName{inputStream}, fEvents2Read{0}, 
    fEventsProcessed{0}, fIsMc{isMc}, 
    fUseHltBranch{useHltBranch}, fUseSkimmingBranch{useSkimmingBranch}, 
    fUsePartFlowJetBranch{usePFJetBranch}, fUseCaloJetBranch{useCaloJetBranch}, 
    fUseTrackBranch{useTrackBranch}, fUseGenTrackBranch{useGenTrackBranch},
    fJEC{nullptr}, fJECFiles{}, fJEU{nullptr}, fJEUFiles{},
    fCollidingSystem{Form("PbPb")}, fCollidingEnergyGeV{5020},
    fYearOfDataTaking{2018}, fDoJetPtSmearing{kFALSE}, 
    fFixJetArrays{kFALSE}, fEventCut{nullptr}, fJetCut{nullptr} {
    // Initialize many variables
    clearVariables();
}

//_________________
ForestAODReader::~ForestAODReader() {
    if (fEvent) delete fEvent;
    if (fHltTree) delete fHltTree;
    if (fSkimTree) delete fSkimTree;
    if (fEventTree) delete fEventTree;
    if (fCaloJetTree) delete fCaloJetTree;
    if (fPartFlowJetTree) delete fPartFlowJetTree;
    if (fTrkTree) delete fTrkTree;
    if (fGenTrkTree) delete fGenTrkTree;
    if (fJEC) delete fJEC;
    if (fJEU) delete fJEU;
    if (fEventCut) delete fEventCut;
    if (fJetCut) delete fJetCut;
}

//_________________
void ForestAODReader::clearVariables() {
    fRunId = {0};
    fEventId = {0};
    fLumi = {0};
    fVertexZ = {-999.f};
    fHiBin = {-1};
    fPtHatWeight = {-1.f};
    fPtHat = {-1.f};

    // bad jets and multiplicity to be added

    fNPFRecoJets = {0};
    fNPFGenJets = {0};
    fNCaloRecoJets = {0};
    fNCaloGenJets = {0};
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

    // Loop over jets and tracks
    for (Short_t i{0}; i<9999; i++) {

        // Jet variables
        if (i<100) {
            fPFRecoJetPt[i] = {0.f};
            fPFRecoJetEta[i] = {0.f};
            fPFRecoJetPhi[i] = {0.f};
            fPFRecoJetWTAEta[i] = {0.f};
            fPFRecoJetWTAPhi[i] = {0.f};
            fPFRecoJetTrackMax[i] = {0.f};
            fPFRefJetPt[i] = {0.f};
            fPFRefJetEta[i] = {0.f};
            fPFRefJetPhi[i] = {0.f};
            fPFRefJetWTAEta[i] = {0.f};
            fPFRefJetWTAPhi[i] = {0.f};
            fPFRefJetPartonFlavor[i] = {-999};
            fPFRefJetPartonFlavorForB[i] = {-99};
            fPFGenJetPt[i] = {0.f};
            fPFGenJetEta[i] = {0.f};
            fPFGenJetPhi[i] = {0.f};
            fPFGenJetWTAEta[i] = {0.f};
            fPFGenJetWTAPhi[i] = {0.f};

            fCaloRecoJetPt[i] = {0.f};
            fCaloRecoJetEta[i] = {0.f};
            fCaloRecoJetPhi[i] = {0.f};
            fCaloRecoJetWTAEta[i] = {0.f};
            fCaloRecoJetWTAPhi[i] = {0.f};
            fCaloRecoJetTrackMax[i] = {0.f};
            fCaloRefJetPt[i] = {0.f};
            fCaloRefJetEta[i] = {0.f};
            fCaloRefJetPhi[i] = {0.f};
            fCaloRefJetWTAEta[i] = {0.f};
            fCaloRefJetWTAPhi[i] = {0.f};
            fCaloRefJetPartonFlavor[i] = {-999};
            fCaloRefJetPartonFlavorForB[i] = {-99};
            fCaloGenJetPt[i] = {0.f};
            fCaloGenJetEta[i] = {0.f};
            fCaloGenJetPhi[i] = {0.f};
            fCaloGenJetWTAEta[i] = {0.f};
            fCaloGenJetWTAPhi[i] = {0.f};
        } // if (i<100)

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
}

//_________________
Int_t ForestAODReader::init() {
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
void ForestAODReader::setupJEC() {
    if ( fJECFiles.empty() ) {
        if (fJECInputFileName.Length()<=0) {
            setJECFileName();
        }
        fJECFiles.push_back(Form("/Users/gnigmat/work/cms/soft/jetAnalysis/aux_files/%s_%i/JEC/%s",fCollidingSystem.Data(),fCollidingEnergyGeV,fJECInputFileName.Data()));
        std::cout << Form("Add JEC file: %s\n", fJECFiles.back().c_str());
    }
	
	fJEC = new JetCorrector( fJECFiles );
}

//________________
Float_t ForestAODReader::jetPtWeight(const Bool_t &isMC, const std::string &system, const Int_t &year, 
                                     const Int_t &energy, float jetpt) const {
    Float_t jetptweight = 1.0;

    // JetPtWeightFunction is derived from MC vs data jet pT spectra.
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *JetPtWeightFunction = new TF1("JetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from all jets above 120 GeV and JECv6
        JetPtWeightFunction->SetParameters(0.79572,0.0021861,-6.35407e-06,6.66435e-09);
        jetptweight = JetPtWeightFunction->Eval(jetpt);
    }
    */
    return jetptweight;
}

//________________
Float_t ForestAODReader::leadJetPtWeight(const Bool_t &isMC, const std::string &system, const Int_t& year, 
                                         const Int_t& energy, const Float_t& leadjetpt) const {
    Float_t leadjetptweight = 1.0;

    // LeadJetPtWeightFunction is derived from MC vs data leading jet pT spectra.
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *LeadJetPtWeightFunction = new TF1("LeadJetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from leading jets above 120 GeV and JECv6
        LeadJetPtWeightFunction->SetParameters(0.876682,0.00131479,-3.90884e-06,4.40358e-09); ;
        leadjetptweight = LeadJetPtWeightFunction->Eval(leadjetpt);
    }
    */
    return leadjetptweight;
}

//________________
Float_t ForestAODReader::subleadJetPtWeight(const Bool_t &isMC, const std::string &system, const Int_t &year, 
                                            const Int_t &energy, const Float_t &subleadjetpt) {
    Float_t subleadjetptweight = 1.0;

    // SubLeadJetPtWeightFunction is derived from MC vs data subleading jet pT spectra.
    /*
    if(isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *SubLeadJetPtWeightFunction = new TF1("SubLeadJetPtWeightFunction", "pol3", 0.0, 500.0); //Derived from leading jets above 120 GeV and JECv6
        SubLeadJetPtWeightFunction->SetParameters(0.876682,0.00131479,-3.90884e-06,4.40358e-09); ;
        subleadjetptweight = SubLeadJetPtWeightFunction->Eval(subleadjetpt);
    }
    */
    return subleadjetptweight;
}

//________________
Float_t ForestAODReader::jetPtSmeringWeight(const Bool_t &isMC, const std::string &system, const Int_t &year, 
                                            const Int_t &energy, const Float_t &jetpt, 
                                            const Bool_t &dosmearing, const Float_t resolutionfactor) const {
    Float_t jetptsmearweight = 1.0;
    if(!dosmearing) return jetptsmearweight;

    // JetPtSmearingWeightFunction is derived from MC vs data jet pT spectra.
    if(!isMC && system == "pp" && energy == 5020 && year == 2017){
        TF1 *JetPtSmearingWeightFunction = new TF1("JetPtSmearingWeightFunction", "pol3", 0.0, 500.0); //Derived from all jets above 120 GeV and JECv6
        JetPtSmearingWeightFunction->SetParameters(0.174881, -0.00091979, 3.50064e-06, -6.52541e-09, 4.64199e-12);
        jetptsmearweight = JetPtSmearingWeightFunction->Eval(jetpt);
        jetptsmearweight = jetptsmearweight*resolutionfactor;
    }

    return jetptsmearweight;
}

//________________
Float_t ForestAODReader::trkEtaMixWeight(const Bool_t &isMC, const std::string &system, const Int_t &year, const Int_t &energy, 
                                         const Float_t &trketa, const Bool_t &reco) const {
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
Float_t ForestAODReader::eventWeight(const Bool_t &isMC, const Bool_t &use_centrality, 
                                     const std::string& system, const Int_t &year, const Int_t &energy, 
                                     const Float_t &vz, const Int_t mult, const Float_t &weighttree, 
                                     const Float_t &leadjetpt) const {

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

	if (isMC && !use_centrality && system == "pp" && 
        energy == 5020 && year == 2017) {

		TF1 *VzWeightFunction = new TF1("VzWeightFunction", "pol6", -15.0, 15.0);
		VzWeightFunction->SetParameters(0.973805, 0.00339418, 0.000757544, -1.37331e-06, -2.82953e-07, -3.06778e-10, 3.48615e-09);
		vzweight = VzWeightFunction->Eval(vz);

		TF1 *MultCentWeightFunction = new TF1("MultCentWeightFunction", "pol0", 0.0, 500.0);
		MultCentWeightFunction->SetParameter(0,1.0);
		multweight = MultCentWeightFunction->Eval(mult);

		TF1 *MultTriggerWeightFunction = new TF1("MultTriggerWeightFunction", "pol0", 0.0, 500.0); // fitted from turn on curves
		MultTriggerWeightFunction->SetParameter(0,1.0);
		Float_t multtrigweight = 1.0;
		multtrigweight = MultTriggerWeightFunction->Eval(mult);
		multefficiency = 1./multtrigweight;

		TF1 *JetTriggerWeightFunction = new TF1("JetTriggerWeightFunction", "pol0", 0.0, 500.0); // fitted from turn on curves
		JetTriggerWeightFunction->SetParameter(0,1.0);
		Float_t jettrigweight = 1.0;
		jettrigweight = JetTriggerWeightFunction->Eval(leadjetpt);
		jetefficiency = 1./jettrigweight;

		evtweight = weighttree;
	}

	totalweight = evtweight * multweight * vzweight * multefficiency * jetefficiency;
	return totalweight;
}

//_________________
void ForestAODReader::finish() {

}

//_________________
Int_t ForestAODReader::setupChains() {

    // Setup chains (0-good, 1-bad)
    Int_t returnStatus = 1;

    // Setup chains to read

    std::cout << "Setting chains... ";

    // Use event branch
    fEventTree = new TChain("hiEvtAnalyzer/HiTree");

    // Use HLT branch
    if ( fUseHltBranch ) {
        fHltTree = new TChain("hltanalysis/HltTree");
    }
    // Use skimming branch
    if ( fUseSkimmingBranch ) {
        fSkimTree = new TChain("skimanalysis/HltTree");
    }
    // Use calo jet branch
    if ( fUseCaloJetBranch ) {
        fCaloJetTree = new TChain("akPu4CaloJetAnalyzer/t");
    }
    // Use particle flow jet branch
    if ( fUsePartFlowJetBranch ) {
        fPartFlowJetTree = new TChain("akCs4PFJetAnalyzer/t");
    }
    // Use reconstructed track branch
    if ( fUseTrackBranch ) {
        fTrkTree = new TChain("ppTrack/trackTree");
    }
    // Use generated track branch
    if ( fIsMc && fUseGenTrackBranch ) {
        fGenTrkTree = new TChain("HiGenParticleAna/hi");
    }
    std::cout << "\t[DONE]\n";

    // Initialize input file name (should switch to const char* processing later)
    TString input(fInFileName);

    // Check input exists
    if (  input.Length()<= 0 ) {
        std::cerr << "No normal inputfile. Terminating." << std::endl;
        returnStatus = 1;
        exit(0);
    }
    // Normail input
    else {
        // If input is a single ROOT file
        if ( input.Index(".root") > 0 ) {
            std::cout << Form( "Adding %s file to chains\n", input.Data() );
            fEventTree->Add( input.Data() );
            if ( fUseHltBranch ) fHltTree->Add( input.Data() );
            if ( fUseSkimmingBranch ) fSkimTree->Add( input.Data() );
            if ( fUseCaloJetBranch ) fCaloJetTree->Add(input.Data() );
            if ( fUsePartFlowJetBranch ) fPartFlowJetTree->Add( input.Data() );
            if ( fUseTrackBranch ) fTrkTree->Add( input.Data() );
            if ( fIsMc && fUseGenTrackBranch ) fGenTrkTree->Add( input.Data() );

            fEvents2Read = fEventTree->GetEntries();
            std::cout << Form("Total number of events to read: %lld\n", fEvents2Read );
            Long64_t fEvents2Read2 = fPartFlowJetTree->GetEntries();
            std::cout << Form("Total number of events to read2: %lld\n", fEvents2Read2 );
        }
        // Assuming that list of files is provided instead of a single file
        else {
            std::ifstream inputStream( input.Data() );

            if ( !inputStream ) std::cout << Form( "ERROR: Cannot open file list: %s\n", input.Data() );
            Int_t nFiles = 0;
            std::string file;
            size_t pos;
            while ( getline( inputStream, file ) ) {
                // NOTE: our external formatters may pass "file NumEvents"
                //       Take only the first part
                //cout << "DEBUG found " <<  file << endl;
                pos = file.find_first_of(" ");
                if ( pos != std::string::npos ) file.erase( pos, file.length() - pos );
                //cout << "DEBUG found [" <<  file << "]" << endl;

                // Check that file is of a correct name
                if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos &&
                     file.find("AOD") != std::string::npos ) {
                    
                    // Open file
                    TFile* ftmp = TFile::Open(file.c_str());

                    // Check file is not zombie and contains information
                    if ( ftmp && !ftmp->IsZombie() && ftmp->GetNkeys() ) {
                        std::cout << Form("Adding file to chain: %s\n", file.c_str() );
                        // Adding file to chains
                        fEventTree->Add( file.c_str() );
                        if ( fUseHltBranch ) fHltTree->Add( file.c_str() );
                        if ( fUseSkimmingBranch ) fSkimTree->Add( file.c_str() );
                        if ( fUseCaloJetBranch ) fCaloJetTree->Add( file.c_str() );
                        if ( fUsePartFlowJetBranch ) fPartFlowJetTree->Add( file.c_str() );
                        if ( fUseTrackBranch ) fTrkTree->Add( file.c_str() );
                        if ( fIsMc && fUseGenTrackBranch ) fGenTrkTree->Add( file.c_str() );
                        ++nFiles;
                    } //if(ftmp && !ftmp->IsZombie() && ftmp->GetNkeys())

                    if (ftmp) {
                        ftmp->Close();
                    } //if (ftmp)
                } //if ( file.find(".root") != std::string::npos && file.find("Forest") != std::string::npos && file.find("AOD") != std::string::npos )
            } //while ( getline( inputStream, file ) )

            std::cout << Form("Total number of files in chain: %d\n", nFiles);
            fEvents2Read = fEventTree->GetEntries();
            std::cout << Form("Total number of events to read: %lld\n", fEvents2Read );
        } // else {   if file list
        returnStatus = 0;
    } // else {   if normal input
    return returnStatus;
}

//_________________
void ForestAODReader::setupBranches() {

    // Disable all branches - this is important while reading big files
    fEventTree->SetBranchStatus("*", 0);
    if ( fUseHltBranch ) fHltTree->SetBranchStatus("*", 0);
    if ( fUseSkimmingBranch ) fSkimTree->SetBranchStatus("*", 0);
    if ( fUseCaloJetBranch) fCaloJetTree->SetBranchStatus("*", 0);
    if ( fUsePartFlowJetBranch ) fPartFlowJetTree->SetBranchStatus("*", 0);
    if ( fUseTrackBranch ) fTrkTree->SetBranchStatus("*", 0);
    if ( fUseGenTrackBranch && fIsMc ) fGenTrkTree->SetBranchStatus("*", 0);


    // enable branches of interest -> see definition of each variables above

    // Event quantities
    fEventTree->SetBranchStatus("run", 1);
    fEventTree->SetBranchStatus("evt", 1);
    fEventTree->SetBranchStatus("lumi", 1);
    fEventTree->SetBranchStatus("vz", 1);
    fEventTree->SetBranchStatus("hiBin", 1); //centrality only for PbPb and XeXe
    fEventTree->SetBranchAddress("run", &fRunId);
    fEventTree->SetBranchAddress("evt", &fEventId);
    fEventTree->SetBranchAddress("lumi", &fLumi);
    fEventTree->SetBranchAddress("vz", &fVertexZ);
    fEventTree->SetBranchAddress("hiBin", &fHiBin); 

    if( fIsMc ) {
        fEventTree->SetBranchStatus("weight", 1);
        fEventTree->SetBranchStatus("pthat", 1); 
        fEventTree->SetBranchAddress("weight", &fPtHatWeight);
        fEventTree->SetBranchAddress("pthat", &fPtHat);
    }

    if ( fUseHltBranch ) {
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
    if ( fUseSkimmingBranch ) {
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
    } // if ( fUseSkimmingBranch )

    // Jet quantities
    if ( fUsePartFlowJetBranch ) {

        fPartFlowJetTree->SetBranchStatus("nref", 1);
        fPartFlowJetTree->SetBranchStatus("rawpt", 1);
        fPartFlowJetTree->SetBranchStatus("trackMax", 1);
        fPartFlowJetTree->SetBranchStatus("jteta", 1);
        fPartFlowJetTree->SetBranchStatus("jtphi", 1);
        fPartFlowJetTree->SetBranchStatus("WTAeta", 1);
        fPartFlowJetTree->SetBranchStatus("WTAphi", 1);

        fPartFlowJetTree->SetBranchAddress("nref", &fNPFRecoJets);
        fPartFlowJetTree->SetBranchAddress("rawpt", &fPFRecoJetPt);
        fPartFlowJetTree->SetBranchAddress("trackMax", &fPFRecoJetTrackMax);
        fPartFlowJetTree->SetBranchAddress("jteta", &fPFRecoJetEta);
        fPartFlowJetTree->SetBranchAddress("jtphi", &fPFRecoJetPhi);
        fPartFlowJetTree->SetBranchAddress("WTAeta", &fPFRecoJetWTAEta);
        fPartFlowJetTree->SetBranchAddress("WTAphi", &fPFRecoJetWTAPhi);

        // Gen jet quantities
        if ( fIsMc ) {
            fPartFlowJetTree->SetBranchStatus("ngen", 1);
            fPartFlowJetTree->SetBranchStatus("genpt", 1);
            fPartFlowJetTree->SetBranchStatus("geneta", 1);
            fPartFlowJetTree->SetBranchStatus("genphi", 1);
            fPartFlowJetTree->SetBranchStatus("WTAgeneta", 1);
            fPartFlowJetTree->SetBranchStatus("WTAgenphi", 1);
            fPartFlowJetTree->SetBranchAddress("ngen", &fNPFGenJets);
            fPartFlowJetTree->SetBranchAddress("genpt", &fPFGenJetPt);
            fPartFlowJetTree->SetBranchAddress("geneta", &fPFGenJetEta);
            fPartFlowJetTree->SetBranchAddress("genphi", &fPFGenJetPhi);
            fPartFlowJetTree->SetBranchAddress("WTAgeneta", &fPFGenJetWTAEta);
            fPartFlowJetTree->SetBranchAddress("WTAgenphi", &fPFGenJetWTAPhi);
        }

        // Jet-matching quantities
        if ( fIsMc ) {
            fPartFlowJetTree->SetBranchStatus("refpt", 1);
            fPartFlowJetTree->SetBranchStatus("refeta", 1);
            fPartFlowJetTree->SetBranchStatus("refphi", 1);
            fPartFlowJetTree->SetBranchStatus("refWTAeta", 1);
            fPartFlowJetTree->SetBranchStatus("refWTAphi", 1);
            fPartFlowJetTree->SetBranchStatus("refparton_flavor", 1);
            fPartFlowJetTree->SetBranchStatus("refparton_flavorForB", 1);
            fPartFlowJetTree->SetBranchAddress("refpt", &fPFRefJetPt);
            fPartFlowJetTree->SetBranchAddress("refeta", &fPFRefJetEta);
            fPartFlowJetTree->SetBranchAddress("refphi", &fPFRefJetPhi);
            fPartFlowJetTree->SetBranchAddress("refWTAeta", &fPFRefJetWTAEta);
            fPartFlowJetTree->SetBranchAddress("refWTAphi", &fPFRefJetWTAPhi);
            fPartFlowJetTree->SetBranchAddress("refparton_flavor", &fPFRefJetPartonFlavor);
            fPartFlowJetTree->SetBranchAddress("refparton_flavorForB", &fPFRefJetPartonFlavorForB);
        }

    } // if ( fUsePartFlowJetBranch )

    if ( fUseCaloJetBranch ) {
        fCaloJetTree->SetBranchStatus("nref", 1);
        fCaloJetTree->SetBranchStatus("rawpt", 1);
        fCaloJetTree->SetBranchStatus("trackMax", 1);
        fCaloJetTree->SetBranchStatus("jteta", 1);
        fCaloJetTree->SetBranchStatus("jtphi", 1);
        fCaloJetTree->SetBranchStatus("WTAeta", 1);
        fCaloJetTree->SetBranchStatus("WTAphi", 1);

        fCaloJetTree->SetBranchAddress("nref", &fNCaloRecoJets);
        fCaloJetTree->SetBranchAddress("rawpt", &fCaloRecoJetPt);
        fCaloJetTree->SetBranchAddress("trackMax", &fCaloRecoJetTrackMax);
        fCaloJetTree->SetBranchAddress("jteta", &fCaloRecoJetEta);
        fCaloJetTree->SetBranchAddress("jtphi", &fCaloRecoJetPhi);
        fCaloJetTree->SetBranchAddress("WTAeta", &fCaloRecoJetWTAEta);
        fCaloJetTree->SetBranchAddress("WTAphi", &fCaloRecoJetWTAPhi);

        // Gen jet quantities
        if ( fIsMc ) {
            fCaloJetTree->SetBranchStatus("ngen", 1);
            fCaloJetTree->SetBranchStatus("genpt", 1);
            fCaloJetTree->SetBranchStatus("geneta", 1);
            fCaloJetTree->SetBranchStatus("genphi", 1);
            fCaloJetTree->SetBranchStatus("WTAgeneta", 1);
            fCaloJetTree->SetBranchStatus("WTAgenphi", 1);

            fCaloJetTree->SetBranchAddress("ngen", &fNCaloGenJets);
            fCaloJetTree->SetBranchAddress("genpt", &fCaloGenJetPt);
            fCaloJetTree->SetBranchAddress("geneta", &fCaloGenJetEta);
            fCaloJetTree->SetBranchAddress("genphi", &fCaloGenJetPhi);
            fCaloJetTree->SetBranchAddress("WTAgeneta", &fCaloGenJetWTAEta);
            fCaloJetTree->SetBranchAddress("WTAgenphi", &fCaloGenJetWTAPhi);
        }

        // Jet-matching quantities
        if ( fIsMc ) {
            fCaloJetTree->SetBranchStatus("refpt", 1);
            fCaloJetTree->SetBranchStatus("refeta", 1);
            fCaloJetTree->SetBranchStatus("refphi", 1);
            fCaloJetTree->SetBranchAddress("refWTAeta", &fCaloRefJetWTAEta);
            fCaloJetTree->SetBranchAddress("refWTAphi", &fCaloRefJetWTAPhi);
            fCaloJetTree->SetBranchStatus("refparton_flavor", 1);
            fCaloJetTree->SetBranchStatus("refparton_flavorForB", 1);
            fCaloJetTree->SetBranchAddress("refpt", &fCaloRefJetPt);
            fCaloJetTree->SetBranchAddress("refeta", &fCaloRefJetEta);
            fCaloJetTree->SetBranchAddress("refphi", &fCaloRefJetPhi);
            fCaloJetTree->SetBranchAddress("refWTAeta", &fCaloRefJetWTAEta);
            fCaloJetTree->SetBranchAddress("refWTAphi", &fCaloRefJetWTAPhi);
            fCaloJetTree->SetBranchAddress("refparton_flavor", &fCaloRefJetPartonFlavor);
            fCaloJetTree->SetBranchAddress("refparton_flavorForB", &fCaloRefJetPartonFlavorForB);
        }
    } // if ( fUsePartFlowJetBranch )

    // Track quantities
    if ( fUseTrackBranch ) {
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
    if( fIsMc && fUseGenTrackBranch ) {
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
void ForestAODReader::report() {
    std::cout << "ForestAODReader::reporting" << std::endl;
}

//_________________
void ForestAODReader::readEvent() {

    // Or one can call the clearVariables() function (will take more time)
    if (fUseGenTrackBranch) {
        fGenTrackPt.clear();
        fGenTrackEta.clear();
        fGenTrackPhi.clear();
        fGenTrackCharge.clear();
        fGenTrackPid.clear();
        fGenTrackSube.clear();
    }

    if ( fEventsProcessed >= fEvents2Read ) { 
        std::cerr << "ForestAODReader::readEvent() out of entry numbers\n"; 
        fReaderStatus = 2; // End of input stream
    }
    fEventTree->GetEntry( fEventsProcessed );
    if (fUseHltBranch) fHltTree->GetEntry(fEventsProcessed);
    if (fUseSkimmingBranch) fSkimTree->GetEntry(fEventsProcessed);
    if (fUseCaloJetBranch) fCaloJetTree->GetEntry(fEventsProcessed);
    if (fUsePartFlowJetBranch) fPartFlowJetTree->GetEntry(fEventsProcessed);
    if (fUseTrackBranch) fTrkTree->GetEntry(fEventsProcessed);
    if (fUseGenTrackBranch) fGenTrkTree->GetEntry(fEventsProcessed);
    fEventsProcessed++;

    //std::cout << "Events processed: " << fEventsProcessed << std::endl;
}

//________________
void ForestAODReader::fixIndices() {

    //std::cout << "Fix indices" << std::endl;
    if (fUsePartFlowJetBranch) {
        //std::cout << "Will fix PF jets. nRecoJets: " << fNPFRecoJets << " nGenJet: " << fNPFGenJets << std::endl;
        // Loop over reconstructed jets
        for (Int_t iRecoJet=0; iRecoJet<fNPFRecoJets; iRecoJet++) {
            // Must have a gen-matched jet
            if ( fPFRefJetPt[iRecoJet] < 0) continue;
            for (Int_t iGenJet=0; iGenJet<fNPFGenJets; iGenJet++) {
                // Skip Ref and Gen jets that do not match on pT within computational precision
                //std::cout << "|Gen pT - Ref pT| = " << TMath::Abs(fPFGenJetPt[iGenJet] - fPFRefJetPt[iRecoJet]) << std::endl;
                if ( TMath::Abs(fPFGenJetPt[iGenJet] - fPFRefJetPt[iRecoJet]) > 2.f * FLT_EPSILON ) continue;
                fPFRefJetEta[iRecoJet] = fPFGenJetEta[iGenJet];
                fPFRefJetPhi[iRecoJet] = fPFGenJetPhi[iGenJet];
                fPFRefJetWTAEta[iRecoJet] = fPFGenJetWTAEta[iGenJet];
                fPFRefJetWTAPhi[iRecoJet] = fPFGenJetWTAPhi[iGenJet];
            }
        } //for (Int_t iRecoJet=0; iRecoJet<fNPFRecoJets; iRecoJet++)
    } // if (fUsePartFlowJetBranch)

    if (fUseCaloJetBranch) {
        //std::cout << "Will fix Calo jets. nRecoJets: " << fNCaloRecoJets << " nGenJet: " << fNCaloGenJets << std::endl;

        // Loop over reconstructed jets
        for (Int_t iRecoJet=0; iRecoJet<fNCaloRecoJets; iRecoJet++) {
            // Must have a gen-matched jet
            if ( fCaloRefJetPt[iRecoJet] < 0) continue;
            for (Int_t iGenJet=0; iGenJet<fNCaloGenJets; iGenJet++) {
                // Skip Ref and Gen jets that do not match on pT within computational precision
                //std::cout << "|Gen pT - Ref pT| = " << TMath::Abs(fCaloGenJetPt[iGenJet] - fCaloRefJetPt[iRecoJet]) << std::endl;
                if ( TMath::Abs(fCaloGenJetPt[iGenJet] - fCaloRefJetPt[iRecoJet]) > 2.f * FLT_EPSILON ) continue;
                fCaloRefJetEta[iRecoJet] = fCaloGenJetEta[iGenJet];
                fCaloRefJetPhi[iRecoJet] = fCaloGenJetPhi[iGenJet];
                fCaloRefJetWTAEta[iRecoJet] = fCaloGenJetWTAEta[iGenJet];
                fCaloRefJetWTAPhi[iRecoJet] = fCaloGenJetWTAPhi[iGenJet];
            }
        } //for (Int_t iRecoJet=0; iRecoJet<fNCaloRecoJets; iRecoJet++)
    } // if (fUseCaloJetBranch)

    //std::cout << "Indices fixed" << std::endl;
}

//_________________
Event* ForestAODReader::returnEvent() {

    //std::cout << "ForestAODReader::returnEvent" << std::endl;
    readEvent();

    Int_t nBadPFJets{0};
    Int_t nBadCaloJets{0};

    if (fFixJetArrays) {
        fixIndices();
    }

    fEvent = new Event();

    fEvent->setRunId( fRunId );
    fEvent->setEventId( fEventId );
    fEvent->setLumi( fLumi );
    fEvent->setVz( fVertexZ );
    fEvent->setHiBin( fHiBin );
    if ( fIsMc ) {
        fEvent->setPtHat( fPtHat );
        fEvent->setPtHatWeight( fPtHatWeight );
    }

    // Fill HLT branch
    if ( fUseHltBranch ) {
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
    if ( fUseSkimmingBranch ) {
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
    }

    //fEvent->print();
    
    // Create particle flow jet instances
    if ( fUsePartFlowJetBranch ) {

        // Loop over reconstructed jets
        
        //std::cout << "nJets: " << fNPFRecoJets << std::endl;

        for (Int_t iJet{0}; iJet<fNPFRecoJets; iJet++) {

            // Create a new jet instance
            Jet *jet = new Jet();

            if ( fIsMc ) {
                // Count number of reconstructed jets
                // with pT > pThat of the event (wrong )
                if ( fPFRecoJetPt[iJet] > fPtHat ) {
                    nBadPFJets++;
                }
            } // if ( fIsMc )

            // Reco
            jet->setRecoJetPt( fPFRecoJetPt[iJet] );
            jet->setRecoJetEta( fPFRecoJetEta[iJet] );
            jet->setRecoJetPhi( fPFRecoJetPhi[iJet] );
            if ( fJEC ) {
                fJEC->SetJetPT( fPFRecoJetPt[iJet] );
                fJEC->SetJetEta( fPFRecoJetEta[iJet] );
                fJEC->SetJetPhi( fPFRecoJetPhi[iJet] );
                double pTcorr = fJEC->GetCorrectedPT();
                //std::cout << "pTCorr: " << pTcorr << std::endl; 
                jet->setRecoJetPtJECCorr( fJEC->GetCorrectedPT() );
                //std::cout << "pTCorr: " << jet->recoJetPtJECCorr() << std::endl; 
                
                jet->setRecoJetPtWeight( jetPtWeight( fIsMc, fCollidingSystem.Data(), fYearOfDataTaking, 
                                                      fCollidingEnergyGeV, jet->recoJetPtJECCorr() ) );
                
                //resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
                jet->setRecoJetPtSmearingWeight( jetPtSmeringWeight( fIsMc, fCollidingSystem.Data(), fYearOfDataTaking, 
                                                                     fCollidingEnergyGeV, jet->recoJetPtJECCorr(), 
                                                                     fDoJetPtSmearing, 0.663) );
            }
            else { // If no JEC available
                jet->setRecoJetPtJECCorr( -999.f );
                jet->setRecoJetPtWeight( -999.f );
                jet->setRecoJetPtSmearingWeight( -999.f );
            }
            jet->setRecoJetWTAEta( fPFRecoJetWTAEta[iJet] );
            jet->setRecoJetWTAPhi( fPFRecoJetWTAPhi[iJet] );

            if ( fIsMc ) {
                // Ref
                jet->setRefJetPt( fPFRefJetPt[iJet] );
                jet->setRefJetEta( fPFRefJetEta[iJet] );
                jet->setRefJetPhi( fPFRefJetPhi[iJet] );
                jet->setRefJetWTAEta( fPFRefJetWTAEta[iJet] );
                jet->setRefJetWTAPhi( fPFRefJetWTAPhi[iJet] );
                jet->setRefFlavor( fPFRefJetPartonFlavor[iJet] );
                jet->setRefFlavorForB( fPFRefJetPartonFlavorForB[iJet] );
                jet->setRefJetPtWeight( jetPtWeight(fIsMc, fCollidingSystem.Data(), fYearOfDataTaking, 
                                                    fCollidingEnergyGeV, fPFRefJetPt[iJet]) );
            }

            // jet->print();

            // Check fronе-loaded cut
            if ( fJetCut && !fJetCut->pass(jet) ) {
                delete jet;
                continue;
            }

            fEvent->pfJetCollection()->push_back( jet );
        } // for (Int_t iJet{0}; iJet<fNPFRecoJets; iJet++)
    } // if ( fUsePartFlowJetBranch )


    // Create jet instances
    if ( fUseCaloJetBranch ) {

        // Loop over reconstructed jets
        for (Int_t iJet{0}; iJet<fNCaloRecoJets; iJet++) {

            // Create a new jet instance
            Jet *jet = new Jet();

            if ( fIsMc ) {
                // Count number of reconstructed jets
                // with pT > pThat of the event (wrong )
                if ( fCaloRecoJetPt[iJet] > fPtHat ) {
                    nBadCaloJets++;
                }
            } // if ( fIsMc )

            // Reco
            jet->setRecoJetPt( fCaloRecoJetPt[iJet] );
            jet->setRecoJetEta( fCaloRecoJetEta[iJet] );
            jet->setRecoJetPhi( fCaloRecoJetPhi[iJet] );
            jet->setRecoJetPtJECCorr( -999.f );
            jet->setRecoJetWTAEta( fCaloRecoJetWTAEta[iJet] );
            jet->setRecoJetWTAPhi( fCaloRecoJetWTAPhi[iJet] );

            if ( fJEC ) {
                fJEC->SetJetPT( fCaloRecoJetPt[iJet] );
                fJEC->SetJetEta( fCaloRecoJetEta[iJet] );
                fJEC->SetJetPhi( fCaloRecoJetPhi[iJet] );
                jet->setRecoJetPtJECCorr( fJEC->GetCorrectedPT() );
                jet->setRecoJetPtWeight( jetPtWeight( fIsMc, fCollidingSystem.Data(), fYearOfDataTaking, 
                                                      fCollidingEnergyGeV, jet->recoJetPtJECCorr() ) );
                
                //resolutionfactor: Worsening resolution by 20%: 0.663, by 10%: 0.458 , by 30%: 0.831
                jet->setRecoJetPtSmearingWeight( jetPtSmeringWeight( fIsMc, fCollidingSystem.Data(), fYearOfDataTaking, 
                                                                     fCollidingEnergyGeV, jet->recoJetPtJECCorr(), 
                                                                     fDoJetPtSmearing, 0.663) );
            }
            else { // If no JEC available
                jet->setRecoJetPtJECCorr( -999.f );
                jet->setRecoJetPtWeight( -999.f );
                jet->setRecoJetPtSmearingWeight( -999.f );
            } 

            if ( fIsMc) {
                // Ref
                jet->setRefJetPt( fCaloRefJetPt[iJet] );
                jet->setRefJetEta( fCaloRefJetEta[iJet] );
                jet->setRefJetPhi( fCaloRefJetPhi[iJet] );
                jet->setRefJetWTAEta( fCaloRefJetWTAEta[iJet] );
                jet->setRefJetWTAPhi( fCaloRefJetWTAPhi[iJet] );
                jet->setRefFlavor( fCaloRefJetPartonFlavor[iJet] );
                jet->setRefFlavorForB( fCaloRefJetPartonFlavorForB[iJet] );
            }

            // Check fron-loaded cut
            if ( fJetCut && !fJetCut->pass(jet) ) {
                delete jet;
                continue;
            }

            fEvent->caloJetCollection()->push_back( jet );
        } // for (Int_t iJet{0}; iJet<fNCaloRecoJets; iJet++)
    } // if ( fUseCaloJetBranch )

    if ( fEventCut && !fEventCut->pass(fEvent) ) {
        delete fEvent;
        fEvent = nullptr;
    }

    return fEvent;
}