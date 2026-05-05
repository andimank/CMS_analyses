// -*- C++ -*-
//
// Package:    OOAnalysisFlow/SPAnalsis
// Class:      SPAnalsis
//
/**\class SPAnalsis SPAnalsis.cc OOAnalysisFlow/SPAnalsis/plugins/SPAnalsis.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Shengquan Tuo
//         Created:  Thu, 01 May 2025 12:40:17 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/CentralityBins.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "OOAnalysisFlow/SPAnalysis/plugins/SPAnalysis.h"


using namespace std;
using namespace reco;
using namespace edm;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<>
// This will improve performance in multithreaded jobs.

using reco::TrackCollection;

class SPAnalysis : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit SPAnalysis(const edm::ParameterSet&);
  ~SPAnalysis() override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override;

  // ----------member data ---------------------------
  edm::EDGetTokenT< edm::View< pat::PackedCandidate > > tracksToken_;
  edm::EDGetTokenT< pat::PackedCandidateCollection > pfToken_;
  edm::EDGetTokenT< reco::VertexCollection > vertexToken_;
  edm::EDGetTokenT<reco::GenParticleCollection> genTracksToken_;
  //edm::EDGetTokenT< edm::View< pat::PackedCandidate > > tracksTokenMiniAOD_;
  //edm::EDGetTokenT< edm::ValueMap< float > > chi2Map_;
//  edm::EDGetTokenT< reco::Centrality > centralityTag_;
  edm::EDGetTokenT< int > centralityBinTag_;
  //edm::EDGetTokenT< std::vector<reco::EvtPlane> > eventPlane_;
  bool isPixelTracks_;

  double etaCutHF_;
  double etaCutTrk_;
  double etaCutTrkNonC_;

  edm::Service<TFileService> fs;


  TComplex QhfPlus[nCentBin][nHarmonics];
  TComplex QhfMinus[nCentBin][nHarmonics];
  TComplex QhfMinusPlus[nCentBin][nHarmonics];
  TComplex QDEMinusPlus[nCentBin][nHarmonics];
  TComplex QtrkPtPlus[nCentBin][nHarmonics][nPtBin];
  TComplex QtrkPtMinus[nCentBin][nHarmonics][nPtBin];
  TComplex QtrkMid[nCentBin][nHarmonics];
  TComplex QtrkPtMid[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnPlusQnMinusStar[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnMinusQnPlusStar[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnPlusQnMinusStar[nCentBin][nHarmonics]; //Minus:A; Plus:B; Mid:C
  TH1D* hQnMinusQnPlusStar[nCentBin][nHarmonics];
  TH1D* hQnMinusQnMidStar[nCentBin][nHarmonics];
  TH1D* hQnPlusQnMidStar[nCentBin][nHarmonics];
  TH1D* hqnQnPlusStar[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnQnMinusStar[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnQnMinusPlusStar[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnMinusPlusQnDEStar[nCentBin][nHarmonics];
  TH1D* hQnMinusPlusQnMidStar[nCentBin][nHarmonics];
  TH1D* hQnDEQnMidStar[nCentBin][nHarmonics];

  TH1D* hQnPlusQnMinusStarPtTrkPlus[nCentBin][nHarmonics][nPtBin]; //tuos
  TH1D* hQnPlusQnMinusStarPtTrkMinus[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnMinusQnPlusStarPtTrkPlus[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnMinusQnPlusStarPtTrkMinus[nCentBin][nHarmonics][nPtBin];
  TH1D* hSumCosn[nCentBin][nHarmonics][nPtBin];

//tuos-sub
  TComplex QA1[nCentBin][nHarmonics]; //minus, trk (-2.4,-1) 
  TComplex QA2[nCentBin][nHarmonics]; //minus, hcal (-2.4, -1)
  TComplex QA3[nCentBin][nHarmonics]; //minus, hf (-5.2,-2)
  TComplex QA4[nCentBin][nHarmonics]; //minus, hf (-5.2,-3)
  TComplex QA4x[nCentBin][nHarmonics]; //minus, hf (-5.2,-3)
  TComplex QA5[nCentBin][nHarmonics]; //minus, hf (-5.2,-4)
  TComplex QA6[nCentBin][nHarmonics]; //minus, hf (-5.2,-4.5)
  TComplex QA6x[nCentBin][nHarmonics]; //minus, hf (-5.2,-4.5)
  TComplex QB1[nCentBin][nHarmonics]; //plus, trk (1,2.4)
  TComplex QB2[nCentBin][nHarmonics]; //plus, hcal (1, 2.4)
  TComplex QB3[nCentBin][nHarmonics]; //plus, hf (2, 5.2)
  TComplex QB4[nCentBin][nHarmonics]; //plus, hf (3, 5.2)
  TComplex QB4x[nCentBin][nHarmonics]; //plus, hf (3, 5.2)
  TComplex QB5[nCentBin][nHarmonics]; //plus, hf (4, 5.2)
  TComplex QB6[nCentBin][nHarmonics]; //plus, hf (4.5, 5.2)
  TComplex QB6x[nCentBin][nHarmonics]; //plus, hf (4.5, 5.2)
  double sumwA1[nCentBin];
  double sumwA2[nCentBin];
  double sumwA3[nCentBin];
  double sumwA4[nCentBin];
  double sumwA4x[nCentBin];
  double sumwA5[nCentBin];
  double sumwA6[nCentBin];
  double sumwA6x[nCentBin];
  double sumwB1[nCentBin];
  double sumwB2[nCentBin];
  double sumwB3[nCentBin];
  double sumwB4[nCentBin];
  double sumwB4x[nCentBin];
  double sumwB5[nCentBin];
  double sumwB6[nCentBin];
  double sumwB6x[nCentBin];
  TH1D* hQnA1QnB1[nCentBin][nHarmonics];
  TH1D* hQnA2QnB2[nCentBin][nHarmonics];
  TH1D* hQnA3QnB3[nCentBin][nHarmonics];
  TH1D* hQnA4QnB4[nCentBin][nHarmonics];
  TH1D* hQnA4xQnB4x[nCentBin][nHarmonics];
  TH1D* hQnA5QnB5[nCentBin][nHarmonics];
  TH1D* hQnA6QnB6[nCentBin][nHarmonics];
  TH1D* hQnA6xQnB6x[nCentBin][nHarmonics];
  TH1D* hqnPlusQA1Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnPlusQA2Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnPlusQA3Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnPlusQA4Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnPlusQA4xStar[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnPlusQA5Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnPlusQA6Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnPlusQA6xStar[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnMinusQB1Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnMinusQB2Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnMinusQB3Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnMinusQB4Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnMinusQB4xStar[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnMinusQB5Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnMinusQB6Star[nCentBin][nHarmonics][nPtBin];
  TH1D* hqnMinusQB6xStar[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA1QnB1ptp[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA1QnB1ptm[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA2QnB2ptp[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA2QnB2ptm[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA3QnB3ptp[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA3QnB3ptm[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA4QnB4ptp[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA4QnB4ptm[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA4xQnB4xptp[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA4xQnB4xptm[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA5QnB5ptp[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA5QnB5ptm[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA6QnB6ptp[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA6QnB6ptm[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA6xQnB6xptp[nCentBin][nHarmonics][nPtBin];
  TH1D* hQnA6xQnB6xptm[nCentBin][nHarmonics][nPtBin];


  TH1D* hMeanPt[nCentBin][nPtBin];
  TH1D* hMeanPtMid[nCentBin][nPtBin];
  TH1D* hPtCent[nCentBin];
  TH1D* hPtCentCorr[nCentBin];
  TH1D *histcent;
  TH1D* hSumHF3[nCentBin];
  TH1D* hSumHF4[nCentBin];
  TH1D* hSumHF3Pt[nCentBin][nPtBin]; // pT dependent
  TH1D* hSumHF4Pt[nCentBin][nPtBin];
  TH1D* hNtrk10[nCentBin];
  TH1D* hNtrk24[nCentBin];
  TH1D* hNtrkoff[nCentBin];
  TH1D* hNtrk24Corr[nCentBin];
  //TH1D* hNtrk10Pt[nCentBin][nPtBin];
  //TH1D* hNtrk24Pt[nCentBin][nPtBin];
  TH1D* hSumTrkPt10[nCentBin];
  TH1D* hSumTrkPt24[nCentBin];
  //TH1D* hSumTrkPt10Pt[nCentBin][nPtBin];
  //TH1D* hSumTrkPt24Pt[nCentBin][nPtBin];
  TH2D *hNtrkCent;
  TH2D *hNtrkCorrCent;

  double sumwPlus[nCentBin];
  double sumwMinus[nCentBin];
  double sumwMinusPlus[nCentBin];
  double sumwDEMinusPlus[nCentBin];
  double sumwPlusTrk[nCentBin][nPtBin];
  double sumwMinusTrk[nCentBin][nPtBin];
  double sumwMidTrk[nCentBin];
  double sumwMidTrkPt[nCentBin][nPtBin];

  double vSumCosn[nCentBin][nHarmonics][nPtBin]; //tuos
  double vSumTrk10[nCentBin][nPtBin]; //tuos

  // ## track selection ##
  double etamin_;    //min eta of the tracks
  double etamax_;    //max eta of the tracks
  double ptmin_;     //min pt of the tracks
  double ptmax_;     //max pt of the tracks
  double dzdzerror_; //cut on dz/dzerror of the tracks
  double d0d0error_; //cut on d0/d0error of the tracks
  double pterrorpt_; //cut on pterror/pt of the tracks
  // ## vertex selection ##
  double  minvz_;         //minimum z distance wrt (0,0,0) for the vertex
  double  maxvz_;         //maximum z distance wrt (0,0,0) for the vertex
  // ## file acc & eff & fake ##
  bool cweight_; //use particle weight to correct from acc X eff
  edm::InputTag fname_;   //file name that contains acc X eff corrections
  TFile* feff_;                 //TFile that contains 2D histos (pt, eta) with eff/(1-fak)
  //std::vector<TH2D*> heff_;     //TH2D histograms used for correction
  TH3F* hEff_3D;
  TH3F* hFak_3D;
  TH3F* hMul_3D;
  TH3F* hSec_3D;
  bool doOO_;
  bool doNeNe_;
  bool runreco_;
  bool rungen_;

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  edm::ESGetToken<SetupData, SetupRecord> setupToken_;
#endif
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SPAnalysis::SPAnalysis(const edm::ParameterSet& iConfig) :
      tracksToken_( consumes< edm::View< pat::PackedCandidate > >( iConfig.getParameter< edm::InputTag >( "tracks" ) ) ),
      pfToken_(consumes<pat::PackedCandidateCollection>(iConfig.getParameter<edm::InputTag>("pfCands"))),
      vertexToken_( consumes< reco::VertexCollection >( iConfig.getParameter< edm:: InputTag >( "vertices" ) ) ),
      genTracksToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genTracks"))),
      etaCutHF_(iConfig.getUntrackedParameter<double>("etaCutHF")),
      etaCutTrk_(iConfig.getUntrackedParameter<double>("etaCutTrk")),
      etaCutTrkNonC_(iConfig.getUntrackedParameter<double>("etaCutTrkNonC")),

  //track selection
  etamin_(iConfig.getUntrackedParameter<double>("etamin")),
  etamax_(iConfig.getUntrackedParameter<double>("etamax")),
  ptmin_(iConfig.getUntrackedParameter<double>("ptmin")),
  ptmax_(iConfig.getUntrackedParameter<double>("ptmax")),
  dzdzerror_(iConfig.getUntrackedParameter<double>("dzdzerror")),
  d0d0error_(iConfig.getUntrackedParameter<double>("d0d0error")),
  pterrorpt_(iConfig.getUntrackedParameter<double>("pterrorpt")),
  //vertex selection
  minvz_(iConfig.getUntrackedParameter<double>("minvz")),
  maxvz_(iConfig.getUntrackedParameter<double>("maxvz")),
  //file acc & eff & fake
  cweight_(iConfig.getUntrackedParameter<bool>("cweight")),
  fname_(iConfig.getUntrackedParameter<edm::InputTag>("fname")),
  doOO_(iConfig.getUntrackedParameter<bool>("doOO")),
  doNeNe_(iConfig.getUntrackedParameter<bool>("doNeNe")),
  runreco_(iConfig.getUntrackedParameter<bool>("runreco")),
  rungen_(iConfig.getUntrackedParameter<bool>("rungen"))
{

   TString filename(fname_.label().c_str());
   feff_ = 0x0;
   if(cweight_ && !filename.IsNull()){
      edm::FileInPath fip(Form("OOAnalysisFlow/SPAnalysis/data/%s",filename.Data()));
      feff_ = new TFile(fip.fullPath().c_str(),"READ");
      //heff_.resize(feff_->GetNkeys());
      //for(unsigned short ik = 0; ik < heff_.size(); ++ik)
      //{
         //heff_[ik] = (TH2D*) feff_->Get(feff_->GetListOfKeys()->At(ik)->GetName());
         //cout<<"hist eff: "<<feff_->GetListOfKeys()->At(ik)->GetName()<<endl;
      //}
      
      hEff_3D = (TH3F*)feff_->Get("hEff_3D");
      hFak_3D = (TH3F*)feff_->Get("hFak_3D");
      hMul_3D = (TH3F*)feff_->Get("hMul_3D");
      hSec_3D = (TH3F*)feff_->Get("hSec_3D");      

      feff_->Close();

   }

  TH1::SetDefaultSumw2();

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  setupDataToken_ = esConsumes<SetupData, SetupRecord>();
#endif
  //now do what ever initialization is needed
}

SPAnalysis::~SPAnalysis() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  //
  // please remove this method altogether if it would be left empty
}

//
// member functions
//

// ------------ method called for each event  ------------
void SPAnalysis::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //re-centering parameters
    for(int ibin=0; ibin<nCentBin; ibin++){
      for(int iH=0; iH<nHarmonics; iH++){
        if(doOO_){ 
          meanQhfPlusX[ibin][iH] = OOmeanQhfPlusX[ibin][iH];
          meanQhfPlusY[ibin][iH] = OOmeanQhfPlusY[ibin][iH];
          meanQhfMinusX[ibin][iH] = OOmeanQhfMinusX[ibin][iH];
          meanQhfMinusY[ibin][iH] = OOmeanQhfMinusY[ibin][iH];
          meanQhfMinusPlusX[ibin][iH] = OOmeanQhfMinusPlusX[ibin][iH];
          meanQhfMinusPlusY[ibin][iH] = OOmeanQhfMinusPlusY[ibin][iH];
          meanQDEMinusPlusX[ibin][iH] = OOmeanQDEMinusPlusX[ibin][iH];
          meanQDEMinusPlusY[ibin][iH] = OOmeanQDEMinusPlusY[ibin][iH];
          meanQtrkMidX[ibin][iH] = OOmeanQtrkMidX[ibin][iH];
          meanQtrkMidY[ibin][iH] = OOmeanQtrkMidY[ibin][iH];
          for(int iPt=0; iPt<nPtBin; iPt++){
            meanQtrkPtPlusX[ibin][iH][iPt] = OOmeanQtrkPtPlusX[ibin][iH][iPt];
            meanQtrkPtPlusY[ibin][iH][iPt] = OOmeanQtrkPtPlusY[ibin][iH][iPt];
            meanQtrkPtMinusX[ibin][iH][iPt] = OOmeanQtrkPtMinusX[ibin][iH][iPt];
            meanQtrkPtMinusY[ibin][iH][iPt] = OOmeanQtrkPtMinusY[ibin][iH][iPt];
            meanQtrkPtMidX[ibin][iH][iPt] = OOmeanQtrkPtMidX[ibin][iH][iPt];
            meanQtrkPtMidY[ibin][iH][iPt] = OOmeanQtrkPtMidY[ibin][iH][iPt];
          }
        } //doOO
        if(doNeNe_){
          meanQhfPlusX[ibin][iH] = NeNemeanQhfPlusX[ibin][iH];
          meanQhfPlusY[ibin][iH] = NeNemeanQhfPlusY[ibin][iH];
          meanQhfMinusX[ibin][iH] = NeNemeanQhfMinusX[ibin][iH];
          meanQhfMinusY[ibin][iH] = NeNemeanQhfMinusY[ibin][iH];
          meanQhfMinusPlusX[ibin][iH] = NeNemeanQhfMinusPlusX[ibin][iH];
          meanQhfMinusPlusY[ibin][iH] = NeNemeanQhfMinusPlusY[ibin][iH];
          meanQDEMinusPlusX[ibin][iH] = NeNemeanQDEMinusPlusX[ibin][iH];
          meanQDEMinusPlusY[ibin][iH] = NeNemeanQDEMinusPlusY[ibin][iH];
          meanQtrkMidX[ibin][iH] = NeNemeanQtrkMidX[ibin][iH];
          meanQtrkMidY[ibin][iH] = NeNemeanQtrkMidY[ibin][iH];
          for(int iPt=0; iPt<nPtBin; iPt++){
            meanQtrkPtPlusX[ibin][iH][iPt] = NeNemeanQtrkPtPlusX[ibin][iH][iPt];
            meanQtrkPtPlusY[ibin][iH][iPt] = NeNemeanQtrkPtPlusY[ibin][iH][iPt];
            meanQtrkPtMinusX[ibin][iH][iPt] = NeNemeanQtrkPtMinusX[ibin][iH][iPt];
            meanQtrkPtMinusY[ibin][iH][iPt] = NeNemeanQtrkPtMinusY[ibin][iH][iPt];
            meanQtrkPtMidX[ibin][iH][iPt] = NeNemeanQtrkPtMidX[ibin][iH][iPt];
            meanQtrkPtMidY[ibin][iH][iPt] = NeNemeanQtrkPtMidY[ibin][iH][iPt];
          }
        } //doNeNe
      } //iH
    } //ibin

//hfpf for centrality
    double etHFtowerSumPlus = 0;
    double etHFtowerSumMinus = 0;
    double etHFtowerSum = 0;
    double etHFtowerSumPlus4 = 0;
    double etHFtowerSumMinus4 = 0;
    double etHFtowerSum4 = 0;
  edm::Handle<pat::PackedCandidateCollection> pfs;
  iEvent.getByToken(pfToken_, pfs);
  for(unsigned int i = 0, n = pfs->size(); i < n; ++i) {
    const pat::PackedCandidate &pf = (*pfs)[i];
    // pdfId 1 and 2 for hadronic and em particles in HF
    if(pf.pdgId() == 1 || pf.pdgId() == 2){
      //if(pf.et()<0.01 || pf.et()>30) continue;

      if(pf.eta()>3.0 && pf.eta()<5.2){
        etHFtowerSumPlus += pf.et();
        if(pf.eta()>4.0) etHFtowerSumPlus4 += pf.et();
      }
      if(pf.eta()<-3.0 && pf.eta()>-5.2){
        etHFtowerSumMinus += pf.et();
        if(pf.eta()<-4.0) etHFtowerSumMinus4 += pf.et();
      }

    } // end pf hf
  } // end pf
    etHFtowerSum = etHFtowerSumPlus + etHFtowerSumMinus;
    etHFtowerSum4 = etHFtowerSumPlus4 + etHFtowerSumMinus4;
    //cout<<"etHFtowerSum = "<<etHFtowerSum<<endl;

    // get centrality bins [0, 200)
    int cent = -1;
    if(doOO_){
      for (int ic=0; ic < nCentTableBins; ic++){
        if(etHFtowerSum>oo_binBoundaries[ic] && etHFtowerSum<=oo_binBoundaries[ic+1]){
          cent = nCentTableBins -1 - ic;
        }
      }
   } else if (doNeNe_){
      for (int ic=0; ic < nCentTableBins; ic++){
        if(etHFtowerSum>nene_binBoundaries[ic] && etHFtowerSum<=nene_binBoundaries[ic+1]){
          cent = nCentTableBins -1 - ic;
        }
      }
   }

   int ibin=-1;
   for(int j=0;j<nCentBin;j++){
     if(cent>=2*cBin[j]&&cent<2*cBin[j+1]){
     //if(cent>=2*cBinlow[j]&&cent<2*cBinhigh[j]){
       ibin=j;
     } 
   }
   //if(ibin<0 || ibin==nCentBin) return;
   //cout<<"etHFtowerSum = "<<etHFtowerSum<<"  cent = "<<cent<<"  ibin = "<<ibin<<endl;
   if(cent<0 || cent>119) return; //may need to change if the ranges are changed; 119; 159; 175
   histcent->Fill(cent);
   hSumHF3[ibin]->Fill(etHFtowerSum);
   hSumHF4[ibin]->Fill(etHFtowerSum4);

  for(int iC=0;iC<nCentBin;iC++){
     sumwPlus[iC] = 0.0;
     sumwMinus[iC] = 0.0;
     sumwMinusPlus[iC] = 0.0;
     sumwDEMinusPlus[iC] = 0.0;
     sumwMidTrk[iC] = 0.0;

//tuos-sub
     sumwA1[iC] = 0.0;
     sumwA2[iC] = 0.0;
     sumwA3[iC] = 0.0;
     sumwA4[iC] = 0.0;
     sumwA4x[iC] = 0.0;
     sumwA5[iC] = 0.0;
     sumwA6[iC] = 0.0;
     sumwA6x[iC] = 0.0;
     sumwB1[iC] = 0.0;
     sumwB2[iC] = 0.0;
     sumwB3[iC] = 0.0;
     sumwB4[iC] = 0.0;
     sumwB4x[iC] = 0.0;
     sumwB5[iC] = 0.0;
     sumwB6[iC] = 0.0;
     sumwB6x[iC] = 0.0;

     for(int iPt=0; iPt<nPtBin; iPt++){
       sumwPlusTrk[iC][iPt] = 0.0;
       sumwMinusTrk[iC][iPt] = 0.0;
       sumwMidTrkPt[iC][iPt] = 0.0;

       vSumTrk10[iC][iPt] = 0.0; //tuos
     }
     for(int iH=0; iH<nHarmonics; iH++){
       QhfPlus[iC][iH] = 0;
       QhfMinus[iC][iH] = 0;
       QhfMinusPlus[iC][iH] = 0;
       QDEMinusPlus[iC][iH] = 0;
       QtrkMid[iC][iH] = 0;
//tuos-sub
       QA1[iC][iH] = 0;
       QA2[iC][iH] = 0;
       QA3[iC][iH] = 0;
       QA4[iC][iH] = 0;
       QA4x[iC][iH] = 0;
       QA5[iC][iH] = 0;
       QA6[iC][iH] = 0;
       QA6x[iC][iH] = 0;
       QB1[iC][iH] = 0;
       QB2[iC][iH] = 0;
       QB3[iC][iH] = 0;
       QB4[iC][iH] = 0;
       QB4x[iC][iH] = 0;
       QB5[iC][iH] = 0;
       QB6[iC][iH] = 0;
       QB6x[iC][iH] = 0;

       for(int iPt=0; iPt<nPtBin; iPt++){
         QtrkPtPlus[iC][iH][iPt] = 0;
         QtrkPtMinus[iC][iH][iPt] = 0;
         QtrkPtMid[iC][iH][iPt] = 0;

         vSumCosn[iC][iH][iPt] = 0; //tuos
       }
     }
   }


   //track collection
   auto trks = iEvent.getHandle( tracksToken_ );

   //access tracks chi2/ndf
   //auto chi2Map = iEvent.getHandle( chi2Map_ );

   //vtx collection
   auto pvs = iEvent.getHandle( vertexToken_ );

   //best vertex
   double bestvzError;
   math::XYZPoint bestvtx;
   math::Error<3>::type vtx_cov;
   if ( !pvs->empty() ) {
      const reco::Vertex& vtx = (*pvs)[0];
      bestvzError = vtx.zError();
      bestvtx = vtx.position();
      vtx_cov = vtx.covariance();
   }else {
      return;
   }

   //cout<<"nVertex = "<<pvs->size()<<endl;
   if (pvs->size() > 1) return;

   if ( std::abs( bestvtx.z() ) > 15. ) return;
   //if (  std::sqrt(bestvtx.x()*bestvtx.x() + bestvtx.y()*bestvtx.y())  > 0.2 ) return;
   if ( std::abs( bestvtx.z() ) < minvz_ || std::abs( bestvtx.z() ) > maxvz_ ) return;
   //cout<<"vz = "<<bestvtx.z()<<endl;


int indexHF3Pt = 0; //tuos
int indexHF4Pt = 0; //tuos
int indexNtrk10 = 0;
int indexNtrkoff = 0;
int indexNtrk24 = 0;
double indexNtrk24Corr = 0;
//int indexNtrk10Pt = 0;
//int indexNtrk24Pt = 0;
double sumTrkPt10 = 0;
double sumTrkPt24 = 0;
if(runreco_){
   for (auto const& trk : *trks) {
      if ( !trk.hasTrackDetails() ) continue;
      auto iter_tk = trk.pseudoTrack();

     gRandom->SetSeed(0);
     double randomAngle_tk = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
     //cout << randomAngle_tk << endl;

     //cout << iter_tk.phi() << endl;
     if ( iter_tk.charge() == 0 ) continue;
     if ( !iter_tk.quality(reco::TrackBase::highPurity) ) continue;
     if ( iter_tk.eta() > etamax_ ) continue;
     if ( iter_tk.eta() < etamin_ ) continue;
     if ( iter_tk.pt() > ptmax_ ) continue;
     if ( iter_tk.pt() < ptmin_ ) continue;
      double dzvtx = iter_tk.dz( bestvtx );
      double dxyvtx = iter_tk.dxy( bestvtx );
      double dzerror = std::hypot( iter_tk.dzError(), bestvzError );
      double dxyerror = iter_tk.dxyError( bestvtx, vtx_cov );
      if ( iter_tk.pt() > 10 &&  std::abs( iter_tk.ptError() ) / iter_tk.pt() >= pterrorpt_ ) continue; //default 0.1 at pt>10
      if ( std::abs( dzvtx / dzerror ) >= dzdzerror_ ) continue; //default 3/loose 5
      if ( std::abs( dxyvtx / dxyerror ) >= d0d0error_ ) continue; //default 3/loose 5

      int cpt=-1;
      for(int ipt=0; ipt<nPtBin; ipt++){
        if(iter_tk.pt()>ptBins[ipt] && iter_tk.pt()<=ptBins[ipt+1])
          cpt=ipt;
      }
      if(cpt<0) cpt=nPtBin-1;

      //if(iter_tk.pt()>9)  cout<<" cpt = "<<cpt<<" pt = "<<iter_tk.pt()<<"  pterr/pt = "<<std::abs( iter_tk.ptError() ) / iter_tk.pt()<<endl;
   
     double weight = 1.0;
     if( cweight_ ){
       /*int idx = 0;
       if( feff_ && heff_[idx] ) {
         double eff = heff_[idx]->GetBinContent(heff_[idx]->FindBin(iter_tk.eta(),iter_tk.pt()));
         if( eff != 0. ){
           weight = 1./eff;
         }
       }*/
       double eff_pt_tmp = iter_tk.pt();
       if(eff_pt_tmp>6.1) eff_pt_tmp=6.1; // fix eff fluctuation vs pt // 6.1; 9.9
       //double eff = hEff_3D->GetBinContent(hEff_3D->GetXaxis()->FindBin(iter_tk.eta()), hEff_3D->GetYaxis()->FindBin(iter_tk.pt()), hEff_3D->GetZaxis()->FindBin(cent));
       //double fake = hFak_3D->GetBinContent(hFak_3D->GetXaxis()->FindBin(iter_tk.eta()), hFak_3D->GetYaxis()->FindBin(iter_tk.pt()), hFak_3D->GetZaxis()->FindBin(cent)) ;
       //double mult = hMul_3D->GetBinContent(hMul_3D->GetXaxis()->FindBin(iter_tk.eta()), hMul_3D->GetYaxis()->FindBin(iter_tk.pt()), hMul_3D->GetZaxis()->FindBin(cent));
       //double sec = hSec_3D->GetBinContent(hSec_3D->GetXaxis()->FindBin(iter_tk.eta()), hSec_3D->GetYaxis()->FindBin(iter_tk.pt()), hSec_3D->GetZaxis()->FindBin(cent));    
       double eff = hEff_3D->GetBinContent(hEff_3D->GetXaxis()->FindBin(iter_tk.eta()), hEff_3D->GetYaxis()->FindBin(eff_pt_tmp), hEff_3D->GetZaxis()->FindBin(cent));
       double fake = hFak_3D->GetBinContent(hFak_3D->GetXaxis()->FindBin(iter_tk.eta()), hFak_3D->GetYaxis()->FindBin(eff_pt_tmp), hFak_3D->GetZaxis()->FindBin(cent));
       double mult = hMul_3D->GetBinContent(hMul_3D->GetXaxis()->FindBin(iter_tk.eta()), hMul_3D->GetYaxis()->FindBin(eff_pt_tmp), hMul_3D->GetZaxis()->FindBin(cent));
       double sec = hSec_3D->GetBinContent(hSec_3D->GetXaxis()->FindBin(iter_tk.eta()), hSec_3D->GetYaxis()->FindBin(eff_pt_tmp), hSec_3D->GetZaxis()->FindBin(cent));
     if(eff > 0.001){
       weight = (1 - fake)*(1 - sec)/eff/(1 + mult) ;
     }else{
       weight = 1.0 ;
     }
     if((iter_tk.pt() < 0 || iter_tk.pt() > 500) || (fabs(iter_tk.eta()) > 2.4)){
      weight = 0 ;
     }   
     //cout<<"cent ="<<cent<<" , pt ="<<iter_tk.pt()<<"  eff ="<<eff<<"  fake ="<<fake<<" , weight ="<<weight<<endl;
    }

      if(iter_tk.pt()>0.3){
        indexNtrk24 += 1;
        indexNtrk24Corr += weight;
        sumTrkPt24 += iter_tk.pt();
        if(iter_tk.pt()>0.4){
          indexNtrkoff += 1;
        }
        if(iter_tk.eta()>-1 && iter_tk.eta()<1){
          indexNtrk10 += 1;
          sumTrkPt10 += iter_tk.pt();
        }
      }
      //hMeanPt[ibin][cpt]->Fill(iter_tk.pt(), weight);

      if(indexHF3Pt==0){
        hSumHF3Pt[ibin][cpt]->Fill(etHFtowerSum);//tuos
        indexHF3Pt += 1;
      }
      if(indexHF4Pt==0){
        hSumHF4Pt[ibin][cpt]->Fill(etHFtowerSum4);//tuos
        indexHF4Pt += 1;
      }
      //if(indexNtrk10Pt==0){
      //  hNtrk10Pt[ibin][cpt]->Fill();//tuos
      //  indexNtrk10Pt += 1;
      //}
/*
//tuos : Sum
if(iter_tk.eta()>-1 && iter_tk.eta()<1){
  vSumTrk10[ibin][cpt] += 1;

   for (auto const& trk2 : *trks) {
      if ( !trk2.hasTrackDetails() ) continue;
      auto iter_tk2 = trk2.pseudoTrack();

     if ( iter_tk2.charge() == 0 ) continue;
     if ( !iter_tk2.quality(reco::TrackBase::highPurity) ) continue;
     if ( iter_tk2.eta() > etamax_ ) continue;
     if ( iter_tk2.eta() < etamin_ ) continue;
     if ( iter_tk2.pt() > ptmax_ ) continue;
     if ( iter_tk2.pt() < ptmin_ ) continue;
      double dzvtx = iter_tk2.dz( bestvtx );
      double dxyvtx = iter_tk2.dxy( bestvtx );
      double dzerror = std::hypot( iter_tk2.dzError(), bestvzError );
      double dxyerror = iter_tk2.dxyError( bestvtx, vtx_cov );
      if ( iter_tk2.pt() > 10 && std::abs( iter_tk2.ptError() ) / iter_tk2.pt() >= pterrorpt_ ) continue; //default 0.1 at pt>10
      if ( std::abs( dzvtx / dzerror ) >= dzdzerror_ ) continue; //default 3/loose 5
      if ( std::abs( dxyvtx / dxyerror ) >= d0d0error_ ) continue; //default 3/loose 5

     if(iter_tk2.eta()>1 || iter_tk2.eta()<-1) continue;
     if(iter_tk2.pt()>0.3){
       for(int iH=0; iH<nHarmonics; iH++){
         vSumCosn[ibin][iH][cpt] += cos( (iH+1)*(iter_tk.phi() - iter_tk2.phi()) ); //tuos
       }
     }

  }
}
//tuos
*/

      double pt_weight = 1.0; 
      if(iter_tk.eta() > 1.0*etaCutTrkNonC_ && iter_tk.eta() < 1.0*etaCutTrk_){
        sumwPlusTrk[ibin][cpt]+= weight; // pt_weight->weight
        for(int iH=0; iH<nHarmonics; iH++){
          QtrkPtPlus[ibin][iH][cpt]+= weight * TComplex::Exp(TComplex(0,(iH+1)*randomAngle_tk)); //iter_tk.phi()));
        }
      }
      if(iter_tk.eta() < -1.0*etaCutTrkNonC_ && iter_tk.eta() > -1.0*etaCutTrk_){
        sumwMinusTrk[ibin][cpt]+= weight; // pt_weight->weight
        for(int iH=0; iH<nHarmonics; iH++){
          QtrkPtMinus[ibin][iH][cpt]+= weight * TComplex::Exp(TComplex(0,(iH+1)*randomAngle_tk)); //iter_tk.phi()));
        }
      }
      if(iter_tk.eta() > -0.75 && iter_tk.eta() < 0.75 && iter_tk.pt() > 0.3 && iter_tk.pt() < 3.0){ //Q_nC 3sub
        sumwMidTrk[ibin] += iter_tk.pt(); // keep pt weight for now
        for(int iH=0; iH<nHarmonics; iH++){
          QtrkMid[ibin][iH]+= iter_tk.pt() * TComplex::Exp(TComplex(0,(iH+1)*randomAngle_tk)); //iter_tk.phi()));
        }
      }
//tuos-sub
      if(iter_tk.eta() > -2.4 && iter_tk.eta() < -1.0){
        sumwA1[ibin] += weight; // pt_weight->weight
        for(int iH=0; iH<nHarmonics; iH++){
          QA1[ibin][iH]+= weight * TComplex::Exp(TComplex(0,(iH+1)*randomAngle_tk)); //iter_tk.phi()));
        }
      }
      if(iter_tk.eta() > 1 && iter_tk.eta() < 2.4){
        sumwB1[ibin] += weight; // pt_weight->weight
        for(int iH=0; iH<nHarmonics; iH++){
          QB1[ibin][iH]+= weight * TComplex::Exp(TComplex(0,(iH+1)*randomAngle_tk)); //iter_tk.phi()));
        }
      }
      //pT, Trk, Qvector
      //if(iter_tk.eta() > -1.0*etaCutTrk_ && iter_tk.eta() < 1.0*etaCutTrk_ ){
      if( (iter_tk.eta() < -1.0*etaCutTrkNonC_ && iter_tk.eta() > -1.0*etaCutTrk_)   ||   (iter_tk.eta() > 1.0*etaCutTrkNonC_ && iter_tk.eta() < 1.0*etaCutTrk_) ){
        hMeanPt[ibin][cpt]->Fill(iter_tk.pt(), weight);
        hPtCent[ibin]->Fill(iter_tk.pt());
        hPtCentCorr[ibin]->Fill(iter_tk.pt(), weight);
        hMeanPtMid[ibin][cpt]->Fill(iter_tk.pt(), weight);
        sumwMidTrkPt[ibin][cpt]+= pt_weight;
        for(int iH=0; iH<nHarmonics; iH++){
          QtrkPtMid[ibin][iH][cpt]+= pt_weight * TComplex::Exp(TComplex(0,(iH+1)*randomAngle_tk)); //iter_tk.phi()));
        }
      }

  } // end trk
}

hNtrk10[ibin]->Fill(indexNtrk10); //tuos
hNtrk24[ibin]->Fill(indexNtrk24);
hNtrkoff[ibin]->Fill(indexNtrkoff);
hNtrk24Corr[ibin]->Fill(indexNtrk24Corr);
hNtrkCent->Fill(cent, indexNtrk24);
hNtrkCorrCent->Fill(cent, indexNtrk24Corr);
hSumTrkPt10[ibin]->Fill(sumTrkPt10);
hSumTrkPt24[ibin]->Fill(sumTrkPt24);
//for(int iH=0; iH<nHarmonics; iH++){
//  for(int iPt=0; iPt<nPtBin; iPt++){
//    if(vSumTrk10[ibin][iPt]>0){
//      hSumCosn[ibin][iH][iPt]->Fill( vSumCosn[ibin][iH][iPt] - 1 );
//    }
//  }
//}
//cout<<"etHFtowerSum = "<<etHFtowerSum<<"  etHFtowerSum4 = "<<etHFtowerSum4<<"  indexHF3Pt: "<<indexHF3Pt<<"  indexHF4Pt:"<<indexHF4Pt<<"  indexNtrk10:"<<indexNtrk10<<"  indexNtrk24:"<<indexNtrk24<<"  sumTrkPt10="<<sumTrkPt10<<"  sumTrkPt24="<<sumTrkPt24<<endl;

/*
if(rungen_){
   edm::Handle< reco::GenParticleCollection > tracks;
   iEvent.getByToken(genTracksToken_, tracks);
   for( reco::GenParticleCollection::const_iterator itTrk = tracks->begin() ; itTrk != tracks->end() ; ++itTrk ) {
     double eta    = itTrk->eta();
     double pt     = itTrk->pt();
     double phi     = itTrk->phi();

      if(eta>2.4 || eta<-2.4) continue;

      int cpt=-1;
      for(int ipt=0; ipt<nPtBin; ipt++){
        if(pt>ptBins[ipt] && pt<ptBins[ipt+1])
          cpt=ipt;
      }
      if(cpt<0) cpt=nPtBin-1;
      //cout<<" cpt = "<<cpt<<" pt = "<<track.pt()<<endl;

      hMeanPt[ibin][cpt]->Fill(pt);

      double pt_weight = 1.0;
      if(eta > 0 && eta < 1.0*etaCutTrk_){
        sumwPlusTrk[ibin][cpt]+= pt_weight;
        for(int iH=0; iH<nHarmonics; iH++){
          QtrkPtPlus[ibin][iH][cpt]+= pt_weight * TComplex::Exp(TComplex(0,(iH+1)*phi));
        }
      }
      if(eta < 0 && eta > -1.0*etaCutTrk_){
        sumwMinusTrk[ibin][cpt]+= pt_weight;
        for(int iH=0; iH<nHarmonics; iH++){
          QtrkPtMinus[ibin][iH][cpt]+= pt_weight * TComplex::Exp(TComplex(0,(iH+1)*phi));
        }
      }
      if(eta > -0.75 && eta < 0.75 && pt > 0.3 && pt < 3.0){
        sumwMidTrk[ibin] += pt;
        for(int iH=0; iH<nHarmonics; iH++){
          QtrkMid[ibin][iH]+= pt * TComplex::Exp(TComplex(0,(iH+1)*phi));
        }
      }

   } // gen trk
}
*/

    for(int iH=0; iH<nHarmonics; iH++){
      for(int iPt=0; iPt<nPtBin; iPt++){
        if(sumwPlusTrk[ibin][iPt]>0) { QtrkPtPlus[ibin][iH][iPt]=QtrkPtPlus[ibin][iH][iPt]/sumwPlusTrk[ibin][iPt] - TComplex(meanQtrkPtPlusX[ibin][iH][iPt], meanQtrkPtPlusY[ibin][iH][iPt]); }
        else {QtrkPtPlus[ibin][iH][iPt]=0; }
        if(sumwMinusTrk[ibin][iPt]>0) { QtrkPtMinus[ibin][iH][iPt]=QtrkPtMinus[ibin][iH][iPt]/sumwMinusTrk[ibin][iPt] - TComplex(meanQtrkPtMinusX[ibin][iH][iPt], meanQtrkPtMinusY[ibin][iH][iPt]); }
        else { QtrkPtMinus[ibin][iH][iPt]=0; }
        if(sumwMidTrkPt[ibin][iPt]>0) {
          QtrkPtMid[ibin][iH][iPt] = QtrkPtMid[ibin][iH][iPt]/sumwMidTrkPt[ibin][iPt] - TComplex(meanQtrkPtMidX[ibin][iH][iPt], meanQtrkPtMidY[ibin][iH][iPt]); 
        } else { QtrkPtMid[ibin][iH][iPt] = 0; }
      }
      if(sumwMidTrk[ibin]>0){
        QtrkMid[ibin][iH] = QtrkMid[ibin][iH]/sumwMidTrk[ibin] - TComplex(meanQtrkMidX[ibin][iH], meanQtrkMidY[ibin][iH]);
      } else { QtrkMid[ibin][iH]=0; }
    }//iH
 //cout<<"end of trk"<<endl;

    // hfpf
  for(unsigned int i = 0, n = pfs->size(); i < n; ++i) {
    const pat::PackedCandidate &pf = (*pfs)[i];

    gRandom->SetSeed(0);
    double randomAngleHF = gRandom->Uniform(-TMath::Pi(), TMath::Pi());
    //cout << randomAngleHF << endl;

    // pdfId 1 and 2 for hadronic and em particles in HF
    if(pf.pdgId() == 1 || pf.pdgId() == 2){
       double eta = pf.eta();
       //cout<<"eta = "<<eta<<endl;
          if(eta>3.0 && eta<5.2 && eta > 1.0*etaCutHF_){
            sumwPlus[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QhfPlus[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if(eta<-3.0 && eta>-5.2 && eta < -1.0*etaCutHF_){
            sumwMinus[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QhfMinus[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>-5.2 && eta < -1.0*etaCutHF_) || (eta > 1.0*etaCutHF_ && eta<5.2) ){
            sumwMinusPlus[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QhfMinusPlus[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
//tuos-sub
          if( (eta>-5.2 && eta < -3.0) ){
            sumwA4[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QA4[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>3 && eta < 5.2) ){
            sumwB4[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QB4[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }

          if( (eta>-5.2 && eta < -3.5) ){//tuosx
            sumwA4x[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QA4x[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>3.5 && eta < 5.2) ){
            sumwB4x[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QB4x[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }

          if( (eta>-5.2 && eta < -4.0) ){
            sumwA5[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QA5[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>4 && eta < 5.2) ){
            sumwB5[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QB5[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>-5.2 && eta < -4.5) ){
            sumwA6[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QA6[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>4.5 && eta < 5.2) ){
            sumwB6[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QB6[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }

          if( (eta>-5.2 && eta < -4.8) ){
            sumwA6x[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QA6x[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>4.8 && eta < 5.2) ){
            sumwB6x[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QB6x[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }

      } // end hf
          double eta = pf.eta();
          if( (eta>-3.0 && eta < -1.0) || (eta > 1.0 && eta<3.0) ){
            sumwDEMinusPlus[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QDEMinusPlus[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
//tuos-sub
          if( (eta>-2.4 && eta < -1.0) ){
            sumwA2[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QA2[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>1 && eta < 2.4) ){
            sumwB2[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QB2[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>-3 && eta < -2.0) ){
            sumwA3[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QA3[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
          if( (eta>2 && eta < 3) ){
            sumwB3[ibin]+= pf.et();
            for(int iH=0; iH<nHarmonics; iH++){
              QB3[ibin][iH]+= pf.et() * TComplex::Exp(TComplex(0,(iH+1)*randomAngleHF)); //pf.phi()));
            }
          }
      //} // end hf
    } // end hfpf
    for(int iH=0; iH<nHarmonics; iH++){
      if(sumwPlus[ibin]>0) { QhfPlus[ibin][iH]=QhfPlus[ibin][iH]/sumwPlus[ibin] - TComplex(meanQhfPlusX[ibin][iH], meanQhfPlusY[ibin][iH]); }
      else { QhfPlus[ibin][iH]=0; }
      if(sumwMinus[ibin]>0) { QhfMinus[ibin][iH]=QhfMinus[ibin][iH]/sumwMinus[ibin] - TComplex(meanQhfMinusX[ibin][iH], meanQhfMinusY[ibin][iH]); }
      else { QhfMinus[ibin][iH]=0; }
      if(sumwMinusPlus[ibin]>0) { QhfMinusPlus[ibin][iH]=QhfMinusPlus[ibin][iH]/sumwMinusPlus[ibin] - TComplex(meanQhfMinusPlusX[ibin][iH], meanQhfMinusPlusY[ibin][iH]); }
      else { QhfMinusPlus[ibin][iH]=0; }
      if(sumwDEMinusPlus[ibin]>0) { QDEMinusPlus[ibin][iH]=QDEMinusPlus[ibin][iH]/sumwDEMinusPlus[ibin] - TComplex(meanQDEMinusPlusX[ibin][iH], meanQDEMinusPlusY[ibin][iH]); }
      else { QDEMinusPlus[ibin][iH]=0; }
      //cout<<"sumwDEMinusPlus = "<<sumwPlus[ibin]<<"  "<<QhfPlus[ibin][iH]<<"  "<<sumwMinusPlus[ibin]<<"  "<<sumwDEMinusPlus[ibin]<<"  "<<QDEMinusPlus[ibin][iH]<<endl;
    }//iH

//tuos-sub, no normalization for QAx/QBx needed; 
    //cout<<"sumwA1="<<sumwA1[ibin]<<"  sumwB1="<<sumwB1[ibin]<<endl;
    //cout<<"sumwA3="<<sumwA3[ibin]<<"  sumwB3="<<sumwB3[ibin]<<endl;
    //cout<<"sumwA5="<<sumwA5[ibin]<<"  sumwB5="<<sumwB5[ibin]<<endl;
    for(int iH=0; iH<nHarmonics; iH++){
      if(sumwA1[ibin]>0 && sumwB1[ibin]>0) {
        hQnA1QnB1[ibin][iH]->Fill((QA1[ibin][iH]*TComplex::Conjugate(QB1[ibin][iH])).Re());
        //if(iH==1) cout<<"hQnA1QnB1="<<fabs((QA1[ibin][iH]*TComplex::Conjugate(QB1[ibin][iH])).Re())<<endl;
      }
      if(sumwA2[ibin]>0 && sumwB2[ibin]>0) {
        hQnA2QnB2[ibin][iH]->Fill((QA2[ibin][iH]*TComplex::Conjugate(QB2[ibin][iH])).Re());
        //if(iH==1) cout<<"hQnA2QnB2="<<fabs((QA2[ibin][iH]*TComplex::Conjugate(QB2[ibin][iH])).Re())<<endl;
      }
      if(sumwA3[ibin]>0 && sumwB3[ibin]>0) {
        hQnA3QnB3[ibin][iH]->Fill((QA3[ibin][iH]*TComplex::Conjugate(QB3[ibin][iH])).Re());
        //if(iH==1) cout<<"hQnA3QnB3="<<fabs((QA3[ibin][iH]*TComplex::Conjugate(QB3[ibin][iH])).Re())<<endl;
      }
      if(sumwA4[ibin]>0 && sumwB4[ibin]>0) {
        hQnA4QnB4[ibin][iH]->Fill((QA4[ibin][iH]*TComplex::Conjugate(QB4[ibin][iH])).Re());
        //if(iH==1) cout<<"hQnA4QnB4="<<fabs((QA4[ibin][iH]*TComplex::Conjugate(QB4[ibin][iH])).Re())<<endl;
      }
      if(sumwA4x[ibin]>0 && sumwB4x[ibin]>0) {
        hQnA4xQnB4x[ibin][iH]->Fill((QA4x[ibin][iH]*TComplex::Conjugate(QB4x[ibin][iH])).Re());
      }
      if(sumwA5[ibin]>0 && sumwB5[ibin]>0) {
        hQnA5QnB5[ibin][iH]->Fill((QA5[ibin][iH]*TComplex::Conjugate(QB5[ibin][iH])).Re());
        //if(iH==1) cout<<"hQnA5QnB5="<<fabs((QA5[ibin][iH]*TComplex::Conjugate(QB5[ibin][iH])).Re())<<endl;
      }
      if(sumwA6[ibin]>0 && sumwB6[ibin]>0) {
        hQnA6QnB6[ibin][iH]->Fill((QA6[ibin][iH]*TComplex::Conjugate(QB6[ibin][iH])).Re());
      }
      if(sumwA6x[ibin]>0 && sumwB6x[ibin]>0) {
        hQnA6xQnB6x[ibin][iH]->Fill((QA6x[ibin][iH]*TComplex::Conjugate(QB6x[ibin][iH])).Re());
      }

    //cout<<"sumwA6x="<<sumwA6x[ibin]<<"  sumwB6x="<<sumwB6x[ibin]<<endl;

      for(int iPt=0; iPt<nPtBin; iPt++){ //do we need to fill 0 if there is no particles in a pT bin? A: No
        if(sumwPlusTrk[ibin][iPt]>0 && sumwA1[ibin]>0){
          hqnPlusQA1Star[ibin][iH][iPt]->Fill((QtrkPtPlus[ibin][iH][iPt]*TComplex::Conjugate(QA1[ibin][iH])).Re());
        }
        if(sumwPlusTrk[ibin][iPt]>0 && sumwA2[ibin]>0){
          hqnPlusQA2Star[ibin][iH][iPt]->Fill((QtrkPtPlus[ibin][iH][iPt]*TComplex::Conjugate(QA2[ibin][iH])).Re());
        }
        if(sumwPlusTrk[ibin][iPt]>0 && sumwA3[ibin]>0){
          hqnPlusQA3Star[ibin][iH][iPt]->Fill((QtrkPtPlus[ibin][iH][iPt]*TComplex::Conjugate(QA3[ibin][iH])).Re());
        }
        if(sumwPlusTrk[ibin][iPt]>0 && sumwA4[ibin]>0){
          hqnPlusQA4Star[ibin][iH][iPt]->Fill((QtrkPtPlus[ibin][iH][iPt]*TComplex::Conjugate(QA4[ibin][iH])).Re());
        }
        if(sumwPlusTrk[ibin][iPt]>0 && sumwA4x[ibin]>0){
          hqnPlusQA4xStar[ibin][iH][iPt]->Fill((QtrkPtPlus[ibin][iH][iPt]*TComplex::Conjugate(QA4x[ibin][iH])).Re());
        }
        if(sumwPlusTrk[ibin][iPt]>0 && sumwA5[ibin]>0){
          hqnPlusQA5Star[ibin][iH][iPt]->Fill((QtrkPtPlus[ibin][iH][iPt]*TComplex::Conjugate(QA5[ibin][iH])).Re());
        }
        if(sumwPlusTrk[ibin][iPt]>0 && sumwA6[ibin]>0){
          hqnPlusQA6Star[ibin][iH][iPt]->Fill((QtrkPtPlus[ibin][iH][iPt]*TComplex::Conjugate(QA6[ibin][iH])).Re());
        }
        if(sumwPlusTrk[ibin][iPt]>0 && sumwA6x[ibin]>0){
          hqnPlusQA6xStar[ibin][iH][iPt]->Fill((QtrkPtPlus[ibin][iH][iPt]*TComplex::Conjugate(QA6x[ibin][iH])).Re());
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwB1[ibin]>0){
          hqnMinusQB1Star[ibin][iH][iPt]->Fill((QtrkPtMinus[ibin][iH][iPt]*TComplex::Conjugate(QB1[ibin][iH])).Re());
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwB2[ibin]>0){
          hqnMinusQB2Star[ibin][iH][iPt]->Fill((QtrkPtMinus[ibin][iH][iPt]*TComplex::Conjugate(QB2[ibin][iH])).Re());
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwB3[ibin]>0){
          hqnMinusQB3Star[ibin][iH][iPt]->Fill((QtrkPtMinus[ibin][iH][iPt]*TComplex::Conjugate(QB3[ibin][iH])).Re());
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwB4[ibin]>0){
          hqnMinusQB4Star[ibin][iH][iPt]->Fill((QtrkPtMinus[ibin][iH][iPt]*TComplex::Conjugate(QB4[ibin][iH])).Re());
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwB4x[ibin]>0){
          hqnMinusQB4xStar[ibin][iH][iPt]->Fill((QtrkPtMinus[ibin][iH][iPt]*TComplex::Conjugate(QB4x[ibin][iH])).Re());
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwB5[ibin]>0){
          hqnMinusQB5Star[ibin][iH][iPt]->Fill((QtrkPtMinus[ibin][iH][iPt]*TComplex::Conjugate(QB5[ibin][iH])).Re());
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwB6[ibin]>0){
          hqnMinusQB6Star[ibin][iH][iPt]->Fill((QtrkPtMinus[ibin][iH][iPt]*TComplex::Conjugate(QB6[ibin][iH])).Re());
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwB6x[ibin]>0){
          hqnMinusQB6xStar[ibin][iH][iPt]->Fill((QtrkPtMinus[ibin][iH][iPt]*TComplex::Conjugate(QB6x[ibin][iH])).Re());
        }
        //tuospt
        if(sumwA1[ibin]>0 && sumwB1[ibin]>0) {
          if(sumwPlusTrk[ibin][iPt]>0) hQnA1QnB1ptp[ibin][iH][iPt]->Fill((QA1[ibin][iH]*TComplex::Conjugate(QB1[ibin][iH])).Re());
          if(sumwMinusTrk[ibin][iPt]>0) hQnA1QnB1ptm[ibin][iH][iPt]->Fill((QA1[ibin][iH]*TComplex::Conjugate(QB1[ibin][iH])).Re());
        }
        if(sumwA2[ibin]>0 && sumwB2[ibin]>0) {
          if(sumwPlusTrk[ibin][iPt]>0) hQnA2QnB2ptp[ibin][iH][iPt]->Fill((QA2[ibin][iH]*TComplex::Conjugate(QB2[ibin][iH])).Re());
          if(sumwMinusTrk[ibin][iPt]>0) hQnA2QnB2ptm[ibin][iH][iPt]->Fill((QA2[ibin][iH]*TComplex::Conjugate(QB2[ibin][iH])).Re());
        }
        if(sumwA3[ibin]>0 && sumwB3[ibin]>0) {
          if(sumwPlusTrk[ibin][iPt]>0) hQnA3QnB3ptp[ibin][iH][iPt]->Fill((QA3[ibin][iH]*TComplex::Conjugate(QB3[ibin][iH])).Re());
          if(sumwMinusTrk[ibin][iPt]>0) hQnA3QnB3ptm[ibin][iH][iPt]->Fill((QA3[ibin][iH]*TComplex::Conjugate(QB3[ibin][iH])).Re());
        }
        if(sumwA4[ibin]>0 && sumwB4[ibin]>0) {
          if(sumwPlusTrk[ibin][iPt]>0) hQnA4QnB4ptp[ibin][iH][iPt]->Fill((QA4[ibin][iH]*TComplex::Conjugate(QB4[ibin][iH])).Re());
          if(sumwMinusTrk[ibin][iPt]>0) hQnA4QnB4ptm[ibin][iH][iPt]->Fill((QA4[ibin][iH]*TComplex::Conjugate(QB4[ibin][iH])).Re());
        }
        if(sumwA4x[ibin]>0 && sumwB4x[ibin]>0) {
          if(sumwPlusTrk[ibin][iPt]>0) hQnA4xQnB4xptp[ibin][iH][iPt]->Fill((QA4x[ibin][iH]*TComplex::Conjugate(QB4x[ibin][iH])).Re());
          if(sumwMinusTrk[ibin][iPt]>0) hQnA4xQnB4xptm[ibin][iH][iPt]->Fill((QA4x[ibin][iH]*TComplex::Conjugate(QB4x[ibin][iH])).Re());
        }
        if(sumwA5[ibin]>0 && sumwB5[ibin]>0) {
          if(sumwPlusTrk[ibin][iPt]>0) hQnA5QnB5ptp[ibin][iH][iPt]->Fill((QA5[ibin][iH]*TComplex::Conjugate(QB5[ibin][iH])).Re());
          if(sumwMinusTrk[ibin][iPt]>0) hQnA5QnB5ptm[ibin][iH][iPt]->Fill((QA5[ibin][iH]*TComplex::Conjugate(QB5[ibin][iH])).Re());
        }
        if(sumwA6[ibin]>0 && sumwB6[ibin]>0) {
          if(sumwPlusTrk[ibin][iPt]>0) hQnA6QnB6ptp[ibin][iH][iPt]->Fill((QA6[ibin][iH]*TComplex::Conjugate(QB6[ibin][iH])).Re());
          if(sumwMinusTrk[ibin][iPt]>0) hQnA6QnB6ptm[ibin][iH][iPt]->Fill((QA6[ibin][iH]*TComplex::Conjugate(QB6[ibin][iH])).Re());
        }
        if(sumwA6x[ibin]>0 && sumwB6x[ibin]>0) {
          if(sumwPlusTrk[ibin][iPt]>0) hQnA6xQnB6xptp[ibin][iH][iPt]->Fill((QA6x[ibin][iH]*TComplex::Conjugate(QB6x[ibin][iH])).Re());
          if(sumwMinusTrk[ibin][iPt]>0) hQnA6xQnB6xptm[ibin][iH][iPt]->Fill((QA6x[ibin][iH]*TComplex::Conjugate(QB6x[ibin][iH])).Re());
        }

      }//iPt
    }//iH


    // final histograms for flow
    for(int iH=0; iH<nHarmonics; iH++){
      if(sumwPlus[ibin]>0 && sumwMinus[ibin]>0) {
        hQnPlusQnMinusStar[ibin][iH]->Fill((QhfPlus[ibin][iH]*TComplex::Conjugate(QhfMinus[ibin][iH])).Re());
        hQnMinusQnPlusStar[ibin][iH]->Fill((QhfMinus[ibin][iH]*TComplex::Conjugate(QhfPlus[ibin][iH])).Re());
      }
      if(sumwMinus[ibin]>0 && sumwMidTrk[ibin]>0) {
        hQnMinusQnMidStar[ibin][iH]->Fill((QhfMinus[ibin][iH]*TComplex::Conjugate(QtrkMid[ibin][iH])).Re());
      }
      if(sumwPlus[ibin]>0 && sumwMidTrk[ibin]>0) {
        hQnPlusQnMidStar[ibin][iH]->Fill((QhfPlus[ibin][iH]*TComplex::Conjugate(QtrkMid[ibin][iH])).Re());
      }
      if(sumwMinusPlus[ibin]>0 && sumwDEMinusPlus[ibin]>0) {
        hQnMinusPlusQnDEStar[ibin][iH]->Fill((QhfMinusPlus[ibin][iH]*TComplex::Conjugate(QDEMinusPlus[ibin][iH])).Re());
      }
      if(sumwMinusPlus[ibin]>0 && sumwMidTrk[ibin]>0) {
        hQnMinusPlusQnMidStar[ibin][iH]->Fill((QhfMinusPlus[ibin][iH]*TComplex::Conjugate(QtrkMid[ibin][iH])).Re());
      }
      if(sumwDEMinusPlus[ibin]>0 && sumwMidTrk[ibin]>0) {
        hQnDEQnMidStar[ibin][iH]->Fill((QDEMinusPlus[ibin][iH]*TComplex::Conjugate(QtrkMid[ibin][iH])).Re());
      }

      for(int iPt=0; iPt<nPtBin; iPt++){ //do we need to fill 0 if there is no particles in a pT bin? A: No
        if(sumwPlusTrk[ibin][iPt]>0 && sumwMinus[ibin]>0){
          hqnPlusQnMinusStar[ibin][iH][iPt]->Fill((QtrkPtPlus[ibin][iH][iPt]*TComplex::Conjugate(QhfMinus[ibin][iH])).Re());
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwPlus[ibin]>0){
          hqnMinusQnPlusStar[ibin][iH][iPt]->Fill((QtrkPtMinus[ibin][iH][iPt]*TComplex::Conjugate(QhfPlus[ibin][iH])).Re());
        }
        if(sumwMidTrkPt[ibin][iPt]>0 && sumwPlus[ibin]>0) {
          hqnQnPlusStar[ibin][iH][iPt]->Fill((QtrkPtMid[ibin][iH][iPt]*TComplex::Conjugate(QhfPlus[ibin][iH])).Re());
        }
        if(sumwMidTrkPt[ibin][iPt]>0 && sumwMinus[ibin]>0) {
          hqnQnMinusStar[ibin][iH][iPt]->Fill((QtrkPtMid[ibin][iH][iPt]*TComplex::Conjugate(QhfMinus[ibin][iH])).Re());
        }
        if(sumwMidTrkPt[ibin][iPt]>0 && sumwMinusPlus[ibin]>0) {
          hqnQnMinusPlusStar[ibin][iH][iPt]->Fill((QtrkPtMid[ibin][iH][iPt]*TComplex::Conjugate(QhfMinusPlus[ibin][iH])).Re());
        }

        if(sumwPlusTrk[ibin][iPt]>0 && sumwMinus[ibin]>0 && sumwPlus[ibin]>0){ //tuos
          hQnPlusQnMinusStarPtTrkPlus[ibin][iH][iPt]->Fill( (QhfPlus[ibin][iH]*TComplex::Conjugate(QhfMinus[ibin][iH])).Re() );
          hQnMinusQnPlusStarPtTrkPlus[ibin][iH][iPt]->Fill( (QhfMinus[ibin][iH]*TComplex::Conjugate(QhfPlus[ibin][iH])).Re() );
        }
        if(sumwMinusTrk[ibin][iPt]>0 && sumwMinus[ibin]>0 && sumwPlus[ibin]>0){ //tuos
          hQnPlusQnMinusStarPtTrkMinus[ibin][iH][iPt]->Fill( (QhfPlus[ibin][iH]*TComplex::Conjugate(QhfMinus[ibin][iH])).Re() );
          hQnMinusQnPlusStarPtTrkMinus[ibin][iH][iPt]->Fill( (QhfMinus[ibin][iH]*TComplex::Conjugate(QhfPlus[ibin][iH])).Re() );
        }

      }
    }//iH

  //Histograms for recentering
    for(int iH=0; iH<nHarmonics; iH++){
      if(sumwPlus[ibin]>0) {
        hQhfPlusX[ibin][iH]->Fill(QhfPlus[ibin][iH].Re());      
        hQhfPlusY[ibin][iH]->Fill(QhfPlus[ibin][iH].Im());      
      }
      if(sumwMinus[ibin]>0) {
        hQhfMinusX[ibin][iH]->Fill(QhfMinus[ibin][iH].Re());      
        hQhfMinusY[ibin][iH]->Fill(QhfMinus[ibin][iH].Im());      
      }
      if(sumwMinusPlus[ibin]>0) {
        hQhfMinusPlusX[ibin][iH]->Fill(QhfMinusPlus[ibin][iH].Re());      
        hQhfMinusPlusY[ibin][iH]->Fill(QhfMinusPlus[ibin][iH].Im());      
      }
      if(sumwDEMinusPlus[ibin]>0) {
        hQDEMinusPlusX[ibin][iH]->Fill(QDEMinusPlus[ibin][iH].Re());      
        hQDEMinusPlusY[ibin][iH]->Fill(QDEMinusPlus[ibin][iH].Im());      
      }
      for(int iPt=0; iPt<nPtBin; iPt++){
        if(sumwPlusTrk[ibin][iPt]>0){
          hQtrkPtPlusX[ibin][iH][iPt]->Fill(QtrkPtPlus[ibin][iH][iPt].Re());
          hQtrkPtPlusY[ibin][iH][iPt]->Fill(QtrkPtPlus[ibin][iH][iPt].Im());
        }
        if(sumwMinusTrk[ibin][iPt]>0){
          hQtrkPtMinusX[ibin][iH][iPt]->Fill(QtrkPtMinus[ibin][iH][iPt].Re());
          hQtrkPtMinusY[ibin][iH][iPt]->Fill(QtrkPtMinus[ibin][iH][iPt].Im());
        }
        if(sumwMidTrkPt[ibin][iPt]>0) {
          hQtrkPtMidX[ibin][iH][iPt]->Fill(QtrkPtMid[ibin][iH][iPt].Re());
          hQtrkPtMidY[ibin][iH][iPt]->Fill(QtrkPtMid[ibin][iH][iPt].Im());
        }
      }
      if(sumwMidTrk[ibin]>0){
        hQtrkMidX[ibin][iH]->Fill(QtrkMid[ibin][iH].Re());
        hQtrkMidY[ibin][iH]->Fill(QtrkMid[ibin][iH].Im());
      }
    }//iH


#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  // if the SetupData is always needed
  auto setup = iSetup.getData(setupToken_);
  // if need the ESHandle to check if the SetupData was there or not
  auto pSetup = iSetup.getHandle(setupToken_);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void SPAnalysis::beginJob() {
  // please remove this method if not needed

  TH1::SetDefaultSumw2(true);
  TH1::StatOverflows(true);

  for(int iC=0;iC<nCentBin;iC++){
    for(int iH=0; iH<nHarmonics; iH++){
      hQnPlusQnMinusStar[iC][iH] = fs->make<TH1D>(Form("hQnPlusQnMinusStar_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQnMinusQnPlusStar[iC][iH] = fs->make<TH1D>(Form("hQnMinusQnPlusStar_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQnMinusQnMidStar[iC][iH] = fs->make<TH1D>(Form("hQnMinusQnMidStar_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQnPlusQnMidStar[iC][iH] = fs->make<TH1D>(Form("hQnPlusQnMidStar_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQnMinusPlusQnDEStar[iC][iH] = fs->make<TH1D>(Form("hQnMinusPlusQnDEStar_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQnMinusPlusQnMidStar[iC][iH] = fs->make<TH1D>(Form("hQnMinusPlusQnMidStar_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQnDEQnMidStar[iC][iH] = fs->make<TH1D>(Form("hQnDEQnMidStar_iC%d_iH%d",iC,iH),"",10,-1000,1000);
     
      hQhfPlusX[iC][iH] = fs->make<TH1D>(Form("hQhfPlusX_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQhfPlusY[iC][iH] = fs->make<TH1D>(Form("hQhfPlusY_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQhfMinusX[iC][iH] = fs->make<TH1D>(Form("hQhfMinusX_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQhfMinusY[iC][iH] = fs->make<TH1D>(Form("hQhfMinusY_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQhfMinusPlusX[iC][iH] = fs->make<TH1D>(Form("hQhfMinusPlusX_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQhfMinusPlusY[iC][iH] = fs->make<TH1D>(Form("hQhfMinusPlusY_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQDEMinusPlusX[iC][iH] = fs->make<TH1D>(Form("hQDEMinusPlusX_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQDEMinusPlusY[iC][iH] = fs->make<TH1D>(Form("hQDEMinusPlusY_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQtrkMidX[iC][iH] = fs->make<TH1D>(Form("hQtrkMidX_iC%d_iH%d",iC,iH),"",10,-1000,1000);
      hQtrkMidY[iC][iH] = fs->make<TH1D>(Form("hQtrkMidY_iC%d_iH%d",iC,iH),"",10,-1000,1000);

//tuos-sub
      hQnA1QnB1[iC][iH] = fs->make<TH1D>(Form("hQnA1QnB1_iC%d_iH%d",iC,iH),"",10,-9999,9999);
      hQnA2QnB2[iC][iH] = fs->make<TH1D>(Form("hQnA2QnB2_iC%d_iH%d",iC,iH),"",10,-9999,9999);
      hQnA3QnB3[iC][iH] = fs->make<TH1D>(Form("hQnA3QnB3_iC%d_iH%d",iC,iH),"",10,-9999,9999);
      hQnA4QnB4[iC][iH] = fs->make<TH1D>(Form("hQnA4QnB4_iC%d_iH%d",iC,iH),"",10,-9999,9999);
      hQnA4xQnB4x[iC][iH] = fs->make<TH1D>(Form("hQnA4xQnB4x_iC%d_iH%d",iC,iH),"",10,-9999,9999);
      hQnA5QnB5[iC][iH] = fs->make<TH1D>(Form("hQnA5QnB5_iC%d_iH%d",iC,iH),"",10,-9999,9999);
      hQnA6QnB6[iC][iH] = fs->make<TH1D>(Form("hQnA6QnB6_iC%d_iH%d",iC,iH),"",10,-9999,9999);
      hQnA6xQnB6x[iC][iH] = fs->make<TH1D>(Form("hQnA6xQnB6x_iC%d_iH%d",iC,iH),"",10,-9999,9999);


      for(int iPt=0; iPt<nPtBin; iPt++){
        hqnPlusQnMinusStar[iC][iH][iPt] = fs->make<TH1D>(Form("hqnPlusQnMinusStar_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hqnMinusQnPlusStar[iC][iH][iPt] = fs->make<TH1D>(Form("hqnMinusQnPlusStar_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hqnQnPlusStar[iC][iH][iPt] = fs->make<TH1D>(Form("hqnQnPlusStar_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000); 
        hqnQnMinusStar[iC][iH][iPt] = fs->make<TH1D>(Form("hqnQnMinusStar_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000); 
        hqnQnMinusPlusStar[iC][iH][iPt] = fs->make<TH1D>(Form("hqnQnMinusPlusStar_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000); 

        hQnPlusQnMinusStarPtTrkPlus[iC][iH][iPt] = fs->make<TH1D>(Form("hQnPlusQnMinusStarPtTrkPlus_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000); //tuos
        hQnPlusQnMinusStarPtTrkMinus[iC][iH][iPt] = fs->make<TH1D>(Form("hQnPlusQnMinusStarPtTrkMinus_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hQnMinusQnPlusStarPtTrkPlus[iC][iH][iPt] = fs->make<TH1D>(Form("hQnMinusQnPlusStarPtTrkPlus_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hQnMinusQnPlusStarPtTrkMinus[iC][iH][iPt] = fs->make<TH1D>(Form("hQnMinusQnPlusStarPtTrkMinus_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hSumCosn[iC][iH][iPt] = fs->make<TH1D>(Form("hSumCosn_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-10000,10000);
   
        hQtrkPtPlusX[iC][iH][iPt] = fs->make<TH1D>(Form("hQtrkPtPlusX_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hQtrkPtPlusY[iC][iH][iPt] = fs->make<TH1D>(Form("hQtrkPtPlusY_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hQtrkPtMinusX[iC][iH][iPt] = fs->make<TH1D>(Form("hQtrkPtMinusX_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hQtrkPtMinusY[iC][iH][iPt] = fs->make<TH1D>(Form("hQtrkPtMinusY_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hQtrkPtMidX[iC][iH][iPt] = fs->make<TH1D>(Form("hQtrkPtMidX_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
        hQtrkPtMidY[iC][iH][iPt] = fs->make<TH1D>(Form("hQtrkPtMidY_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-1000,1000);
//tuos-sub
        hqnPlusQA1Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnPlusQA1Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnPlusQA2Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnPlusQA2Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnPlusQA3Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnPlusQA3Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnPlusQA4Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnPlusQA4Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnPlusQA4xStar[iC][iH][iPt] = fs->make<TH1D>(Form("hqnPlusQA4xStar_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnPlusQA5Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnPlusQA5Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnPlusQA6Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnPlusQA6Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnPlusQA6xStar[iC][iH][iPt] = fs->make<TH1D>(Form("hqnPlusQA6xStar_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnMinusQB1Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnMinusQB1Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnMinusQB2Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnMinusQB2Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnMinusQB3Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnMinusQB3Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnMinusQB4Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnMinusQB4Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnMinusQB4xStar[iC][iH][iPt] = fs->make<TH1D>(Form("hqnMinusQB4xStar_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnMinusQB5Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnMinusQB5Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnMinusQB6Star[iC][iH][iPt] = fs->make<TH1D>(Form("hqnMinusQB6Star_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hqnMinusQB6xStar[iC][iH][iPt] = fs->make<TH1D>(Form("hqnMinusQB6xStar_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
//tuospt
        hQnA1QnB1ptp[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA1QnB1ptp_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA1QnB1ptm[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA1QnB1ptm_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA2QnB2ptp[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA2QnB2ptp_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA2QnB2ptm[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA2QnB2ptm_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA3QnB3ptp[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA3QnB3ptp_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA3QnB3ptm[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA3QnB3ptm_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA4QnB4ptp[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA4QnB4ptp_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA4QnB4ptm[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA4QnB4ptm_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA4xQnB4xptp[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA4xQnB4xptp_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA4xQnB4xptm[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA4xQnB4xptm_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA5QnB5ptp[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA5QnB5ptp_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA5QnB5ptm[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA5QnB5ptm_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA6QnB6ptp[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA6QnB6ptp_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA6QnB6ptm[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA6QnB6ptm_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA6xQnB6xptp[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA6xQnB6xptp_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);
        hQnA6xQnB6xptm[iC][iH][iPt] = fs->make<TH1D>(Form("hQnA6xQnB6xptm_iC%d_iH%d_iPt%d",iC,iH,iPt),"",10,-9999,9999);

      }
    }//iH
    for(int iPt=0; iPt<nPtBin; iPt++){
      hMeanPt[iC][iPt] = fs->make<TH1D>(Form("hMeanPt_iC%d_iPt%d",iC,iPt),"",20,0,200);
      hMeanPtMid[iC][iPt] = fs->make<TH1D>(Form("hMeanPtMid_iC%d_iPt%d",iC,iPt),"",20,0,200);

      hSumHF3Pt[iC][iPt] = fs->make<TH1D>(Form("hSumHF3Pt_iC%d_iPt%d",iC,iPt),"",20,0,20000);
      hSumHF4Pt[iC][iPt] = fs->make<TH1D>(Form("hSumHF4Pt_iC%d_iPt%d",iC,iPt),"",20,0,20000);
      //hNtrk10Pt[iC][iPt] = fs->make<TH1D>(Form("hNtrk10Pt_iC%d_iPt%d",iC,iPt),"",20,0,20000);
      //hNtrk24Pt[iC][iPt] = fs->make<TH1D>(Form("hNtrk24Pt_iC%d_iPt%d",iC,iPt),"",20,0,20000);
      //hSumTrkPt10Pt[iC][iPt] = fs->make<TH1D>(Form("hSumTrkPt10Pt_iC%d_iPt%d",iC,iPt),"",20,0,20000);
      //hSumTrkPt24Pt[iC][iPt] = fs->make<TH1D>(Form("hSumTrkPt24Pt_iC%d_iPt%d",iC,iPt),"",20,0,20000);

    }
    hPtCent[iC] = fs->make<TH1D>(Form("hPtCent_iC%d",iC),"",500,0,50);
    hPtCentCorr[iC] = fs->make<TH1D>(Form("hPtCentCorr_iC%d",iC),"",500,0,50);
    hSumHF3[iC] = fs->make<TH1D>(Form("hSumHF3_iC%d",iC),"",200,0,20000);
    hSumHF4[iC] = fs->make<TH1D>(Form("hSumHF4_iC%d",iC),"",200,0,20000);
    hNtrk10[iC] = fs->make<TH1D>(Form("hNtrk10_iC%d",iC),"",20,0,20000);
    hNtrk24[iC] = fs->make<TH1D>(Form("hNtrk24_iC%d",iC),"",20,0,20000);
    hNtrkoff[iC] = fs->make<TH1D>(Form("hNtrkoff_iC%d",iC),"",500,0,1000);
    hNtrk24Corr[iC] = fs->make<TH1D>(Form("hNtrk24Corr_iC%d",iC),"",20,0,20000);
    hSumTrkPt10[iC] = fs->make<TH1D>(Form("hSumTrkPt10_iC%d",iC),"",20,0,20000);
    hSumTrkPt24[iC] = fs->make<TH1D>(Form("hSumTrkPt24_iC%d",iC),"",20,0,20000);
  }
  histcent = fs->make<TH1D>(Form("histcent"),"",200,0,200); 
  hNtrkCent = fs->make<TH2D>(Form("hNtrkCent"),"",200,0,200, 1000,0,10000);
  hNtrkCorrCent = fs->make<TH2D>(Form("hNtrkCorrCent"),"",200,0,200, 1000,0,10000);
}

// ------------ method called once each job just after ending the event loop  ------------
void SPAnalysis::endJob() {
  // please remove this method if not needed
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void SPAnalysis::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //edm::ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks", edm::InputTag("ctfWithMaterialTracks"));
  //descriptions.addWithDefaultLabel(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SPAnalysis);
