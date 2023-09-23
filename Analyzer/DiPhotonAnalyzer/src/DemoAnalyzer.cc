// Skeleton of this code was taken from: https://twiki.cern.ch/twiki/pub/CMSPublic/WorkBookWriteFrameworkModule/DemoAnalyzer_cmssw1200_v1.cc

// system include files
#include <memory>
#include <cmath>
#include <algorithm> 
#include <utility>  

//----Extra
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/ESHandle.h" 
#include "FWCore/ServiceRegistry/interface/Service.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

// for Root histograms
#include "TH1.h"
#include "TTree.h"
#include "TH2.h"
#include "TLorentzVector.h" //

//classes to extract Muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include<vector>

//classes to extract Photon information
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "EgammaAnalysis/ElectronTools/interface/PFIsolationEstimator.h"

//classes to extract PFJet information
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/BTauReco/interface/JetTag.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/SimpleJetCorrectionUncertainty.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfo.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourInfoMatching.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"




//
// class declaration
//

class DemoAnalyzer : public edm::EDAnalyzer {
  public:
    explicit DemoAnalyzer(const edm::ParameterSet&);
    ~DemoAnalyzer();

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


  private:
    virtual void beginJob() ;
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob() ;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void endRun(edm::Run const&, edm::EventSetup const&);
    virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
    virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;



    //declare histograms
    TH1D *hist1;
    TH1D *hist2;


      // ----------member data ---------------------------

  //photon information
  int numphoton; 
  std::vector<float> photon_e; //energy
  std::vector<float> photon_pt; //transverse momentum
  std::vector<float> photon_px; //momentum in x direction
  std::vector<float> photon_py; //momentum in y direction
  std::vector<float> photon_pz; //momentum in z direction
  std::vector<float> photon_eta; //eta angle
  std::vector<float> photon_phi; //phi angle
  std::vector<float> photon_ch; // charge
  std::vector<float> photon_ET; // transverse energy

  /// jets information
  int numjet; //number of jets
  std::vector<float> jet_e; //energy
  std::vector<float> jet_pt; //transverse momentum
  std::vector<float> jet_px; //momentum in x direction
  std::vector<float> jet_py; //momentum in y direction
  std::vector<float> jet_pz; //momentum in z direction
  std::vector<float> jet_eta; //eta angle
  std::vector<float> jet_phi; //phi angle

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
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
   edm::Service<TFileService> fs; 
  
  //histogram with 100 bins in range 0 to 600
   hist1 = fs->make<TH1D>("PHOTONS", "DiPhoton Mass", 100, 0., 600.);
   hist1->GetXaxis()->SetTitle("Invariant Mass (in GeV/c^2)");
   hist1->GetYaxis()->SetTitle("Number of Events");

  //histogram with 100 bins in range 105 to 160
   hist2 = fs->make<TH1D>("PHOTONS_lower_range", "DiPhoton Mass", 100, 105., 160.);
   hist2->GetXaxis()->SetTitle("Invariant Mass (in GeV/c^2)");
   hist2->GetYaxis()->SetTitle("Number of Events");

}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//





// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
  using namespace edm;
  using namespace std;
  
  // clear containers
  photon_e.clear();
  numphoton = 0;
  photon_pt.clear();
  photon_px.clear();
  photon_py.clear();
  photon_pz.clear();
  photon_eta.clear();
  photon_phi.clear();
  photon_ch.clear();


  numjet = 0;
  jet_e.clear();
  jet_pt.clear();
  jet_px.clear();
  jet_py.clear();
  jet_pz.clear();
  jet_eta.clear();
  jet_phi.clear();

  bool jetConditionsMet = false;

  /// import jets from event
  Handle<reco::PFJetCollection> myjets;
  iEvent.getByLabel("ak5PFJets", myjets);

  // indices for sorting jets
  int jet1_index = 0;
  int jet2_index = 1;

  
  numjet = myjets->size(); // number of jets in event

  // if there are at least two jets in event, choose two with highest pt
  if (myjets.isValid() && numjet >= 2)
  {
    for (reco::PFJetCollection::const_iterator itjet = myjets->begin(); itjet != myjets->end(); ++itjet)
    {
      //update jets containers
      jet_e.push_back(itjet->energy());
      jet_pt.push_back(itjet->pt());
      jet_px.push_back(itjet->px());
      jet_py.push_back(itjet->py());
      jet_pz.push_back(itjet->pz());
      jet_eta.push_back(itjet->eta());
      jet_phi.push_back(itjet->phi());
    }

    // initiate variables
    double jet1_pT = 0;
    double jet2_pT = 0;

    double dijet_invariant_mass = 0;

    double jet1_E = 0;
    double jet2_E = 0;

    double jet1_px = 0;
    double jet1_py = 0;
    double jet1_pz = 0;
    double jet2_px = 0;
    double jet2_py = 0;
    double jet2_pz = 0;

    bool jet1_index_met = false;
    bool jet2_index_met = false;

    //iterate over jets and choose two with highest pt
    for (int i = 0; i < numjet; i++)
    {
      double current_jet = jet_pt.at(i);
      double current_jet_eta = jet_eta.at(i);

      //only take jets within abs(eta)<4.7 range
      if (abs(current_jet_eta) < 4.7)
      {
        if (current_jet > jet1_pT)
        {
          jet1_pT = current_jet;
          jet1_index = i;
          jet1_index_met = true;
        }
        else if (current_jet > jet2_pT )
        {
          jet2_pT = current_jet;
          jet2_index = i;
          jet2_index_met = true;
        }
      }
    }

    // update variables with two chosen jets information
    jet1_E = jet_e.at(jet1_index);
    jet2_E = jet_e.at(jet2_index);

    jet1_px = jet_px.at(jet1_index);
    jet1_py = jet_py.at(jet1_index);
    jet1_pz = jet_pz.at(jet1_index);

    jet2_px = jet_px.at(jet2_index);
    jet2_py = jet_py.at(jet2_index);
    jet2_pz = jet_pz.at(jet2_index);

    // calculate invariant mass if two jets were chosen
    if (jet1_index_met == true && jet2_index_met == true)
    {
      dijet_invariant_mass = (jet1_E + jet2_E) * (jet1_E + jet2_E) - (jet1_px + jet2_px) * (jet1_px + jet2_px) - (jet1_py + jet2_py) * (jet1_py + jet2_py) - (jet1_pz + jet2_pz) * (jet1_pz + jet2_pz);
    }

    // leading and subleading jets must have pt higher than 30 and 20 respectively
    if (jet_pt.at(jet1_index) > 30 && jet_pt.at(jet2_index) > 20)
    {
      // dijet invariant mass must be higher than 250
      if (dijet_invariant_mass > 250)
      {
        // jets eta separation must be greater than 3.5
        if (abs(jet_eta.at(jet1_index) - jet_eta.at(jet2_index)) > 3.5)
        {

          jetConditionsMet = true;
        }
        else return;
      }
    }
  }
  
  // if conditions regarding jets are not met, break the loop and go to next event
  if(jetConditionsMet == false){
    return;
  }



  // import photons information
  Handle<reco::PhotonCollection> myphotons;
  iEvent.getByLabel("photons", myphotons);


  numphoton = myphotons->size(); // number of photons in event

  // update photon information only if there are at least 2 photons in event
  if (myphotons.isValid() && numphoton >= 2){
    for (reco::PhotonCollection::const_iterator itphoton = myphotons->begin(); itphoton != myphotons->end(); ++itphoton)
    {
      photon_e.push_back(itphoton->energy());
      photon_pt.push_back(itphoton->pt());
      photon_px.push_back(itphoton->px());
      photon_py.push_back(itphoton->py());
      photon_pz.push_back(itphoton->pz());
      photon_eta.push_back(itphoton->eta());
      photon_phi.push_back(itphoton->phi());
      photon_ch.push_back(itphoton->charge());
    }
  }
  else return; // if there are less than 2 photons, break the loop

  // variables for sorting photons
  double max_pT = 0;
  double second_pT = 0;
  int first_pT_index = 0;
  int second_pT_index = 0;

  // proceed only if there are at least 2 photons
  if (numphoton >= 2)
  {
    //iterate over photons in event
    for (int i = 0; i < numphoton; i++)
    {
      double curr_photon = photon_pt.at(i);
      double curr_photon_eta = abs(photon_eta.at(i));

      // take only photons within fiducial eta range
      if ( curr_photon_eta < 2.5){
        if(curr_photon_eta > 1.57 || curr_photon_eta < 1.44){


          if (curr_photon > max_pT)
          {
            max_pT = curr_photon;
            first_pT_index = i;
          }

          else if (curr_photon > second_pT)
          {
            second_pT = curr_photon;
            second_pT_index = i;
          }
        }
      }
    }

    // renaming
    int aa = first_pT_index;
    int bb = second_pT_index;

    // define variables for photons
    double Energy1, Energy2;
    double p1x, p1y, p1z;
    double p2x, p2y, p2z;

    // initially set high values for invariant masses
    double s1 = 9999.99;
    double s2 = 9999.99;
    double s = 9999.99;

    // leading and subleading photons information
    Energy1 = photon_e.at(aa);
    Energy2 = photon_e.at(bb);

    p1x = photon_px.at(aa);
    p1y = photon_py.at(aa);
    p1z = photon_pz.at(aa);

    p2x = photon_px.at(bb);
    p2y = photon_py.at(bb);
    p2z = photon_pz.at(bb);

    // calculate diphoton invariant mass
    s1 = (Energy1 + Energy2) * (Energy1 + Energy2);
    s2 = (p1x + p2x) * (p1x + p2x) + (p1y + p2y) * (p1y + p2y) + (p1z + p2z) * (p1z + p2z);
    s = sqrt(s1 - s2);

    // We are not interrested in values below about 105
    if (s < 85)
    {
      return;
    }

    // calculate diphoton eta angle
    double E_tot = Energy1 + Energy2;
    double pz_tot = p1z + p2z;
    double eta_diphoton = 0.5 * log((E_tot + pz_tot) / (E_tot - pz_tot));
    //double eta_diphoton = (photon_eta.at(aa) + photon_eta.at(bb))/2;

    // calculate dijet eta angle
    double eta_dijet = (jet_eta.at(jet1_index) + jet_eta.at(jet2_index)) / 2;

    // difference in eta between diphoton and dijet must be greater than 2.5
    if (abs(eta_diphoton - eta_dijet) >= 2.5)
    {
      return;
    }

    // leading and subleading photons pt thresholds
    if (max_pT > s/3 && second_pT > s/4)
    {
    std::cout << "INVARIANT MASS = " << s << std::endl;

    // if all conditions passed, fill histogram with diphoton invariant mass
    hist1->Fill(s);
    hist2->Fill(s);
    }
  }

///////////////////////////////////////////////////////////////
#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
DemoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
DemoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
DemoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
DemoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
