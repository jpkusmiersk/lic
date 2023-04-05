#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/Muon.h"

#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/Jet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "CommonTools/UtilAlgos/interface/ParameterAdapter.h"

//#include "DataFormats/JetReco/interface/PFJet.h"
//#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <sstream>
#include "cmath"


using namespace std;


//object definition
class Lic : public edm::one::EDAnalyzer<> {
public:

  //constructor, function is called when new object is created
  explicit Lic(const edm::ParameterSet& conf);

  //destructor, function is called when object is destroyed
  ~Lic();

  //edm filter plugin specific functions
  virtual void beginJob();
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

private:

  edm::ParameterSet theConfig;
  unsigned int theEventCount;
  TH1D *histo;
  TH1D *histo1;
  TH1D *histo2;
  TH1D *histo3;
  //TH2D *histo2D;
  //TH2D *histo2D1;

  edm::EDGetTokenT< vector<pat::Muon> > theMuonToken;
  edm::EDGetTokenT<l1t::MuonBxCollection> theGmtToken;
  edm::EDGetTokenT< vector<pat::Jet> > theJetToken;
  edm::EDGetTokenT<l1t::JetBxCollection> theGjtToken;

};


Lic::Lic(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;
  theMuonToken = consumes< vector<pat::Muon> >( theConfig.getParameter<edm::InputTag>("muonSrc"));
  theGmtToken  = consumes<l1t::MuonBxCollection>( theConfig.getParameter<edm::InputTag>("gmtSrc"));
  theJetToken = consumes< vector<pat::Jet> >( theConfig.getParameter<edm::InputTag>("jetSrc"));
  theGjtToken  = consumes<l1t::JetBxCollection>( theConfig.getParameter<edm::InputTag>("gjtSrc"));
}

Lic::~Lic()
{
  cout <<" DTOR" << endl;
}

void Lic::beginJob()
{
  //create a histogram
  histo =new TH1D("histo","test; #GMT; #events",100, -6., 6.);
  histo1 =new TH1D("histo1","test; #GMT; #events",100, -6., 6.);
  histo2 =new TH1D("histo2","test; #GMT; #events",100, -6., 6.);
  histo3 =new TH1D("histo3","test; #GMT; #events",100, -6., 6.);
  //histo1 =new TH1D("histo1","test; #GMT; #events",100, -1., 12.);
  //histo2 =new TH1D("histo2","test; #GMT; #events",100, 0., 20.);
  //histo2D = new TH2D("histo2D","y,x,#entries", 100, 0., 1000., 100, 0., 4.);
  //histo2D1 = new TH2D("histo2D1","y,x,#entries", 100, -4., 4, 100, -4., 4.);
  cout << "HERE Lic::beginJob()" << endl;
}

void Lic::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");
  //write histogram data
  histo->Write();
  histo1->Write();
  histo2->Write();
  histo3->Write();
  //histo2D->Write();
  //histo2D1->Write();
  myRootFile.Close();
  delete histo;
  cout << "HERE Cwiczenie::endJob()" << endl;
}


void Lic::analyze(
    const edm::Event& ev, const edm::EventSetup& es)
{
  int wiadomosci = 100;
  //vector<double> x;
  //vector<double> y;
  std::cout << " -------------------------------- HERE Cwiczenie::analyze "<< std::endl;
  bool debug = true;
  const vector<pat::Muon> & muons = ev.get(theMuonToken);
  //const vector<l1t::MuonBxCollection> & muons = ev.get(theGmtToken);
  //if (debug && theEventCount%wiadomosci==0) std::cout <<" number of muons: " << muons.size() <<std::endl;
  for (const auto & muon : muons) {
    //if (debug) std::cout <<" reco muon pt: "<<muon.pt()<<std::endl;
    //if (muon.phi()>1 || muon.phi()<-1){
      //histo->Fill(muon.phi());
      //histo1->Fill(muon.energy());
      //histo2->Fill(muon.p());
      //h_2dgaus->Fill(muon.eta(),0.);
      //x.push_back(muon.eta());
    //}
    //histo->Fill(muon.pt());
    //histo1->Fill(muon.vy());
    //histo2->Fill(muon.vz());
  }

  const l1t::MuonBxCollection & gmts = ev.get(theGmtToken); 
  int bxNumber = 0;
  
  for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
    if (theEventCount%wiadomosci==0) {
    //std::cout <<"GMT: "<<it->phi()<<std::endl;
    }
    //if (it->phi()>1 || it->phi()<-1){
        //histo1->Fill(it->phi());
      //h_2dgaus->Fill(0.,it->eta()); 
      //y.push_back(it->eta());
    //}
    //histo->Fill(it->hwEta());
  }
  //for(long unsigned int i = 1; i<y.size(); i++){
    //h_2dgaus->Fill(x[i],y[i]);
  //}
  /////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////
  //Tu robione sa Jesty
  const vector<pat::Jet> & jets = ev.get(theJetToken);
  //const vector<l1t::MuonBxCollection> & muons = ev.get(theGmtToken);
  if (debug && theEventCount%wiadomosci==0) std::cout <<" number of muons: " << jets.size() <<std::endl;
  for (const auto & jet : jets) {
    if (debug && theEventCount%wiadomosci==0) std::cout <<" reco jet pt: "<<jet.pt()<<std::endl;
    if (jet.pt()>10){
      histo->Fill(jet.eta());
      //histo1->Fill(muon.energy());
      //histo2->Fill(muon.p());
      //h_2dgaus->Fill(muon.eta(),0.);
      //x.push_back(muon.eta());
    }
    if(jet.pt()>20){
      histo1->Fill(jet.eta());
    }
    if(jet.pt()>30){
      histo2->Fill(jet.eta());
    }
    if(jet.pt()>40){
      histo3->Fill(jet.eta());
    }
    //histo->Fill(muon.pt());
    //histo1->Fill(muon.vy());
    //histo2->Fill(muon.vz());
  }

  const l1t::JetBxCollection & gjts = ev.get(theGjtToken); 
  int bx1Number = 0;
  
  for (l1t::JetBxCollection::const_iterator it = gjts.begin(bx1Number); it != gjts.end(bx1Number); ++it) {
    if (theEventCount%wiadomosci==0) {
      std::cout <<"GJT: "<<it->pt()<<std::endl;
    }
    //if (it->phi()>1 || it->phi()<-1){
        //histo1->Fill(it->phi());
      //h_2dgaus->Fill(0.,it->eta()); 
      //y.push_back(it->eta());
    //}
    //histo->Fill(it->pt());
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //TUTAJ HISTOGRAMY 2D
  //const l1t::JetBxCollection & gjts = ev.get(theGjtToken); 
  //int bx1Number = 0;
  
  //for (l1t::JetBxCollection::const_iterator it = gjts.begin(bx1Number); it != gjts.end(bx1Number); ++it) {
  /*
  for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
    if (theEventCount%wiadomosci==0) {
      std::cout <<"GJT: "<<it->eta()<<std::endl;
    }
    double eta = it->eta();
    double deltaetamin = 1000;
    double etamin = 0;
    //const vector<pat::Jet> & jets = ev.get(theJetToken);
    //const vector<l1t::MuonBxCollection> & muons = ev.get(theGmtToken);
    //if (debug && theEventCount%wiadomosci==0) std::cout <<" number of muons: " << jets.size() <<std::endl;
    //for (const auto & jet : jets) {
    const vector<pat::Muon> & muons = ev.get(theMuonToken);
    //const vector<l1t::MuonBxCollection> & muons = ev.get(theGmtToken);
    if (debug && theEventCount%wiadomosci==0) std::cout <<" number of muons: " << muons.size() <<std::endl;
    for (const auto & muon : muons) {
      if (debug && theEventCount%wiadomosci==0) std::cout <<" reco jet pt: "<<muon.pt()<<std::endl;
      double etareco = muon.eta();
      double deltaeta = abs(eta - etareco);
      if(deltaeta<deltaetamin){
        etamin = etareco;
        deltaetamin = deltaeta;
      } 
      histo2D -> Fill(eta,etareco);
      
      
    }
    histo2D1 -> Fill(eta,etamin);
  
  }*/
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //TUTAJ DELTA R
  for (const auto & muon : muons) {
    double DeltaRmin = 100000;
    for (const auto & jet : jets) {
      if(muon.pt()>5){
        double DeltaR = reco::deltaR(muon, jet);
        if(DeltaR<DeltaRmin){
          DeltaRmin = DeltaR;
        }
        //histo -> Fill(DeltaR);
      }
    }
    //histo -> Fill(DeltaRmin);
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*for (const auto & jet : jets) {
    double etajet = jet.eta();
    //double deltaetamin = 10000;
    double jetmin = 0;
    double DeltaRmin = 100000;
    for (const auto & muon : muons) {
      if(muon.pt()>5){
        double etamuon = muon.eta();
        double DeltaR = reco::deltaR(muon, jet);
        if(DeltaR<DeltaRmin){
          DeltaRmin = DeltaR;
          jetmin = etamuon;
        }
    
        //histo2D -> Fill(etamuon,etajet);
      }
    }
    //histo2D -> Fill(etajet,jetmin);
  }
  */
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (const auto & jet : jets) {
    double DeltaRmin = 100000;
    double jetpt = jet.pt();
    for (const auto & muon : muons) {
      double DeltaR = reco::deltaR(muon, jet);
      
      if(DeltaR<DeltaRmin){
        DeltaRmin = DeltaR;
      }
      if(muon.pt()>5){
        //histo2D -> Fill(jetpt,DeltaR);
      }
    }
    if(jet.pt()>2){
      //histo2D -> Fill(jetpt, DeltaRmin);
    }
    //histo2D -> Fill(jetpt,DeltaRmin);
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Tu liczymy mase niezmiennicza
  /*if (jets.size()==4){
    double sumpx = 0;
    double sumpy = 0;
    double sumpz = 0;
    double sumenergy = 0;
    double sumpt = 0;
    for (const auto & jet : jets) {
      sumpx += jet.px();
      sumpy += jet.py();
      sumpz += jet.pz();
      sumenergy += jet.energy();
      sumpt += jet.pt();
    }
    if(sumpt>20){
      histo->Fill(sqrt(sumenergy*sumenergy-sumpx*sumpx-sumpy*sumpy-sumpz*sumpz));
    }
  }
  */
  

  ++theEventCount;
  if (theEventCount%wiadomosci==0) {
  //write std io
  cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" analysed event count:"<<theEventCount << endl;
  }
}

DEFINE_FWK_MODULE(Lic);

