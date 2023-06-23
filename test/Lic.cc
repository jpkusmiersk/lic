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
#include "TEfficiency.h"


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
  TH1D *jetpt;
  TH1D *jetpt_fullrange;
  TH1D *muonpt;
  TH1D *muonpt_fullrange;
  TH1D *l1tpt;
  TH1D *l1tpt_fullrange;
  TH1D *jetdphi;
  TH1D *jeteta10;
  TH1D *jeteta50;
  TH1D *jeteta100;
  TH1D *jeteta200;
  TH1D *deltar;
  TH1D *deltarmin;
  TH1D *deltarpik;
  TH1D *deltarml1t;
  TH1D *deltarml1tpik;
  TH1D *deltarml1t20;
  TH1D *deltarml1tpik20;
  TH1D *masan2m5_15_50;
  TH1D *masan2m5_15_100;
  TH1D *masan2m5_15_200;
  TH1D *masan2m0_12;
  TH1D *masan2m0_120;
  TH1D *vertex;
  TH1D *jetcount50;
  TH1D *jetcount100;
  TH1D *jetcount200;
  TH1D *muoncount5;
  TH1D *muoncount10;
  TH1D *muoncount20;
  TH1D *histo1;
  TH1D *histo2;
  TH1D *histo3;
  TH1D *histo4;
  TH1D *histo5;
  TH2D *histo2D;
  TH2D *muonvsjeteta;
  TH2D *muonvsjetetamin;
  TH2D *deltaphivspt;
  TH2D *deltarvspt;
  TH2D *deltarvsptmin;
  //TH2D *histo2D1;
  TEfficiency *peta;
  TEfficiency *peta2;
  TEfficiency *peta3;
  TEfficiency *peta4;
  TEfficiency *peta5;
  TEfficiency *peta6;
  TEfficiency *pjeta5;
  TEfficiency *pjeta10;
  TEfficiency *pjeta20;
  TEfficiency *pjeta50;
  TEfficiency *pjeta100;
  TEfficiency *pjeta200;
  TEfficiency *pjeta200_24;
  TEfficiency *pjeta1_5;
  TEfficiency *pjeta1_10;
  TEfficiency *pjeta1_20;
  TEfficiency *pjeta1_50;
  TEfficiency *pjeta1_100;
  TEfficiency *pjeta1_200;
  TEfficiency *pjeta1_200_24;
  TEfficiency *pmeta50;
  TEfficiency *pmeta100;
  TEfficiency *pmeta200;
  TEfficiency *pmeta1_50;
  TEfficiency *pmeta1_100;
  TEfficiency *pmeta1_200;
  TEfficiency *pml1t5;
  TEfficiency *pl1tm5;
  TEfficiency *pjetamuon50;
  TEfficiency *pjetal1t50;

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
  histo =new TH1D("histo","test; #GMT; #events",10, 0., 10.);
  jetpt =new TH1D("jetpt","test; #GMT; #events",1000, 0., 1000.);
  jetpt_fullrange =new TH1D("jetpt_fullrange","test; #GMT; #events",70, 0., 7000.);
  muonpt =new TH1D("muonpt","test; #GMT; #events",200, 0., 200.);
  muonpt_fullrange =new TH1D("muonpt_fullrange","test; #GMT; #events",100, 0., 10000.);
  l1tpt =new TH1D("l1tpt","test; #GMT; #events",200, 0., 200.);
  l1tpt_fullrange =new TH1D("l1tpt_fullrange","test; #GMT; #events",300, 0., 300.);
  jetdphi =new TH1D("jetdphi","test; #GMT; #events",100, 0., 3.5);
  jeteta10 =new TH1D("jeteta10","test; #GMT; #events",200, -6., 6.);
  jeteta50 =new TH1D("jeteta50","test; #GMT; #events",200, -6., 6.);
  jeteta100 =new TH1D("jeteta100","test; #GMT; #events",200, -6., 6.);
  jeteta200 =new TH1D("jeteta200","test; #GMT; #events",200, -6., 6.);
  deltar =new TH1D("deltar","test; #GMT; #events",100, 0., 8.);
  deltarmin =new TH1D("deltarmin","test; #GMT; #events",100, 0., 6.);
  deltarpik =new TH1D("deltarpik","test; #GMT; #events",100, 0., 4.);
  deltarml1t =new TH1D("deltarml1t","test; #GMT; #events",100, 0., 8.);
  deltarml1tpik =new TH1D("deltarml1tpik","test; #GMT; #events",100, 0., 4.);
  deltarml1t20 =new TH1D("deltarml1t20","test; #GMT; #events",100, 0., 8.);
  deltarml1tpik20 =new TH1D("deltarml1tpik20","test; #GMT; #events",100, 0., 4.);
  masan2m5_15_50 =new TH1D("masan2m5_15_50","test; #GMT; #events",75, 0., 1.5);
  masan2m5_15_100 =new TH1D("masan2m5_15_100","test; #GMT; #events",150, 0., 1.5);
  masan2m5_15_200 =new TH1D("masan2m5_15_200","test; #GMT; #events",300, 0., 1.5);
  masan2m0_12 =new TH1D("masan2m0_12","test; #GMT; #events",240, 0., 12.);
  masan2m0_120 =new TH1D("masan2m0_120","test; #GMT; #events",60., 0., 120.);
  vertex =new TH1D("vertex","test; #GMT; #events",240, 0.01, 12.);
  jetcount50 =new TH1D("jetcount50","test; #GMT; #events",10, 0., 10.);
  jetcount100 =new TH1D("jetcount100","test; #GMT; #events",10, 0., 10.);
  jetcount200 =new TH1D("jetcount200","test; #GMT; #events",10, 0., 10.);
  muoncount5 =new TH1D("muoncount5","test; #GMT; #events",10, 0., 10.);
  muoncount10 =new TH1D("muoncount10","test; #GMT; #events",10, 0., 10.);
  muoncount20 =new TH1D("muoncount20","test; #GMT; #events",10, 0., 10.);
  histo1 =new TH1D("histo1","test; #GMT; #events",10, 0., 10.);
  histo2 =new TH1D("histo2","test; #GMT; #events",10, 0., 10.);
  histo3 =new TH1D("histo3","test; #GMT; #events",10, 0., 10.);
  histo4 =new TH1D("histo4","test; #GMT; #events",10, 0., 10.);
  histo5 =new TH1D("histo5","test; #GMT; #events",10, 0., 10.);
  //histo1 =new TH1D("histo1","test; #GMT; #events",100, -1., 12.);
  //histo2 =new TH1D("histo2","test; #GMT; #events",100, 0., 20.);
  histo2D = new TH2D("histo2D","y,x,#entries", 100, 0., 1000., 100, 0., 4.);
  muonvsjeteta = new TH2D("muonvsjeteta","y,x,#entries", 100, -5., 5., 100, -4., 4.);
  muonvsjetetamin = new TH2D("muonvsjetetamin","y,x,#entries", 100, -5., 5., 100, -4., 4.);
  deltaphivspt = new TH2D("deltaphivspt","y,x,#entries", 100, 0., 1000., 100, 0., 3.5);
  deltarvspt = new TH2D("deltarvspt","y,x,#entries", 100, 0., 1000., 100, 0., 4.);
  deltarvsptmin = new TH2D("deltarvsptmin","y,x,#entries", 100, 0., 1000., 100, 0., 4.);
  //histo2D1 = new TH2D("histo2D1","y,x,#entries", 100, -4., 4, 100, -4., 4.);
  peta = new TEfficiency("peta", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta2 = new TEfficiency("peta2", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta3 = new TEfficiency("peta3", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta4 = new TEfficiency("peta4", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta5 = new TEfficiency("peta5", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta6 = new TEfficiency("peta6", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta5 = new TEfficiency("pjeta5", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta10 = new TEfficiency("pjeta10", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta20 = new TEfficiency("pjeta20", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta1_5 = new TEfficiency("pjeta1_5", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta1_10 = new TEfficiency("pjeta1_10", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta1_20 = new TEfficiency("pjeta1_20", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta50 = new TEfficiency("pjeta50", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta100 = new TEfficiency("pjeta100", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta200 = new TEfficiency("pjeta200", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta200_24 = new TEfficiency("pjeta200_24", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta1_50 = new TEfficiency("pjeta1_50", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta1_100 = new TEfficiency("pjeta1_100", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta1_200 = new TEfficiency("pjeta1_200", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjeta1_200_24 = new TEfficiency("pjeta1_200_24", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pmeta50 = new TEfficiency("pmeta50", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pmeta100 = new TEfficiency("pmeta100", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pmeta200 = new TEfficiency("pmeta200", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pmeta1_50 = new TEfficiency("pmeta1_50", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pmeta1_100 = new TEfficiency("pmeta1_100", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pmeta1_200 = new TEfficiency("pmeta1_200", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pml1t5 = new TEfficiency("pml1t5", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pl1tm5 = new TEfficiency("pl1tm5", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjetamuon50 = new TEfficiency("pjetamuon50", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  pjetal1t50 = new TEfficiency("pjetal1t50", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  cout << "HERE Lic::beginJob()" << endl;
}

void Lic::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");
  //write histogram data
  histo->Write();
  jetpt->Write();
  jetpt_fullrange->Write();
  muonpt->Write();
  muonpt_fullrange->Write();
  l1tpt->Write();
  l1tpt_fullrange->Write();
  jetdphi->Write();
  jeteta10->Write();
  jeteta50->Write();
  jeteta100->Write();
  jeteta200->Write();
  deltar->Write();
  deltarmin->Write();
  deltarpik->Write();
  deltarml1t->Write();
  deltarml1tpik->Write();
  deltarml1t20->Write();
  deltarml1tpik20->Write();
  masan2m5_15_50->Write();
  masan2m5_15_100->Write();
  masan2m5_15_200->Write();
  masan2m0_12->Write();
  masan2m0_120->Write();
  vertex->Write();
  jetcount50->Write();
  jetcount100->Write();
  jetcount200->Write();
  muoncount5->Write();
  muoncount10->Write();
  muoncount20->Write();
  histo1->Write();
  histo2->Write();
  histo3->Write();
  histo4->Write();
  histo5->Write();
  histo2D->Write();
  muonvsjeteta->Write();
  muonvsjetetamin->Write();
  deltaphivspt->Write();
  deltarvspt->Write();
  deltarvsptmin->Write();
  peta->Write();
  peta2->Write();
  peta3->Write();
  peta4->Write();
  peta5->Write();
  peta6->Write();
  pjeta5->Write();
  pjeta10->Write();
  pjeta20->Write();
  pjeta50->Write();
  pjeta100->Write();
  pjeta200->Write();
  pjeta200_24->Write();
  pjeta1_5->Write();
  pjeta1_10->Write();
  pjeta1_20->Write();
  pjeta1_50->Write();
  pjeta1_100->Write();
  pjeta1_200->Write();
  pjeta1_200_24->Write();
  pmeta50->Write();
  pmeta100->Write();
  pmeta200->Write();
  pmeta1_50->Write();
  pmeta1_100->Write();
  pmeta1_200->Write();
  pml1t5->Write();
  pl1tm5->Write();
  pjetamuon50->Write();
  pjetal1t50->Write();
  //histo2D1->Write();
  myRootFile.Close();
  delete histo;
  cout << "HERE Cwiczenie::endJob()" << endl;
}


void Lic::analyze(const edm::Event& ev, const edm::EventSetup& es){
  int wiadomosci = 100;
  //vector<double> x;
  //vector<double> y;
  //std::cout << " -------------------------------- HERE Cwiczenie::analyze "<< std::endl;
  bool debug = true;
  const vector<pat::Muon> & muons = ev.get(theMuonToken);
  //const vector<l1t::MuonBxCollection> & muons = ev.get(theGmtToken);
  //if (debug && theEventCount%wiadomosci==0) std::cout <<" number of muons: " << muons.size() <<std::endl;
  for (const auto & muon : muons) {
    //if (debug) std::cout <<" reco muon pt: "<<muon.pt()<<std::endl;
    //histo -> Fill(muons.size());
    muonpt->Fill(muon.pt());
    muonpt_fullrange->Fill(muon.pt());
    if (muon.pt()>5){
      //muonpt->Fill(muon.py());
      //histo -> Fill(muons.size());
      //histo->Fill(muon.pt());
      //histo1->Fill(muon.energy());
      //histo2->Fill(muon.p());
      //h_2dgaus->Fill(muon.eta(),0.);
      //x.push_back(muon.eta());
    }
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
    l1tpt->Fill(it->pt());
    l1tpt_fullrange->Fill(it->pt());
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
  //if (debug && theEventCount%wiadomosci==0) std::cout <<" number of muons: " << jets.size() <<std::endl;
  
  for (const auto & jet : jets) {
    //cout<<"bat"<<endl;
    //if (debug && theEventCount%wiadomosci==0) std::cout <<" reco jet pt: "<<jet.pt()<<std::endl;
    jetpt->Fill(jet.pt());
    jetpt_fullrange->Fill(jet.pt());
    //if (jet.pt()>50){
      //jetpt->Fill(jet.pt());
      //histo->Fill(jets.size());
      //histo1->Fill(muon.energy());
      //histo2->Fill(muon.p());
      //h_2dgaus->Fill(muon.eta(),0.);
      //x.push_back(muon.eta());
    //}
    if(jet.pt()>10){
      jeteta10->Fill(jet.eta());
    }
    if(jet.pt()>50){
      jeteta50->Fill(jet.eta());
    }
    if(jet.pt()>100){
      jeteta100->Fill(jet.eta());
    }
    if(jet.pt()>200){
      jeteta200->Fill(jet.eta());
    }
    //histo->Fill(muon.pt());
    //histo1->Fill(muon.vy());
    //histo2->Fill(muon.vz());
  }

  const l1t::JetBxCollection & gjts = ev.get(theGjtToken); 
  int bx1Number = 0;
  
  for (l1t::JetBxCollection::const_iterator it = gjts.begin(bx1Number); it != gjts.end(bx1Number); ++it) {
    if (theEventCount%wiadomosci==0) {
      //std::cout <<"GJT: "<<it->pt()<<std::endl;
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
        deltar -> Fill(DeltaR);
      }
    }
    deltarmin -> Fill(DeltaRmin);
    deltarpik->Fill(DeltaRmin);
    /*
    double DeltaRminml1t  = 10000;
    for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it){
      if(muon.pt()>5){
        double DeltaRml1t = reco::deltaR(muon,*it);
        if(DeltaRml1t<DeltaRminml1t){
          DeltaRminml1t = DeltaRml1t;
        }
        deltarml1t->Fill(DeltaRml1t);
      }
    }
    deltarml1tpik->Fill(DeltaRminml1t);
    */
  }
  /*
  for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it){
    double DeltaRminml1t  = 10000;
    for (const auto & jet : jets) {
      if(it->pt()>5){
        double DeltaRml1t = reco::deltaR(*it, jet);
        if(DeltaRml1t<DeltaRminml1t){
          DeltaRminml1t = DeltaRml1t;
        }
        deltarml1t->Fill(DeltaRml1t);
      }
    }
    deltarml1tpik->Fill(DeltaRminml1t);
  }
  */
  for (const auto & jet : jets){
    if (jet.pt()>50){
      double DeltaRminml1t  = 10000;
      for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
        if(it->pt()>5){
          double DeltaRml1t = reco::deltaR(*it, jet);
          if(DeltaRml1t<DeltaRminml1t){
            DeltaRminml1t = DeltaRml1t;
          }
          deltarml1t->Fill(DeltaRml1t);
        }
      }
      deltarml1tpik->Fill(DeltaRminml1t);
    }
  }
  for (const auto & jet : jets){
    if (jet.pt()>50){
      double DeltaRminml1t20  = 10000;
      for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
        if(it->pt()>22){
          double DeltaRml1t20 = reco::deltaR(*it, jet);
          if(DeltaRml1t20<DeltaRminml1t20){
            DeltaRminml1t20 = DeltaRml1t20;
          }
          deltarml1t20->Fill(DeltaRml1t20);
        }
      }
      deltarml1tpik20->Fill(DeltaRminml1t20);
    }
  }
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  for (const auto & jet : jets) {
    double etajet = jet.eta();
    //double deltaetamin = 10000;
    double jetmin = 0;
    double DeltaRmin = 100000;
    for (const auto & muon : muons) {
      if(muon.pt()>5 && muon.isMediumMuon()){
        double etamuon = muon.eta();
        double DeltaR = reco::deltaR(muon, jet);
        if(DeltaR<DeltaRmin){
          DeltaRmin = DeltaR;
          jetmin = etamuon;
        }
    
        muonvsjeteta -> Fill(etamuon,etajet);
      }
    }
    muonvsjetetamin -> Fill(jetmin,etajet);
  }
  
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
        deltarvspt -> Fill(jetpt,DeltaR);
      }
    }
    if(jet.pt()>2){
      //histo2D -> Fill(jetpt, DeltaRmin);
    }
    deltarvsptmin -> Fill(jetpt,DeltaRmin);
  }
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Tu liczymy mase niezmiennicza
  
  if (muons.size()==2){
    double sumpx = 0;
    double sumpy = 0;
    double sumpz = 0;
    double sumenergy = 0;
    double sumpt = 0;
    int ladunek = 0;
    double licznik = 0;
    double vz1 = 0;
    double vz2 = 0;
    for (const auto & muon : muons) {
      licznik++;
      sumpx += muon.px();
      sumpy += muon.py();
      sumpz += muon.pz();
      sumenergy += muon.energy();
      sumpt += muon.pt();
      ladunek += muon.charge();
      if (licznik==1) vz1 = muon.vz();
      if (licznik == 2) vz2 = muon.vz();
    }
    if(sumpt>10 && ladunek == 0 && abs(vz1-vz2)<0.25){
      masan2m5_15_50->Fill(sqrt(sumenergy*sumenergy-sumpx*sumpx-sumpy*sumpy-sumpz*sumpz));
      masan2m5_15_100->Fill(sqrt(sumenergy*sumenergy-sumpx*sumpx-sumpy*sumpy-sumpz*sumpz));
      masan2m5_15_200->Fill(sqrt(sumenergy*sumenergy-sumpx*sumpx-sumpy*sumpy-sumpz*sumpz));
      masan2m0_12->Fill(sqrt(sumenergy*sumenergy-sumpx*sumpx-sumpy*sumpy-sumpz*sumpz));
      masan2m0_120->Fill(sqrt(sumenergy*sumenergy-sumpx*sumpx-sumpy*sumpy-sumpz*sumpz));
    }
  }
  
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Wierzcholki
  int licznik = 0;
  double vz1 = 0;
  double vz2 = 0;
  int ladunek = 0;
  for (const auto & muon : muons){
    licznik++;
    ladunek += muon.charge();
    //cout<<licznik<<endl;
    if(licznik==1){
      vz1 = muon.vz();
    }
    if(licznik == 2){
      vz2 = muon.vz();
      //histo->Fill(abs(vz1-vz2));
    }
    if(ladunek==0){
      vertex->Fill(abs(vz1-vz2));
    }
  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //delta fi miedzy dwoma dzetami
  
  if (jets.size()==2){
    double phi1 = 0;
    double phi2 = 0;
    vector<double> phi;
    for (const auto & jet : jets){
      phi.push_back(jet.phi());
    }
    phi1  = phi[0];
    phi2 = phi[1];
    double delta_phi = abs(phi1 -  phi2);
    if (delta_phi > M_PI) delta_phi = 2*M_PI - delta_phi;
    jetdphi->Fill(delta_phi);
  }
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //delta fi kontra jet pt
  
  
  
  
  int licznikphi = 0;
  double phi1 = 0;
  double phi2 = 0;
  double jetpt = 0;
  double prawda = 0;
  double prawda2 = 0;
  for (const auto & jet : jets){
    licznikphi++;
    //cout<<licznik<<endl;
    if(licznikphi==1 && jet.pt()>50){
      phi1 = jet.phi();
      jetpt = jet.pt();
      prawda = 1;
      //cout<<"bat"<<endl;
    }
    if(licznikphi == 2 && jet.pt()>50 && prawda == 1){
      phi2 = jet.phi();
      prawda = 2;
    }
    if(prawda==2){
      double delta_phi = abs(phi1 -  phi2);
      if (delta_phi > M_PI) delta_phi = 2*M_PI - delta_phi;
      //cout<<"bat"<<endl;
      deltaphivspt -> Fill(jetpt,delta_phi);
    }

  }
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Liczba jet i muon
  
  long unsigned int licznik1 = 0;
  long unsigned int licznik2 = 0;
  long unsigned int licznik3 = 0;
  long unsigned int licznik4 = 0;
  long unsigned int licznik5 = 0;
  long unsigned int licznik6 = 0;
  for (const auto & jet : jets){
    if(jet.pt()>50){
      licznik1++;
      if(licznik1==jets.size()) jetcount50->Fill(jets.size());
    }
    if(jet.pt()>100){
      licznik2++;
      if(licznik2==jets.size()) jetcount100->Fill(jets.size());
    }
    if(jet.pt()>200){
      licznik3++;
      if(licznik3==jets.size()) jetcount200->Fill(jets.size());
    }


  }

  for (const auto & muon : muons){
    if(muon.pt()>5){
      licznik4++;
      if(licznik4==muons.size()) muoncount5->Fill(muons.size());
    }
    if(muon.pt()>10){
      licznik5++;
      if(licznik5==muons.size()) muoncount10->Fill(muons.size());
    }
    if(muon.pt()>20){
      licznik6++;
      if(licznik6==muons.size()) muoncount20->Fill(muons.size());
    }

  }
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Finalny rysunek
  
  for (const auto & jet : jets) {
    if(jet.pt()<50) continue;
    bool bPassed5 = 0;
    bool bPassed10 = 0;
    bool bPassed20 = 0;
    bool bPassed50 = 0;
    bool bPassed1_5 = 0;
    bool bPassed1_10 = 0;
    bool bPassed1_20 = 0;
    bool bPassed1_50 = 0;
    for (const auto & muon : muons) {
      if(!muon.isGlobalMuon()) continue;
      if(!muon.isMediumMuon()) continue;
      double DeltaR = deltaR(muon, jet);
      if(DeltaR>0.3) continue;
      if(muon.pt()>5) bPassed5 = true;
      if(muon.pt()>5) bPassed50 = true;
      if(muon.pt()>10) bPassed10 = true;
      if(muon.pt()>24) bPassed20 = true;
    }
    pjeta5 -> Fill(bPassed5, jet.eta());
    pjeta10 -> Fill(bPassed10, jet.eta());
    pjeta20 -> Fill(bPassed20, jet.eta());
    pjeta50 -> Fill(bPassed50, jet.eta());

    for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it){
      if(it->hwQual()<12) continue;
      double DeltaR  = deltaR(*it, jet);
      if(DeltaR>0.3) continue;
      if(it->pt()>5) bPassed1_5 = true;
      if(it->pt()>10) bPassed1_10 = true;
      if(it->pt()>22) bPassed1_20 = true;
      if(it->pt()>5) bPassed1_50 = true;
    }
    pjeta1_5 -> Fill(bPassed1_5, jet.eta());
    pjeta1_10 -> Fill(bPassed1_10, jet.eta());
    pjeta1_20 -> Fill(bPassed1_20, jet.eta());
    pjeta1_50 -> Fill(bPassed1_50, jet.eta());

    if(jet.pt()<100) continue;
    bool bPassed100 = 0;
    bool bPassed1_100 = 0;
    for (const auto & muon : muons) {
      if(!muon.isGlobalMuon()) continue;
      if(!muon.isMediumMuon()) continue;
      double DeltaR = deltaR(muon, jet);
      if(DeltaR>0.3) continue;
      if(muon.pt()>5) bPassed100 = true;
    }
    pjeta100 -> Fill(bPassed100, jet.eta());

    for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
      if(it->hwQual()<12) continue;
      double DeltaR  = deltaR(*it, jet);
      if(DeltaR>0.3) continue;
      if(it->pt()>5) bPassed1_100 = true;
    }
    pjeta1_100 -> Fill(bPassed1_100, jet.eta());

    if(jet.pt()<200) continue;
    bool bPassed200 = 0;
    bool bPassed1_200 = 0;
    bool bPassed200_24 = 0;
    bool bPassed1_200_24 = 0;
    for (const auto & muon : muons) {
      if(!muon.isGlobalMuon()) continue;
      if(!muon.isMediumMuon()) continue;
      double DeltaR = deltaR(muon, jet);
      if(DeltaR>0.3) continue;
      if(muon.pt()>5) bPassed200 = true;
      if(muon.pt()>24) bPassed200_24 = true;
    }
    pjeta200 -> Fill(bPassed200, jet.eta());
    pjeta200_24 -> Fill(bPassed200_24, jet.eta());

    for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
      if(it->hwQual()<12) continue;
      double DeltaR  = deltaR(*it, jet);
      if(DeltaR>0.3) continue;
      if(it->pt()>5) bPassed1_200 = true;
      if(it->pt()>22) bPassed1_200_24 = true;
    }
    pjeta1_200 -> Fill(bPassed1_200, jet.eta());
    pjeta1_200_24 -> Fill(bPassed1_200_24, jet.eta());

    /**
      bool bPassed = 0;
      if (jet.pt()>50 && muon.pt()>5 && muon.isGlobalMuon() && muon.isMediumMuon()){
          bPassed = DeltaR < 0.2;
          peta -> Fill(bPassed, jet.eta());
          cout<<"Jet  ponad 5"<<endl;
      } 
      
      if (jet.pt()>50 && muon.pt()>10 && muon.isGlobalMuon() && muon.isMediumMuon()){
          bPassed = DeltaR < 0.2;
          peta3 -> Fill(bPassed, jet.eta());
          cout<<"Jet  ponad 10"<<endl;
      } 
      
      if (jet.pt()>50 && muon.pt()>20 && muon.isGlobalMuon() && muon.isMediumMuon()){
          bPassed = DeltaR < 0.2;
          peta4 -> Fill(bPassed, jet.eta());
          cout<<"Jet  ponad 20"<<endl;
      } 
      
    }
    
    
    

    for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
        double DeltaR = reco::deltaR(jet, *it);
        bool bPassed;
        if (jet.pt()>50 && it->pt()>5 && it->hwQual()>=12){
          bPassed = DeltaR < 0.4;
          peta2 -> Fill(bPassed, jet.eta());
        } 
      
    }
    */
    
  }
  

  for (const auto & muon : muons) {
    if(muon.pt()<5) continue;
    bool bPassedm50 = 0;
    bool bPassedm100 = 0;
    bool bPassedm200 = 0;
    for (const auto & jet : jets) {
      if(!muon.isGlobalMuon()) continue;
      if(!muon.isMediumMuon()) continue;
      double DeltaR = deltaR(muon, jet);
      if(DeltaR>0.3) continue;
      if(jet.pt()>50) bPassedm50 = true;
      if(jet.pt()>100) bPassedm100 = true;
      if(jet.pt()>200) bPassedm200 = true;
    }
    pmeta50 -> Fill(bPassedm50, muon.eta());
    pmeta100 -> Fill(bPassedm100, muon.eta());
    pmeta200 -> Fill(bPassedm200, muon.eta()); 
  }

  for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it){
    if(it->pt()<5) continue;
    bool bPassedm1_50 = 0;
    bool bPassedm1_100 = 0;
    bool bPassedm1_200 = 0;
    for(const auto & jet : jets){
      if(it->hwQual()<12) continue;
      double DeltaR  = deltaR(*it, jet);
      if(DeltaR>0.3) continue;
      if(jet.pt()>50) bPassedm1_50 = true;
      if(jet.pt()>100) bPassedm1_100 = true;
      if(jet.pt()>200) bPassedm1_200 = true;
    }
    pmeta1_50 -> Fill(bPassedm1_50, it->eta());
    pmeta1_100 -> Fill(bPassedm1_100, it->eta());
    pmeta1_200 -> Fill(bPassedm1_200, it->eta());
  }

  for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it){
    if(it->pt()<15) continue; //15GeV
    bool bPassedm1t = 0;
    for (const auto & muon : muons) {
      if(it->hwQual()<12) continue;
      if(!muon.isGlobalMuon()) continue;
      if(!muon.isMediumMuon()) continue;
      double DeltaR = deltaR(muon, *it);
      if(DeltaR>0.3) continue;
      if(muon.pt()>5) bPassedm1t = true;
      
    }
    pl1tm5 -> Fill(bPassedm1t, it->eta());
    
  }

  for (const auto & muon : muons) {
    if(muon.pt()<10) continue; //Wieksze niz 10 GeV
    bool bPassedl1tm = 0;
    for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it){
      if(it->hwQual()<12) continue;
      if(!muon.isGlobalMuon()) continue;
      if(!muon.isMediumMuon()) continue;
      double DeltaR = deltaR(muon, *it);
      if(DeltaR>0.3) continue;
      if(it->pt()>5) bPassedl1tm = true;
    }
    pml1t5 -> Fill(bPassedl1tm, muon.eta());
  }

  for (const auto & muon : muons) {
    if(muon.pt()<50) continue;
    bool bPassedmuon50 = 0;
    for (const auto & jet : jets) {
      if(jet.pt()<100) continue;
      double DeltaR = deltaR(muon,  jet);
      if(DeltaR<0.8) continue;
      for (const auto & muon : muons) {
        if(!muon.isGlobalMuon()) continue;
        if(!muon.isMediumMuon()) continue;
        double DeltaR2 = deltaR(muon, jet);
        if(DeltaR2>0.3) continue;
        if(muon.pt()>50) bPassedmuon50 = true;
      }
      pjetamuon50 -> Fill(bPassedmuon50, jet.eta());
    }
  }

  for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
    if(it->pt()<50) continue;
    bool bPassedl1t50 = 0;
    for (const auto & jet : jets) {
      if(jet.pt()<100) continue;
      double DeltaR = deltaR(*it,  jet);
      if(DeltaR<0.3) continue;
      for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
        if(it->hwQual()<12) continue;
        double DeltaR2 = deltaR(*it, jet);
        if(DeltaR2>0.2) continue;
        if(it->pt()>50) bPassedl1t50 = true;
      }
      pjetal1t50 -> Fill(bPassedl1t50, jet.eta());
    }
  }
  


  ++theEventCount;
  if (theEventCount%wiadomosci==0) {
  //write std io
  //cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" analysed event count:"<<theEventCount << endl;
  }
}

DEFINE_FWK_MODULE(Lic);

