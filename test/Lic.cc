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
  TH1D *histo1;
  TH1D *histo2;
  TH1D *histo3;
  TH1D *histo4;
  TH1D *histo5;
  TH2D *histo2D;
  //TH2D *histo2D1;
  TEfficiency *peta;
  TEfficiency *peta2;
  TEfficiency *peta3;
  TEfficiency *peta4;
  TEfficiency *peta5;
  TEfficiency *peta6;

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
  histo1 =new TH1D("histo1","test; #GMT; #events",10, 0., 10.);
  histo2 =new TH1D("histo2","test; #GMT; #events",10, 0., 10.);
  histo3 =new TH1D("histo3","test; #GMT; #events",10, 0., 10.);
  histo4 =new TH1D("histo4","test; #GMT; #events",10, 0., 10.);
  histo5 =new TH1D("histo5","test; #GMT; #events",10, 0., 10.);
  //histo1 =new TH1D("histo1","test; #GMT; #events",100, -1., 12.);
  //histo2 =new TH1D("histo2","test; #GMT; #events",100, 0., 20.);
  histo2D = new TH2D("histo2D","y,x,#entries", 100, 0., 1000., 100, 0., 4.);
  //histo2D1 = new TH2D("histo2D1","y,x,#entries", 100, -4., 4, 100, -4., 4.);
  peta = new TEfficiency("peta", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta2 = new TEfficiency("peta2", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta3 = new TEfficiency("peta3", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta4 = new TEfficiency("peta4", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta5 = new TEfficiency("peta5", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
  peta6 = new TEfficiency("peta6", "my efficiency;x;#epcilon", 250, -2.5, 2.5);
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
  histo4->Write();
  histo5->Write();
  histo2D->Write();
  peta->Write();
  peta2->Write();
  peta3->Write();
  peta4->Write();
  peta5->Write();
  peta6->Write();
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
    if (muon.pt()>5){
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
    if (jet.pt()>50){
      //histo->Fill(jets.size());
      //histo1->Fill(muon.energy());
      //histo2->Fill(muon.p());
      //h_2dgaus->Fill(muon.eta(),0.);
      //x.push_back(muon.eta());
    }
    /*
    if(jet.pt()>50){
      histo1->Fill(jet.eta());
    }
    if(jet.pt()>100){
      histo2->Fill(jet.eta());
    }
    if(jet.pt()>200){
      histo3->Fill(jet.eta());
    }
    //histo->Fill(muon.pt());
    //histo1->Fill(muon.vy());
    //histo2->Fill(muon.vz());
  */}

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
  /*
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
      histo->Fill(sqrt(sumenergy*sumenergy-sumpx*sumpx-sumpy*sumpy-sumpz*sumpz));
      histo1->Fill(sqrt(sumenergy*sumenergy-sumpx*sumpx-sumpy*sumpy-sumpz*sumpz));
      histo2->Fill(sqrt(sumenergy*sumenergy-sumpx*sumpx-sumpy*sumpy-sumpz*sumpz));
    }
  }
  */
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /*int licznik = 0;
  double vz1 = 0;
  double vz2 = 0;
  int ladunek = 0;
  for (const auto & muon : muons){
    licznik++;
    ladunek += muon.charge();
    cout<<licznik<<endl;
    if(licznik==1){
      vz1 = muon.vz();
    }
    if(licznik == 2){
      vz2 = muon.vz();
      //histo->Fill(abs(vz1-vz2));
    }
    if(ladunek==0){
      histo3->Fill(abs(vz1-vz2));
    }
  }
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //delta fi miedzy dwoma dzetami
  /*
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
    histo->Fill(delta_phi);
  }*/
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //delta fi kontra jet pt
  /*
  
  
  
  int licznik = 0;
  double phi1 = 0;
  double phi2 = 0;
  double jetpt = 0;
  double prawda = 0;
  double prawda2 = 0;
  for (const auto & jet : jets){
    licznik++;
    cout<<licznik<<endl;
    if(licznik==1 && jet.pt()>50){
      phi1 = jet.phi();
      jetpt = jet.pt();
      prawda = 1;
      cout<<"bat"<<endl;
    }
    if(licznik == 2 && jet.pt()>50 && prawda == 1){
      phi2 = jet.phi();
      prawda = 2;
    }
    if(prawda==2){
      double delta_phi = abs(phi1 -  phi2);
      if (delta_phi > M_PI) delta_phi = 2*M_PI - delta_phi;
      //cout<<"bat"<<endl;
      histo2D -> Fill(jetpt,delta_phi);
    }

  }*/
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Liczba jet i muon
  /*
  long unsigned int licznik1 = 0;
  long unsigned int licznik2 = 0;
  long unsigned int licznik3 = 0;
  long unsigned int licznik4 = 0;
  long unsigned int licznik5 = 0;
  long unsigned int licznik6 = 0;
  for (const auto & jet : jets){
    if(jet.pt()>50){
      licznik1++;
      if(licznik1==jets.size()) histo->Fill(jets.size());
    }
    if(jet.pt()>100){
      licznik2++;
      if(licznik2==jets.size()) histo1->Fill(jets.size());
    }
    if(jet.pt()>200){
      licznik3++;
      if(licznik3==jets.size()) histo2->Fill(jets.size());
    }


  }

  for (const auto & muon : muons){
    if(muon.pt()>5){
      licznik4++;
      if(licznik4==muons.size()) histo3->Fill(muons.size());
    }
    if(muon.pt()>10){
      licznik5++;
      if(licznik5==muons.size()) histo4->Fill(muons.size());
    }
    if(muon.pt()>20){
      licznik6++;
      if(licznik6==muons.size()) histo5->Fill(muons.size());
    }

  }
  */
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //Finalny rysunek
  
  for (const auto & jet : jets) {
    if(jet.pt()<50) continue;
    bool bPassed5 = 0;
    bool bPassed10 = 0;
    bool bPassed20 = 0;
    for (const auto & muon : muons) {
      if(!muon.isGlobalMuon()) continue;
      if(!muon.isMediumMuon()) continue;
      double DeltaR = deltaR(muon, jet);
      if(DeltaR>0.2) continue;
      if(muon.pt()>5) bPassed5 = true;
      if(muon.pt()>10) bPassed10 = true;
      if(muon.pt()>24) bPassed20 = true;
    }
    peta -> Fill(bPassed5, jet.eta());
    peta3 -> Fill(bPassed10, jet.eta());
    peta4 -> Fill(bPassed20, jet.eta());

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
    */
    
    

    for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
        double DeltaR = reco::deltaR(jet, *it);
        bool bPassed;
        if (jet.pt()>50 && it->pt()>5 && it->hwQual()>=12){
          bPassed = DeltaR < 0.4;
          peta2 -> Fill(bPassed, jet.eta());
        } 
      
    }
    
    
  }
  

  for (const auto & muon : muons) {
    
      
    for (l1t::MuonBxCollection::const_iterator it = gmts.begin(bxNumber); it != gmts.end(bxNumber); ++it) {
      double DeltaR = reco::deltaR(muon, *it);
      bool bPassed;
      if (muon.pt()>5 && muon.isGlobalMuon() && muon.isMediumMuon() && it->pt()>5 && it->hwQual()>=12){
          bPassed = DeltaR < 0.4;
          peta5 -> Fill(bPassed, muon.eta());
          peta6 -> Fill(bPassed, it->eta());
      } 
      
    }
    
    
    
  }
  


  ++theEventCount;
  if (theEventCount%wiadomosci==0) {
  //write std io
  //cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" analysed event count:"<<theEventCount << endl;
  }
}

DEFINE_FWK_MODULE(Lic);

