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

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include <sstream>


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
  //TH1D *histo2;
  //TH2D *h_2dgaus;

  edm::EDGetTokenT< vector<pat::Muon> > theMuonToken;
  edm::EDGetTokenT<l1t::MuonBxCollection> theGmtToken;

};


Lic::Lic(const edm::ParameterSet& conf)
  : theConfig(conf), theEventCount(0)
{
  cout <<" CTORXX" << endl;
  theMuonToken = consumes< vector<pat::Muon> >( theConfig.getParameter<edm::InputTag>("muonSrc"));
  theGmtToken  = consumes<l1t::MuonBxCollection>( theConfig.getParameter<edm::InputTag>("gmtSrc"));
}

Lic::~Lic()
{
  cout <<" DTOR" << endl;
}

void Lic::beginJob()
{
  //create a histogram
  histo =new TH1D("histo","test; #GMT; #events",100, -4., 4.);
  histo1 =new TH1D("histo1","test; #GMT; #events",100, -4., 4.);
  //histo2 =new TH1D("histo2","test; #GMT; #events",100, 0., 20.);
  //h_2dgaus = new TH2D("h_2dqaus","y,x,#entries", 100, -4., 4, 100, -4., 4.);
  cout << "HERE Lic::beginJob()" << endl;
}

void Lic::endJob()
{
  //make a new Root file
  TFile myRootFile( theConfig.getParameter<std::string>("outHist").c_str(), "RECREATE");
  //write histogram data
  histo->Write();
  histo1->Write();
  //histo2->Write();
  //h_2dgaus->Write();
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
  if (debug && theEventCount%wiadomosci==0) std::cout <<" number of muons: " << muons.size() <<std::endl;
  for (const auto & muon : muons) {
    if (debug) std::cout <<" reco muon pt: "<<muon.pt()<<std::endl;
    //if (muon.phi()>1 || muon.phi()<-1){
      histo->Fill(muon.phi());
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
    std::cout <<"GMT: "<<it->phi()<<std::endl;
    }
    //if (it->phi()>1 || it->phi()<-1){
        histo1->Fill(it->phi());
      //h_2dgaus->Fill(0.,it->eta()); 
      //y.push_back(it->eta());
    //}
    //histo->Fill(it->hwEta());
  }
  //for(long unsigned int i = 1; i<y.size(); i++){
    //h_2dgaus->Fill(x[i],y[i]);
  //}
  
  ++theEventCount;
  if (theEventCount%wiadomosci==0) {
  //write std io
  cout <<"*** Cwiczenie, analyze event: " << ev.id()<<" analysed event count:"<<theEventCount << endl;
  }
}

DEFINE_FWK_MODULE(Lic);

