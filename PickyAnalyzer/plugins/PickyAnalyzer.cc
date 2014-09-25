// -*- C++ -*-
//
// Package:    PikyOpt/PickyAnalyzer
// Class:      PickyAnalyzer
// 
/**\class PickyAnalyzer PickyAnalyzer.cc PikyOpt/PickyAnalyzer/plugins/PickyAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Raffaella Radogna
//         Created:  Tue, 23 Sep 2014 09:32:30 GMT
//
//
// system include files
#include <memory>
#include <fstream>
#include <sys/time.h>
#include <string>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <math.h>

// root include files
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "SimDataFormats/Track/interface/SimTrackContainer.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"

//
// class declaration
//
class TFile;
class TH1F;
class TH2F;

using namespace std;
using namespace edm;
using namespace reco;

class PickyAnalyzer : public edm::EDAnalyzer {
   public:
      explicit PickyAnalyzer(const edm::ParameterSet&);
      ~PickyAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
      std::map<std::string,TH1F*> histContainer_;
      std::map<std::string,TH2F*> histContainer2D_;
      std::string theDataType;
      edm::InputTag muonLabel_;
      int numberOfSimTracks;
      int numberOfMuons;
      int numberOfRecTracks;
      int numberOfMatchedTracks;
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
PickyAnalyzer::PickyAnalyzer(const edm::ParameterSet& iConfig):
histContainer_(),
histContainer2D_()
{
	//now do what ever initialization is needed
  	theDataType = iConfig.getUntrackedParameter<string>("DataType");
	muonLabel_ = iConfig.getUntrackedParameter<edm::InputTag>("MuonCollectionLabel");
    	if(theDataType != "RealData" && theDataType != "SimData")
    		cout<<"Error in Data Type!!"<<endl;
    
  	numberOfSimTracks=0;
	numberOfMuons=0;
  	numberOfRecTracks=0;
	numberOfMatchedTracks=0;
}


PickyAnalyzer::~PickyAnalyzer()
{
   	// do anything here that needs to be done at desctruction time
   	// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

int countDrMatching(SimTrackContainer::const_iterator simTrack, Handle<reco::MuonCollection> muons){

	int countMatchingTmp = 0;

	double simEtaTmp = (*simTrack).momentum().eta();
	double simPhiTmp = (*simTrack).momentum().phi();

  	reco::MuonCollection::const_iterator muon;
	for (muon = muons->begin(); muon != muons->end(); ++muon){
		double dR = sqrt(pow((simEtaTmp-muon->eta()),2) + pow((simPhiTmp-muon->phi()),2));

	    	if(!(muon->pt())) continue;
		if(dR > 0.1) continue;
		countMatchingTmp++;
	}

	return countMatchingTmp;
}

struct MyMuon{

	double pt, eta, phi;
	int charge;
	double chi2, ndof;
	double etaSTA, phiSTA;
	int chargeSTA;
	bool isGlobal, isTracker, isStandAlone;	
	int matchedTracks, recHits;
	double globalPt, standAlonePt, trackerPt;
	double glbCharge, trkCharge;
	double dytPt, pickyPt, tpfmsPt;
	double dytCharge, pickyCharge, tpfmsCharge;
	double bestTunePCharge, bestTunePPt, bestTunePChi2, bestTunePNdof;
	double bestCharge, bestPt, bestChi2, bestNdof;
	int muonType; 

};

MyMuon muonMatching(const Event& iEvent, SimTrackContainer::const_iterator simTrack, edm::InputTag muonLabel_){

	int numMatch = 0;

	double simEta = (*simTrack).momentum().eta();
	double simPhi = (*simTrack).momentum().phi();

  	MyMuon tmpMuon;
	tmpMuon.pt = -999;
	tmpMuon.eta = -999;
	tmpMuon.phi = -999;
	tmpMuon.charge = -999;
	tmpMuon.chi2 = -999; 
	tmpMuon.ndof = -999; 
	tmpMuon.etaSTA = -999;
	tmpMuon.phiSTA = -999;
	tmpMuon.chargeSTA = -999;
	tmpMuon.recHits = -999;
	tmpMuon.matchedTracks = -999;
	tmpMuon.isGlobal = -999;
	tmpMuon.isTracker = -999;
	tmpMuon.isStandAlone = -999;
	tmpMuon.globalPt = -999;
	tmpMuon.standAlonePt = -999;
	tmpMuon.trackerPt = -999;
	tmpMuon.dytPt = -999;
	tmpMuon.pickyPt = -999;
	tmpMuon.tpfmsPt = -999; 
	tmpMuon.dytCharge = -999;
	tmpMuon.pickyCharge = -999;
	tmpMuon.tpfmsCharge = -999; 
	tmpMuon.trkCharge = -999;
	tmpMuon.glbCharge = -999;
	tmpMuon.bestCharge = -999;
	tmpMuon.bestPt = -999;
	tmpMuon.bestChi2 = -999;
	tmpMuon.bestNdof = -999;
	tmpMuon.bestTunePCharge = -999;
	tmpMuon.bestTunePPt = -999;
	tmpMuon.bestTunePChi2 = -999;
	tmpMuon.bestTunePNdof = -999;
	tmpMuon.muonType = -999;

  	Handle<reco::MuonCollection> muons;
  	iEvent.getByLabel(muonLabel_, muons);
  	reco::MuonCollection::const_iterator muon;
	for (muon = muons->begin(); muon != muons->end(); ++muon){

		if(!(muon->pt())) continue;
		//if(!(abs(muon->eta()) > minEta_ && abs(muon->eta()) < maxEta_)) continue;
		double dR = sqrt(pow((simEta-muon->eta()),2) + pow((simPhi-muon->phi()),2));
		if(dR > 0.1) continue;

		tmpMuon.pt = muon->pt();
		tmpMuon.eta = muon->eta();
		tmpMuon.phi = muon->phi();
		tmpMuon.charge = muon->charge();
		tmpMuon.isGlobal = muon->isGlobalMuon();
		tmpMuon.isTracker = muon->isTrackerMuon();
		tmpMuon.isStandAlone = muon->isStandAloneMuon();
		reco::Muon::MuonTrackType tmpType = muon->tunePMuonBestTrackType();
		tmpMuon.muonType = tmpType;

		numMatch++;
		
		TrackRef muonRef = muon->combinedMuon(); 	//global muon
		TrackRef trackerRef = muon->innerTrack(); 	// tracker muon: tracker only hits
		TrackRef standAloneRef = muon->outerTrack();	// sta muon: all muon hits
		TrackRef pickyRef = muon->pickyTrack();		//picky algo refit, use only selected muon hits
		TrackRef dytRef = muon->dytTrack();		//redo pattern recognition with dynamic truncation
		TrackRef tpfmsRef = muon->tpfmsTrack();		//tracker + include only first muon hit(s)
		TrackRef bestTunePRef = muon->tunePMuonBestTrack();//tune p algo per scegliere il best
		TrackRef bestRef = muon->muonBestTrack();	

		if(bestTunePRef.isNonnull()){
 			tmpMuon.bestTunePPt = bestTunePRef->pt();
 			tmpMuon.bestTunePCharge = bestTunePRef->charge();
 			tmpMuon.bestTunePChi2 = bestTunePRef->chi2();
 			tmpMuon.bestTunePNdof = bestTunePRef->ndof();
		}
		if(bestRef.isNonnull()){
 			tmpMuon.bestPt = bestRef->pt();
 			tmpMuon.bestCharge = bestRef->charge();
 			tmpMuon.bestChi2 = bestRef->chi2();
 			tmpMuon.bestNdof = bestRef->ndof();
		}
		if(trackerRef.isNonnull()){
 			tmpMuon.trackerPt = trackerRef->pt();
 			tmpMuon.trkCharge = trackerRef->charge();
		}
		if(standAloneRef.isNonnull()){
			tmpMuon.standAlonePt = standAloneRef->pt();
			tmpMuon.etaSTA = standAloneRef->eta();
			tmpMuon.phiSTA = standAloneRef->phi();
			tmpMuon.chargeSTA = standAloneRef->charge();
		}
		if(dytRef.isNonnull()){
			tmpMuon.dytPt = dytRef->pt();
			tmpMuon.dytCharge = dytRef->charge();
		}
		if(pickyRef.isNonnull()){
			tmpMuon.pickyPt = pickyRef->pt();
			tmpMuon.pickyCharge = pickyRef->charge();
		}
		if(tpfmsRef.isNonnull()){
			tmpMuon.tpfmsPt = tpfmsRef->pt();
			tmpMuon.tpfmsCharge = tpfmsRef->charge();
		}
		if(muonRef.isNonnull()){
			tmpMuon.globalPt = muonRef->pt();
			tmpMuon.glbCharge = muonRef->charge();
			tmpMuon.chi2 = muonRef->chi2();
			tmpMuon.ndof = muonRef->ndof();
		}	
	}

	tmpMuon.matchedTracks = numMatch;

	return tmpMuon;

}

// ------------ method called for each event  ------------
void
PickyAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   	using namespace edm;

	//LogTrace("Analyzer") << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
	MuonPatternRecoDumper debug;
	
	// Get the RecTrack collection from the event
	Handle<reco::TrackCollection> staTracks;
	iEvent.getByLabel("globalMuons", staTracks);
	
	Handle<reco::MuonCollection> muons;
  	iEvent.getByLabel(muonLabel_, muons);

	ESHandle<MagneticField> theMGField;
	iSetup.get<IdealMagneticFieldRecord>().get(theMGField);

	ESHandle<GlobalTrackingGeometry> theTrackingGeometry;
	iSetup.get<GlobalTrackingGeometryRecord>().get(theTrackingGeometry);

	double simPt = 0.;
  	double simEta = 0.;
  	double simPhi = 0.;
	
	
	if(theDataType == "SimData"){
	// Get the SimTrack collection from the event	
 	Handle<SimTrackContainer> simTracks;
 	iEvent.getByLabel("g4SimHits",simTracks);
    
 	numberOfRecTracks += staTracks->size();
	numberOfMuons += muons->size();
 	SimTrackContainer::const_iterator simTrack;
	
 	for (simTrack = simTracks->begin(); simTrack != simTracks->end(); ++simTrack){
      		if (!(abs((*simTrack).type()) == 13)) continue; 
  		//if ((*simTrack).noVertex()) continue;
  		//if ((*simTrack).noGenpart()) continue;
		numberOfSimTracks++;
		int drMatching = countDrMatching(simTrack, muons);
		std::cout<<"Matching with: "<<drMatching<<" muons"<<std::endl;
		if(!(drMatching == 1)) continue;
		
			numberOfMatchedTracks++;
			simPt=(*simTrack).momentum().pt();			
			simEta = (*simTrack).momentum().eta();
			simPhi = (*simTrack).momentum().phi();
			int qGen = simTrack->charge();
			//cout<<"Sim pT: "<<(*simTrack).momentum().pt()<<endl;
			//cout<<"Sim Eta: "<<(*simTrack).momentum().eta()<<endl;
			//cout<<"Sim Phi: "<<(*simTrack).momentum().phi()<<endl;
			histContainer_["hSimEta"]->Fill(simEta);
			histContainer_["hSimPhi"]->Fill(simPhi);
			histContainer_["hPtSim"]->Fill(simPt);	
			
			MyMuon muon;
			muon = muonMatching(iEvent, simTrack, muonLabel_);
			
			double qOverP = (muon.charge/muon.pt - qGen/simPt)/(qGen/simPt);
			double qOverPTuneP = (muon.bestTunePCharge/muon.bestTunePPt - qGen/simPt)/(qGen/simPt);
			double qOverPBest = (muon.bestCharge/muon.bestPt - qGen/simPt)/(qGen/simPt);
			
			histContainer_["hRecEta"]->Fill(muon.eta);
			histContainer_["hRecPhi"]->Fill(muon.phi);
			histContainer_["hPtRec"]->Fill(muon.pt);
			
			histContainer_["hPtDiff"]->Fill(muon.pt - simPt);
			histContainer_["hInvPtRes"]->Fill((1/muon.pt - 1/simPt)/(1/simPt));
  			histContainer_["hQoverPtRes"]->Fill((muon.charge/muon.pt - qGen/simPt)/(qGen/simPt));
  			histContainer_["hPtRes"]->Fill((muon.pt - simPt)/simPt);

                        histContainer_["hPtDiff_picky"]->Fill(muon.pickyPt - simPt);
                        histContainer_["hInvPtRes_picky"]->Fill((1/muon.pickyPt - 1/simPt)/(1/simPt));
                        histContainer_["hQoverPtRes_picky"]->Fill((muon.pickyCharge/muon.pickyPt - qGen/simPt)/(qGen/simPt));
                        histContainer_["hPtRes_picky"]->Fill((muon.pickyPt - simPt)/simPt);

			std::cout<<"-----------------------------------------------------------------------------------------"<<std::endl;
			std::cout<<iEvent.id()<<std::endl;
			std::cout<<"q/p RecoMuon: "<<qOverP<<" q/p TuneP: "<<qOverPTuneP<<" q/p Best: "<<qOverPBest<<std::endl;

			std::cout
			<<"SimTrackEta: "<<simPt<<" SimEta: "<<simEta
			<<" RecoMuonPt: "<<setprecision(10)<<muon.pt
			<<" GLBTrackPt: "<<setprecision(10)<<muon.globalPt
			<<" TrkMuonPt: "<<setprecision(10)<<muon.trackerPt
			<<" PickyPt: "<<setprecision(10)<<muon.pickyPt
			<<" DYTPt: "<<setprecision(10)<<muon.dytPt
			<<" TPFMSPt: "<<setprecision(10)<<muon.tpfmsPt
			<<std::endl;

			std::cout<<"Muon type: "<<muon.muonType
			<<" RecoMuonCharge: "<<setprecision(10)<<muon.charge
			<<" GLBTrackCharge: "<<setprecision(10)<<muon.glbCharge
			<<" TrkMuonCharge: "<<setprecision(10)<<muon.trkCharge
			<<" PickyCharge: "<<setprecision(10)<<muon.pickyCharge
			<<" DYTCharge: "<<setprecision(10)<<muon.dytCharge
			<<" TPFMSCharge: "<<setprecision(10)<<muon.tpfmsCharge
			<<std::endl;

			std::cout<<"-----------------------------------------------------------------------------------------"<<std::endl;
		    
 	 } //for simTracks
	} // if SimData

	/*reco::TrackCollection::const_iterator staTrack;
	for (staTrack = staTracks->begin(); staTrack != staTracks->end(); ++staTrack){
    		reco::TransientTrack track(*staTrack,&*theMGField,theTrackingGeometry);
    		recPt = track.impactPointTSCP().momentum().perp();
    		histContainer_["hPtRec"]->Fill(recPt);
    		cout<<"pt: "<<recPt<<" eta: "<<(*staTrack).momentum().eta()<<endl;
    	}*/

} //close analyze


// ------------ method called once each job just before starting event loop  ------------
void 
PickyAnalyzer::beginJob()
{
	// register to the TFileService
	edm::Service<TFileService> fs;
	TH1::SetDefaultSumw2();
	
	histContainer_["hPtSim"] = fs->make<TH1F>("pTSim","p_{T}^{gen} ",461,-2.5,2302.5);
	histContainer_["hSimEta"] = fs->make<TH1F>("PSimEta","SimTrack #eta",100,-2.5,2.5);
	histContainer_["hSimPhi"] = fs->make<TH1F>("PSimPhi","SimTrack #phi",100,-TMath::Pi(),TMath::Pi());
	
	histContainer_["hPtRec"] = fs->make<TH1F>("pTRec","p_{T}^{rec}",461,-2.5,2302.5);
	histContainer_["hRecEta"] = fs->make<TH1F>("PRecEta","Rec muon #eta",100,-2.5,2.5);
	histContainer_["hRecPhi"] = fs->make<TH1F>("PRecPhi","Rec muon #phi",100,-TMath::Pi(),TMath::Pi());

	histContainer_["hPtDiff"] = fs->make<TH1F>("pTDiff","p_{T}^{rec} - p_{T}^{gen} ",250,-120,120);	
	histContainer_["hInvPtRes"] = fs->make<TH1F>("InvPtRes","1/p_{T} Resolution",1000,-5,5);
  	histContainer_["hQoverPtRes"] = fs->make<TH1F>("hQoverPtRes","q/p_{T} Resolution",1000,-5,5);
  	histContainer_["hPtRes"] = fs->make<TH1F>("PtRes","p_{T} Resolution ",1000,-5,5);
	
	histContainer_["hPtDiff_picky"] = fs->make<TH1F>("pTDiff","p_{T}^{rec} - p_{T}^{gen} ",250,-120,120);	
        histContainer_["hInvPtRes_picky"] = fs->make<TH1F>("InvPtRes","1/p_{T} Resolution",1000,-5,5);
        histContainer_["hQoverPtRes_picky"] = fs->make<TH1F>("hQoverPtRes","q/p_{T} Resolution",1000,-5,5);
        histContainer_["hPtRes_picky"] = fs->make<TH1F>("PtRes","p_{T} Resolution ",1000,-5,5);

	//histContainer_["hDeltaPtRec"] = fs->make<TH1F>("DeltapTRec","#Delta p_{T}^{rec}",400,-2,2);


  

}

// ------------ method called once each job just after ending the event loop  ------------
void 
PickyAnalyzer::endJob() 
{
	if(theDataType == "SimData"){
    		cout << endl << endl << "Number of matched tracks: " << numberOfMatchedTracks << endl << "Number of Sim tracks: " << numberOfSimTracks << endl;
  	}
  	cout << "Number of muons: " << numberOfMuons << endl << "Number of Reco tracks (global): " << numberOfRecTracks << endl << endl;
}

// ------------ method called when starting to processes a run  ------------
/*
void 
PickyAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
PickyAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
PickyAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
PickyAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
PickyAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  	//The following says we do not know what parameters are allowed so do no validation
  	// Please change this to state exactly what you do use, even if it is no parameters
  	edm::ParameterSetDescription desc;
  	desc.setUnknown();
  	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(PickyAnalyzer);
