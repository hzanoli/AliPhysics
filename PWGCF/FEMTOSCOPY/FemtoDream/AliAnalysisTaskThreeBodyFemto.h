/*
 * AliAnalysisTaskNanoXioton.h
 *
 *  Created on: May 13, 2019
 *      Author: schmollweger
 */

#ifndef PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTO_H_
#define PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTO_H_
#include "AliAnalysisTaskSE.h"
#include "AliFemtoDreamEventCuts.h"
#include "AliFemtoDreamEvent.h"
#include "AliFemtoDreamTrackCuts.h"
#include "AliFemtoDreamTrack.h"
#include "AliFemtoDreamv0.h"
#include "AliFemtoDreamv0Cuts.h"
#include "AliFemtoDreamCollConfig.h"
#include "AliFemtoDreamPairCleaner.h"
#include "AliFemtoDreamPartCollection.h"
#include "AliFemtoDreamControlSample.h"

class AliAnalysisTaskThreeBodyFemto : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskThreeBodyFemto();
  AliAnalysisTaskThreeBodyFemto(const char* name, bool isMC);
  virtual ~AliAnalysisTaskThreeBodyFemto();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  void ResetGlobalTrackReference();
  void StoreGlobalTrackReference(AliVTrack *track);
  
  void FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes);
  void FillPairDistributionPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, TH1F* sameEventDistributionPL);
  void SetMixedEvent(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer);
  void FillTripletDistributionMEPP(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes);
  // test different mixing
  void SetMixedEventOnlyPLambdaTEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>*fPartContainer);
  void FillTripletDistributionMEPPTEST(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes);
  void FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes);

  void SetRunTaskLightWeight(bool light) {
    fisLightWeight = light;
  }
  void SetEventCuts(AliFemtoDreamEventCuts* evtCuts) {
    fEventCuts = evtCuts;
  }
  void SetProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fProton = trkCuts;
  }
  void SetAntiProtonCuts(AliFemtoDreamTrackCuts* trkCuts) {
    fAntiProton = trkCuts;
  }
  void Setv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fLambda = v0Cuts;
  }
  void SetAntiv0Cuts(AliFemtoDreamv0Cuts* v0Cuts) {
    fAntiLambda = v0Cuts;
  }
  void SetCorrelationConfig(AliFemtoDreamCollConfig* config) {
    fConfig=config;
  }  
  void SetRunThreeBodyHistograms(bool RunThreeBodyHistos) {
    fRunThreeBody=RunThreeBodyHistos;
  }
  static TLorentzVector RelativePairMomentum(TLorentzVector &PartOne, TLorentzVector &PartTwo);
 private:
  AliAnalysisTaskThreeBodyFemto(const AliAnalysisTaskThreeBodyFemto &task);
  AliAnalysisTaskThreeBodyFemto &operator=(const AliAnalysisTaskThreeBodyFemto &task);
  bool fisLightWeight;//
  AliFemtoDreamEvent* fEvent;//!
  AliFemtoDreamEventCuts* fEventCuts;//
  TList* fEvtList;//!
  AliFemtoDreamTrack* fTrack;//!
  AliFemtoDreamTrackCuts* fProton;//
  TList* fProtonList;//!
  TList* fProtonMCList;//!
  AliFemtoDreamTrackCuts* fAntiProton;//
  TList* fAntiProtonList;//!
  TList* fAntiProtonMCList;//!
  AliFemtoDreamv0* fv0;//!
  AliFemtoDreamv0Cuts* fLambda;//
  TList* fLambdaList;
  TList* fLambdaMCList;
  AliFemtoDreamv0Cuts* fAntiLambda;//
  TList* fAntiLambdaList;
  TList* fAntiLambdaMCList;
  AliFemtoDreamCollConfig *fConfig; //
  AliFemtoDreamPairCleaner *fPairCleaner;   //!
  AliFemtoDreamPartCollection *fPartColl;   //!
  TList *fResults;//!
  // Three particles same event
  TList *fResultsThreeBody;//!
  bool fRunThreeBody;
  TH1F* sameEventDistributionPL;
  TH1F* sameEventDistributionPPL;
  TH1F* sameEventDistributionAPAPAL;
  TH1F* sameEventDistributionPPP;
  TH1F* sameEventDistributionAPAPAP;
  TH1F* sameEventDistributionPLL;
  TH1F* sameEventDistributionAPALAL;
  TH1F* sameEventDistributionLLL;
  TH1F* sameEventDistributionALALAL;
  // Three particles mixed events
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainer;
  std::vector<std::vector<std::vector<AliFemtoDreamPartContainer>>> fPartContainerTEST;
  TH1F* mixedEventDistributionPL;
  TH1F* mixedEventDistributionPPL;
  TH1F* mixedEventDistributionAPAPAL;
  TH1F* mixedEventDistributionPPP;
  TH1F* mixedEventDistributionAPAPAP;
  TH1F* mixedEventDistributionPLL;
  TH1F* mixedEventDistributionAPALAL;
  TH1F* mixedEventDistributionLLL;
  TH1F* mixedEventDistributionALALAL;
  TH1F* mixedEventDistributionPPLTEST;
  TH1F* mixedEventDistributionAPAPALTEST;

  TList *fResultsQA;//!
  AliFemtoDreamControlSample *fSample;   //!
  TList *fResultsSample;//!
  TList *fResultsSampleQA;//!
  int fTrackBufferSize;//
  AliVTrack **fGTI;  //!
  ClassDef(AliAnalysisTaskThreeBodyFemto,1)
};

#endif /* PWGCF_FEMTOSCOPY_FEMTODREAM_ALIANALYSISTASKTHREEBODYFEMTO_H_ */


