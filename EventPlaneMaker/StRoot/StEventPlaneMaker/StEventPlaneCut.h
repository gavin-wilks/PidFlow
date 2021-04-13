#ifndef StEventPlaneCut_h
#define StEventPlaneCut_h

#include "TObject.h"
#include "TString.h"

class StPicoDst;
class StPicoEvent;
class StPicoTrack;

class StEventPlaneCut : public TObject
{
  public:
    StEventPlaneCut(int energy);
    virtual ~StEventPlaneCut();

    bool isMinBias(StPicoEvent*);
    bool isBES();
    bool isPileUpEvent(int, int, int); // refmult/grefmult & nTofMatch & nTofHits
    bool passEventCut(StPicoDst*);

    bool passTrackBasic(StPicoTrack*);
    bool passTrackEp(StPicoTrack*, StPicoEvent*);

    bool passTrackChargedFlowEast(StPicoTrack*, StPicoEvent*);
    bool passTrackChargedFlowWest(StPicoTrack*, StPicoEvent*);
    bool passTrackChargedFlowFull(StPicoTrack*, StPicoEvent*);

    double getBeta(StPicoTrack*, StPicoDst*); // return beta of picoTrack (tof || -999)
    double getPrimaryMass2(StPicoTrack*, StPicoDst*); // return m^2 picoTrack (primary || -999)
    double getGlobalMass2(StPicoTrack*, StPicoDst*); // return m^2 of picoTrack (global || -999)

  private:
    // int mMatchedToF;
    // int mN_prim;
    // int mN_non_prim;
    int mEnergy;

    ClassDef(StEventPlaneCut,1)
};
#endif
