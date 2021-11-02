#include <TSystem>
#include "TStopwatch.h"

class StMaker;
class StChain;
class StPicoDstMaker;
class StPicoDst;
class StPicoEvent;

StChain *chain;

void recoEventPlane(const Char_t *inputFile="../FileList/19p6GeV_2019/pico_prod_random_test.list", const string jobId = "1", const Int_t Mode = 0, const Int_t EpdMode = 0, const Int_t energy = 1)
{
  // mBeamEnergy[NumBeamEnergy] = {"200GeV","54GeV","27GeV"};
  // Mode: 
  // 0 => Gain Correction for ZDC-SMD;
  // 1 => ReCnter Correction for ZDC-SMD East/West & TPC East/West/Full;
  // 2 => Shift Correction for ZDC-SMD East/West & TPC East/West/Full;
  // 3 => Chift Correction for Full ZDC-SMD EP;
  // 4 => Resolution Calculation for ZDC-SMD Sub & TPC Sub/Ran;
  // 5 => Charged Particle Flow w.r.t ZDC-SMD Full & TPC Sub;

  TStopwatch *stopWatch = new TStopwatch();
  stopWatch->Start();

  string SL_version = "pro";
  if(energy == 0||energy==1) SL_version = "SL21c";
  string env_SL = getenv ("STAR");
  if (env_SL.find(SL_version)==string::npos) 
  {
    cout<<"Environment Star Library does not match the requested library in runPicoMixedEventMaker.C. Exiting..."<<endl;
    exit(1);
  }

  Int_t nEvents = 1000000000;
  // Int_t nEvents = 10000;

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoEvent");
  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StRefMultCorr");
  gSystem->Load("StEventPlaneMaker");
  gSystem->Load("StEpdUtil");

  chain = new StChain();
  StPicoDstMaker *picoMaker = new StPicoDstMaker(2,inputFile,"picoDst");

  StEventPlaneMaker *EventPlaneMaker = new StEventPlaneMaker("EventPlane",picoMaker,jobId,Mode,EpdMode,energy,0.3,2.0);

  if( chain->Init()==kStErr ){ 
    cout<<"chain->Init();"<<endl;
    return;
  }

  int total = picoMaker->chain()->GetEntries();
  cout << " Total entries = " << total << endl;
  if(nEvents>total) nEvents = total;

  for (Int_t i=0; i<nEvents; i++)
  {
    if(i != 0 && i%50 == 0)
    {
      cout << "." << flush;
    }
    if(i != 0 && i%500 == 0)
    {
      Float_t event_percent = 100.0*i/nEvents;
      cout << " " << i << "(" << event_percent << "%)" << "\n" << " ==> Processing data " << flush;
    }

    chain->Clear();
    int iret = chain->Make(i);

    if (iret)
    { 
      cout << "Bad return code!" << iret << endl;
      break;
    }

    total++;
  }

  cout << "." << flush;
  cout << " " << nEvents << "(" << 100 << "%)";
  cout << endl;
  cout << "****************************************** " << endl;
  cout << "Work done... now its time to close up shop!"<< endl;
  cout << "****************************************** " << endl;
  chain->Finish();
  cout << "****************************************** " << endl;
  cout << "total number of events  " << nEvents << endl;
  cout << "****************************************** " << endl;

  stopWatch->Stop();
  stopWatch->Print();

  delete EventPlaneMaker;
  delete picoMaker;
  delete chain;
}
