
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliExternalTrackParam.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TSystem.h>
#include <TRandom.h>
#include <TGeoGlobalMagField.h>
#include <fstream>
#include "MVertexFinder.h"
#endif

MVertexFinder* vtf = 0;
AliESDVertex diamond;
double bz = 0.;


AliESDEvent *esdEv=0;
AliESDfriend *esdFr=0;
TTree *esdTree = 0;
const int kNBitsPlp = 6;
const int kMaxPlp = (0x1<<kNBitsPlp)-1;

TClonesArray aliceTracks("AliExternalTrackParam");
TClonesArray aliceVtx("AliESDVertex");

int CreatePileUpEvent(float nplp, int startEv=0, int stepEv=0, ULong_t flags=0);

void PrintTrack(Int_t i);
void PrintTracks();
Bool_t LoadESD(const char* path);
TChain* LoadChainESD(const char* inpData);
Int_t LoadEvent(Int_t iev);
Bool_t ResetEvent();
void FetchESDData(int evID, ULong_t flags=0);
void AttachToVertexer();

void InitVertexer(float zRange = 30.f);


void testVtx()
{
  LoadESD("/home/shahoian/Downloads/alice/sim/2016/LHC16j7a/246392/9999/AliESDs.root");
  InitVertexer();
  int icurr = 0;
  icurr = CreatePileUpEvent(3,icurr, 1, 0x44);
  AttachToVertexer();
  vtf->FindVertices();
  //
}

void InitVertexer(float zRange)
{
  if (vtf) delete vtf;
  vtf = new  MVertexFinder();
  vtf->SetZRange(zRange);
}

void AttachToVertexer()
{
  // load tracks to vertexer
  static Bool_t first = kTRUE; 
  if (first) {
    first = kFALSE;
    //
    vtf->SetIRPos(diamond.GetX(),diamond.GetY(),diamond.GetZ());
    vtf->SetIRSig2(diamond.GetXRes()*diamond.GetXRes(),
		   diamond.GetYRes()*diamond.GetYRes(),
		   diamond.GetZRes()*diamond.GetZRes());
  }
  vtf->Reset();
  //
  int ntr = aliceTracks.GetEntriesFast();
  for (int itr=0;itr<ntr;itr++) {
    AliExternalTrackParam* trc = (AliExternalTrackParam*)aliceTracks[itr];
    double dca[2],dcaCov[3];
    if (!trc->PropagateToDCA(&diamond,bz, 0.3, dca, dcaCov)) continue;
    vtf->AddTrack(trc->GetX(),trc->GetY(),trc->GetZ(),
		  trc->GetSigmaY2(), trc->GetSigmaZ2(), trc->GetSigmaZY(), 
		  trc->GetSnp(), trc->GetTgl(), trc->GetAlpha(), trc->GetUniqueID());
  }
}

//________________________________________________________________________
int CreatePileUpEvent(float nuplp, int startEv, int stepEv, ULong_t flags)
{
  // if stepEv<=0 : choose following events randomly, if >0 : use as step
  if (nuplp>0 && nuplp<1) nuplp = 1.;
  int nEvTot = esdTree->GetEntries();
  int nPlp = 0;
  while ( !(nPlp=(nuplp<0 ? gRandom->Poisson(-nuplp) : int(nuplp))) ) {}
  if (nPlp>kMaxPlp) {
    printf("N events to mix (%d) exceeds max allowed %d, overriding\n",nPlp,kMaxPlp);
    nPlp = kMaxPlp;
  }
  printf("Will merge %d events starting from %d, step = %d\n",nPlp, startEv, stepEv);
  //
  int evsRead[kMaxPlp]={0};
  int ev2read = startEv%nEvTot;
  //
  for (int iev=0;iev<nPlp;iev++) {
    if (iev) { // find next event
      Bool_t repeat = kFALSE; // make sure there are no repeating events
      do {
	int step = stepEv>0 ? stepEv : gRandom->Integer(nEvTot);
	ev2read += step;
	ev2read %= nEvTot;
	for (int j=0;j<iev;j++) {if (evsRead[j]==ev2read) repeat = kTRUE; break;}
      } while(repeat);      
    }
    evsRead[iev] = ev2read;
    LoadEvent(ev2read);
    FetchESDData(iev, flags);
  }
  //
  ResetEvent();
  return ++ev2read; // last event read
}

//__________________________________________________
void FetchESDData(int evID, ULong_t flags)
{
  if (!evID) { // reset containers
    aliceVtx.Clear();
    aliceTracks.Clear("C");
  }
  static AliESDVertex vtxDummy;
  const AliESDVertex* vtx = esdEv->GetPrimaryVertex();
  // use 0 vertex if no SPD or Track vertex found
  if (vtx == esdEv->GetPrimaryVertexTPC() || !vtx->GetStatus()) vtx = &vtxDummy;
  new(aliceVtx[evID]) AliESDVertex(*vtx);
  //
  // fetch tracks
  int ntrPlp = aliceTracks.GetEntriesFast();
  int ntr = esdEv->GetNumberOfTracks();
  for (int itr=0;itr<ntr;itr++) {
    const AliESDtrack* trc = esdEv->GetTrack(itr);
    if (flags && (trc->GetStatus()&flags)!=flags) continue;
    // fiducial cut on DCA
    if (TMath::Abs(trc->GetX())>2 || TMath::Abs(trc->GetY())>2) continue;
    AliExternalTrackParam* tradd = new(aliceTracks[ntrPlp++]) AliExternalTrackParam(*trc);
    UInt_t stamp = evID; // to complete
    tradd->SetUniqueID(stamp);
  }
  //
}


//_______________________________________________________
TChain* LoadChainESD(const char* inpData)
{
  //
  const char* chName="esdTree";
  TChain* chain = new TChain(chName);
  //
  TString inpDtStr = inpData;
  if (inpDtStr.EndsWith(".root")) {
    chain->AddFile(inpData);
  }
  else {
    //
    ifstream inpf(inpData);
    if (!inpf.good()) {
      printf("Failed on input filename %s\n",inpData);
      return 0;
    }
    //
    TString flName;
    flName.ReadLine(inpf);
    while ( !flName.IsNull() ) {
      flName = flName.Strip(TString::kBoth,' ');
      if (flName.BeginsWith("//") || flName.BeginsWith("#")) {flName.ReadLine(inpf); continue;}
      flName = flName.Strip(TString::kBoth,',');
      flName = flName.Strip(TString::kBoth,'"');
      printf("Adding %s\n",flName.Data());
      chain->AddFile(flName.Data());
      flName.ReadLine(inpf);
    }
  }
  //
  int n = chain->GetEntries();
  if (n<1) {
    printf("Obtained chain is empty\n");
    return 0;
  }
  else printf("Opened %s chain with %d entries\n",chName,n);
  return chain;
}

//____________________________________
Bool_t LoadESD(const char* path)
{
  if (!(esdTree = LoadChainESD(path))) return 0;
  esdEv = new AliESDEvent();
  esdEv->ReadFromTree(esdTree);
}

Bool_t ResetEvent()
{
  if (!esdEv) return kFALSE; 
  esdEv->Reset();
  return kTRUE;
}

Int_t LoadEvent(Int_t iev)
{
  if (!ResetEvent()) return -1;
  esdTree->GetEntry(iev);
  esdEv->ConnectTracks();
  if (!TGeoGlobalMagField::Instance()->GetField()) esdEv->InitMagneticField();
  //
  static Bool_t first = kTRUE; 
  if (first) {
    first = kFALSE;
    //
    bz = esdEv->GetMagneticField();
    //
    diamond.SetXv(esdEv->GetDiamondX());
    diamond.SetYv(esdEv->GetDiamondY());
    diamond.SetZv(esdEv->GetDiamondZ());
    //
    double cov[6] = {esdEv->GetSigma2DiamondX(),0.,esdEv->GetSigma2DiamondY(),0.,0.,esdEv->GetSigma2DiamondZ()};
    diamond.SetCovarianceMatrix(cov);
  }
  //  
  return 0;
}


void PrintTracks()
{
  for (int i=0;i<esdEv->GetNumberOfTracks();i++) PrintTrack(i);
}

void PrintTrack(Int_t i)
{
  AliESDtrack* tr = esdEv->GetTrack(i);
  if (!tr) return;
  double p[3];
  tr->GetPxPyPz(p);
  printf("%3d: itsRF%d tpcRF%d trdRF%d tofRF%d | P:%6.2f Px:%+6.2f Py:%+6.2f Pz:%+6.2f\n",
	 i,
	 tr->IsOn(AliESDtrack::kITSrefit),
	 tr->IsOn(AliESDtrack::kTPCrefit),
	 tr->IsOn(AliESDtrack::kTRDrefit),
	 tr->IsOn(AliESDtrack::kTOFrefit),
	 tr->GetP(),p[0],p[1],p[2]);
}

