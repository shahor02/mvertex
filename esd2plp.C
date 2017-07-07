
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliHeader.h"
#include "AliGenEventHeader.h"
#include "AliExternalTrackParam.h"
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TArrayF.h>
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
AliHeader* header=0;
AliGenEventHeader* mcHeader = 0;
TTree *esdTree = 0;
TTree *headTree = 0;
const int kNBitsPlp = 6;
const int kMaxPlp = (0x1<<kNBitsPlp)-1;
int evsRead[kMaxPlp]={0};

TClonesArray aliceTracks("AliExternalTrackParam");
TClonesArray aliceVtx("AliESDVertex");
TClonesArray genVtx("AliESDVertex");

void LoadInput(const char* inpData,Bool_t mc=kTRUE);
TChain* LoadChain(const char* inpData, const char* chName="esdTree",const char* altFile=0);
int CreatePileUpEvent(float nplp, int startEv=0, int stepEv=0, ULong_t flags=0, ULong_t flagsRej=0);

Int_t LoadEvent(Int_t iev);
Bool_t ResetEvent();
void FetchESDData(int evID, int evIDOrig, ULong_t flags=0, ULong_t flagsRej=0);
void AttachToVertexer();

void InitVertexer(float zRange = 30.f);

void esd2plp(int ntest, float nevload=-2,const char* inpData="esds.txt", Bool_t mc=kTRUE)
{
  LoadInput(inpData,mc);
  InitVertexer();
  int icurr = 0;

  for (int iev=0;iev<ntest;iev++) {
    icurr = CreatePileUpEvent(nevload,icurr, -1, AliESDtrack::kITSrefit, AliESDtrack::kITSpureSA);
    AttachToVertexer();
    int nv = vtf->FindVertices();
    //    vtf->PrintVertices();
    //
    // print vertices
    printf("Vertices found: %d\n",nv);
    for (int iv=0;iv<nv;iv++) {
      const MVertexFinder::vertex* vtx = vtf->GetVertex(iv);
      int st = vtf->GetCheckStamps() ? vtx->mStamp-1 : (vtx->mStamp)/100-1;
      
      printf("#%2d XYZ: %+e %+e %+e | Ntracks = %4d | Chi2/ntr = %7.2f | Stamp: %d\n",
	     iv,vtx->mXYZ[0],vtx->mXYZ[1],vtx->mXYZ[2],vtx->mNTracks,vtx->mChi2/vtx->mNTracks,st);
    }
    //
    if (mc) {
      printf("GenVertices: %d\n",genVtx.GetEntries());
      for (int ivg=0;ivg<genVtx.GetEntries();ivg++) {
	AliESDVertex* vtg = (AliESDVertex*)genVtx[ivg];
	AliESDVertex* vtr = (AliESDVertex*)aliceVtx[ivg];
	printf("#%2d Ev#%4d XYZ: %+.3f %+.3f %+8.3f : Ngen %3d | OrigRec: XYZ: %+.3f %+.3f %+8.3f : Nrec %3d\n",
	       ivg,evsRead[ivg],
	       vtg->GetX(),vtg->GetY(),vtg->GetZ(), vtg->GetNContributors(),
	       vtr->GetX(),vtr->GetY(),vtr->GetZ(), vtr->GetNIndices());
      }
    }
    vtf->PrintTracks();
    printf("\n");
  }
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
int CreatePileUpEvent(float nuplp, int startEv, int stepEv, ULong_t flags, ULong_t flagsRej)
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
  int ev2read = startEv%nEvTot;
  //
  for (int iev=0;iev<nPlp;iev++) {
    evsRead[iev] = 0;
    if (iev) { // find next event
      Bool_t repeat = kFALSE; // make sure there are no repeating events
      do {
	repeat = kFALSE;
	int step = stepEv>0 ? stepEv : gRandom->Integer(nEvTot);
	ev2read += step;
	ev2read %= nEvTot;
	for (int j=0;j<iev;j++) {if (evsRead[j]==ev2read) repeat = kTRUE; break;}
      } while(repeat);      
    }
    evsRead[iev] = ev2read;
    LoadEvent(ev2read);
    FetchESDData(iev, ev2read, flags, flagsRej);
  }
  //
  ResetEvent();
  return ++ev2read; // last event read
}

//__________________________________________________
void FetchESDData(int evID, int evIDOrig, ULong_t flags, ULong_t flagsRej)
{
  if (!evID) { // reset containers
    aliceVtx.Clear();
    genVtx.Clear();
    aliceTracks.Clear("C");
  }
  static AliESDVertex vtxDummy;
  const AliESDVertex* vtx = esdEv->GetPrimaryVertex();
  // use 0 vertex if no SPD or Track vertex found
  if (vtx == esdEv->GetPrimaryVertexTPC() || !vtx->GetStatus()) vtx = &vtxDummy;
  new(aliceVtx[evID]) AliESDVertex(*vtx);
  if (headTree) { // we stored gen vertex in the TPCvertex
    new(genVtx[evID]) AliESDVertex(*esdEv->GetPrimaryVertexTPC());
  }
  //
  // fetch tracks
  int ntrPlp = aliceTracks.GetEntriesFast();
  int ntr = esdEv->GetNumberOfTracks();
  for (int itr=0;itr<ntr;itr++) {
    const AliESDtrack* trc = esdEv->GetTrack(itr);
    if (flags && (trc->GetStatus()&flags)!=flags) continue;
    if (flagsRej && (trc->GetStatus()&flagsRej)!=0) continue;
    // fiducial cut on DCA
    if (TMath::Abs(trc->GetX())>2 || TMath::Abs(trc->GetY())>2) continue;
    AliExternalTrackParam* tradd = new(aliceTracks[ntrPlp++]) AliExternalTrackParam(*trc);
    UInt_t stamp = evIDOrig+1; // to complete
    tradd->SetUniqueID(stamp);
  }
  //
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
  // MC vertex
  if (headTree) {
    TArrayF mcv(3);
    mcHeader = header->GenEventHeader();
    mcHeader->PrimaryVertex(mcv);
    // we will store the genVertex in the TPC vertex
    AliESDVertex* vtTPC = (AliESDVertex*)esdEv->GetPrimaryVertexTPC();
    vtTPC->SetXv(mcv[0]);
    vtTPC->SetYv(mcv[1]);
    vtTPC->SetZv(mcv[2]);
    vtTPC->SetNContributors(mcHeader->NProduced());
  }
  //
  return 0;
}


//__________________________________________________________________
TChain* LoadChain(const char* inpData, const char* chName, const char* altFile)
{
  TChain* chain = new TChain(chName);
  //
  TString altFileS = altFile; 
  TString inpDtStr = inpData;
  //
  if (inpDtStr.EndsWith(".root")) {
    if (!altFileS.IsNull()) { // alternative file should be used
      inpDtStr.ReplaceAll(gSystem->BaseName(inpDtStr),altFileS);
    }
    chain->AddFile(inpDtStr.Data());
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
      if (!altFileS.IsNull()) { // alternative file should be used
	flName.ReplaceAll(gSystem->BaseName(flName),altFileS);
      }
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

//__________________________________________
void LoadInput(const char* inpData,Bool_t mc)
{
  esdTree = LoadChain(inpData,"esdTree");  
  esdEv = new AliESDEvent();
  esdEv->ReadFromTree(esdTree);
  if (mc) {
    headTree = LoadChain(inpData,"TE","galice.root");
    header = new AliHeader();
    headTree->SetBranchAddress("Header",&header);
    esdTree->AddFriend(headTree,"mchead");
  }
}

