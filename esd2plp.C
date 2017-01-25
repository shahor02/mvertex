
#if !defined(__CINT__) || defined(__MAKECINT__)
#include "AliESDEvent.h"
#include "AliESD.h"
#include "AliESDfriend.h"
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#endif

TFile* flIn=0;
AliESDEvent *esdEv=0;
AliESDfriend *esdFr=0;
TTree *esdTree = 0;

TClonesArray aliceTracks("AliESDtrack");
TClonesArray aliceVtx("AliESDVertex");

void CreatePileUpEvent(float nplp, int startEv=0, int stepEv=0, ULong_t flags=0)

void PrintTrack(Int_t i);
void PrintTracks();
void ConnectFriends();
Int_t LoadESD(const char* path,Bool_t friends=kFALSE);
Int_t LoadEvent(Int_t iev);
Bool_t ResetEvent();
void FetchESDData(int evID, ULong_t flags=0);

void CreatePileUpEvent(float nuplp, int startEv, int stepEv, ULong_t flags)
{
  // if stepEv<=0 : choose following events randomly, if >0 : use as step
  if (nuplp>0 && nuplp<1) nuplp = 1.;
  int nEvTot = esdTree->GetEntries();
  int nPlp = 0;
  while (!nPlp) nuplp<0 ? gRandom->Poisson(-nuplp) : int(nuplp);
  printf("Will merge %d events starting from %d, randomize = %d\n",nPlp, startEv, randomize);
  //
  int evsRead[nPlp];
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
    LoadEvent(ev2read);
    FetchESDData(iev, flags);
  }
  //
  ResetEvent();
  return ev2read; // last event read
}

//__________________________________________________
void FetchESDData(int evID, ULong_t flags)
{
  if (!evID) { // reset containers
    aliceVtx.Clear();
    aliceTracks.Clear("C");
  }
  static AliESDVertex vtxDummy;
  AliESDvertex* vtx = esdEv->GetPrimaryVertex();
  // use 0 vertex if no SPD or Track vertex found
  if (vtx == esdEv->GetPrimaryVertexTPC() || !vtx->GetStatus()) vtx = vtxDummy();
  new(aliceVtx[evID]) AliESDvertex(*vtx);
}

  
//__________________________________________________
Int_t LoadESD(const char* path,Bool_t friends)
{
  flIn = TFile::Open(path);
  if (!flIn) {
    printf("Failed to open %s\n",path);
    return -1;
  }
  //
  esdTree = (TTree*) flIn->Get("esdTree");
  if (!esdTree) {
    printf("No ESDtree found in %s\n",path);
    return -1;
  }
  //
  printf("Loaded esdTree with %d entries\n",(int)esdTree->GetEntries());
  esdEv = new AliESDEvent();
  esdEv->ReadFromTree(esdTree);
  if (friends) ConnectFriends();
  //
  return 0;
}

Bool_t ResetEvent()
{
  if (!esdEv) return kFALSE; 
  esdEv->Reset();
  if (esdFr) esdFr->Reset();
  return kTRUE;
}

Int_t LoadEvent(Int_t iev)
{
  if (!ResetEvent()) return -1;
  esdTree->GetEntry(iev);
  if (esdFr) esdEv->SetESDfriend(esdFr);
  esdEv->ConnectTracks();
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


void ConnectFriends()
{
  // Connect the friends tree as soon as available.
  //
  // Handle the friends first
  //
  if (!esdTree->FindBranch("ESDfriend.")) {
    // Try to add ESDfriend. branch as friend
      TString esdFriendTreeFName;
      esdFriendTreeFName = (esdTree->GetCurrentFile())->GetName();    
      TString basename = gSystem->BaseName(esdFriendTreeFName);
      Int_t index = basename.Index("#")+1;
      basename.Remove(index);
      basename += "AliESDfriends.root";
      TString dirname = gSystem->DirName(esdFriendTreeFName);
      dirname += "/";
      esdFriendTreeFName = dirname + basename;
      //
      TTree* cTree = esdTree->GetTree();
      if (!cTree) cTree = esdTree;      
      cTree->AddFriend("esdFriendTree", esdFriendTreeFName.Data());
      cTree->SetBranchStatus("ESDfriend.", 1);
      esdFr = (AliESDfriend*)(esdEv->FindListObject("AliESDfriend"));
      if (esdFr) cTree->SetBranchAddress("ESDfriend.", &esdFr);
  }
}
