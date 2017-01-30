#include "MVertexFinder.h"
#include "AliSymMatrix.h"
#include <math.h>
#include <algorithm>
//#include <array>

ClassImp(MVertexFinder)

const float  MVertexFinder::kAlmost0F = 1e-12;
const double MVertexFinder::kAlmost0D  = 1e-16;
const float  MVertexFinder::kHugeF = 1.e12; //1.f/MVertexFinder::kAlmost0F;
const float  MVertexFinder::kDefTukey = 5.0f;

using namespace std;

//______________________________________________
MVertexFinder::MVertexFinder() :
  mVtxTracks()
  ,mVertices()
  ,mUseZSorting(true)
  ,mUseConstraint(true)

  ,mRemovedTracks(0)
  
  ,mMaxVtxIter(100)
  ,mMinTracksPerVtx(2)
  ,mMinChangeZ(10e-4f)
  ,mStopScaleChange(0.95)
  ,mSigma2Accept(3.0)
  ,mTukey2I(1./25.f)
  ,mZRange(30.0f)
{
  mIRPos[0]=mIRPos[1]=mIRPos[2]=0.f;
  mIRSig2I[0]=mIRSig2I[1]=mIRSig2I[2]=0.f;
}

//______________________________________________
bool MVertexFinder::FindNextVertex(std::vector<MVertexFinder::vtxTrack> &tracks,float zseed,float sigScale2Ini, float zmin,float zmax)
{
  // find next vertex usin supplied tracks array and extra Z constraints
  static int level = 0;
  level++;
  printf("\nStarting at level %d : Zseed=%+e sg2=%e | %+e<Z<%e\n",level,zseed,sigScale2Ini,zmin,zmax);
  
  vertex vtx;
  vtx.mXYZ[0] = mIRPos[0];
  vtx.mXYZ[1] = mIRPos[1];
  vtx.mXYZ[2] = zseed;
  //
  int nIter = 0;
  int ntr = tracks.size();
  if (ntr<mMinTracksPerVtx) {
    printf("Stopping level %d: ntr=%d\n",level,ntr);
    level--;
    return false;
  }
  //
  float scaleSigma2 = sigScale2Ini;

  bool res = false, finalize = false;
  while(nIter++<mMaxVtxIter) {
    //
    printf(">>#%3d Ntr:%4d Vtx: %+e %+e %+e Sig: %f\n",nIter,ntr, vtx.mXYZ[0],vtx.mXYZ[1],vtx.mXYZ[2], scaleSigma2);
    //
    float scaleSigma2Old = scaleSigma2, oldZ = vtx.mXYZ[2];
    res = FitVertex(tracks,vtx, scaleSigma2, false, zmin, zmax); // fit with current seed and scaling sigma
    if (!res) break;

    float zChange = fabs(vtx.mXYZ[2] - oldZ);
    //
    float sigRat = scaleSigma2/scaleSigma2Old;
    printf("<<#%3d Ntacc:%4d Vtx: %+e %+e %+e Sig:%f -> %f %f (Dz:%+.4f) Stamp: %d\n",
	   nIter,vtx.mNTracks, vtx.mXYZ[0],vtx.mXYZ[1],vtx.mXYZ[2],
	   scaleSigma2Old,scaleSigma2,sigRat,zChange,vtx.mStamp);
    //
    if (scaleSigma2<1.f) {
      printf("sigmaScale went below min., stop iterations\n");
      finalize=true;
      break;
    }
    //
    if ( (sigRat<1.0f && sigRat>mStopScaleChange) || zChange<mMinChangeZ) { // sigma does not drop enough anymore, check convergence
      if (zChange<mMinChangeZ && scaleSigma2<mSigma2Accept) { // converged, finalize the vertex	
	printf("Converged in dZ or dSigmaChange, stop iterations\n");
	finalize=true;
	break;
      }
      else if (vtx.mNTracks>mMinTracksPerVtx) { // stuck between 2 attractors?
	float zLR[2]={0.f};
	res = false;
	LRAttractors(tracks,zmin,zmax,vtx.mXYZ[2],zLR); // side to be discarded has zLR == kHugeF (>> than zmax)
	if (zLR[0]<zmax) res |= FindNextVertex(tracks,zLR[0],scaleSigma2, zmin, vtx.mXYZ[2]); // left: zmin<ZL<curZ
	if (zLR[1]<zmax) res |= FindNextVertex(tracks,zLR[1],scaleSigma2, vtx.mXYZ[2],zmax);  // right: curZ<ZR<zmax	
	break; // stop looping around current seed
      }
      else { // did not converge, disable contributors
	printf("Did not converge, ntrAcc=%d\n",vtx.mNTracks);
	res = false;
	break;
      }
    }
    //
  }
  if (finalize) {
    scaleSigma2 = 1.f;
    if ((res=FitVertex(tracks,vtx, scaleSigma2, true, zmin,zmax)) ) { // final fit with error extraction
      int vID = mVertices.size();
      mVertices.push_back(vtx);
      for (int itr=ntr;itr--;) {	
	vtxTrack &trc = tracks[itr];
	if (trc.CanUse(zmin,zmax) && trc.mWgh>0.f) {
	  trc.mVtxID = vID;
	  mRemovedTracks++;
	}
      }
    }
  }
  else if (!res) { // disable tracks which failed to make a vertex
    int ntrAcc = 0;
    for (int itr=ntr;itr--;) {
      vtxTrack &trc = tracks[itr];
      if (trc.CanUse(zmin,zmax) && trc.mWgh>0.f) {
	trc.mVtxID = vtxTrack::kDiscarded; // flag as discarded
	trc.mWgh = 0.f;
	ntrAcc++; // tmp
	mRemovedTracks++;
      }
    }
    printf("Disabling %d tracks after failure\n",ntrAcc);
  }
  //
  printf("Stopping level %d: res=%d\n",level,res);
  PrintVertices();
  level--;
 
  return res;
}

//______________________________________________
bool MVertexFinder::FindNextVertex(std::vector<int> &trcIDs,float zseed,float sigScale2Ini)
{
  // find next vertex usin supplied track indices and extra Z constraints
  static int level = 0;
  level++;
  printf("\nStarting at level %d : Zseed=%+e sg2=%e\n",level,zseed,sigScale2Ini);
  
  vertex vtx;
  vtx.mXYZ[0] = mIRPos[0];
  vtx.mXYZ[1] = mIRPos[1];
  vtx.mXYZ[2] = zseed;
  //
  int nIter = 0;
  int ntr = trcIDs.size();
  if (ntr<mMinTracksPerVtx) {
    printf("Stopping level %d: ntr=%d\n",level,ntr);
    level--;
    return false;
  }
  //
  float scaleSigma2 = sigScale2Ini;

  bool res = false, finalize = false;
  while(nIter++<mMaxVtxIter) {
    //
    printf(">>#%3d Ntr:%4d Vtx: %+e %+e %+e Sig: %f\n",nIter,ntr, vtx.mXYZ[0],vtx.mXYZ[1],vtx.mXYZ[2], scaleSigma2);
    //
    float scaleSigma2Old = scaleSigma2, oldZ = vtx.mXYZ[2];
    res = FitVertex(trcIDs,vtx, scaleSigma2, false); // fit with current seed and scaling sigma
    if (!res) break;

    float zChange = fabs(vtx.mXYZ[2] - oldZ);
    //
    float sigRat = scaleSigma2/scaleSigma2Old;
    printf("<<#%3d Ntacc:%4d Vtx: %+e %+e %+e Sig:%f -> %f %f (Dz:%+.4f) | Stamp: %d\n",nIter,vtx.mNTracks, vtx.mXYZ[0],vtx.mXYZ[1],vtx.mXYZ[2],
	   scaleSigma2Old,scaleSigma2,sigRat,zChange,vtx.mStamp);
    //
    if (scaleSigma2<1.f) {
      printf("sigmaScale went below min., stop iterations\n");
      finalize=true;
      break;
    }
    //
    // sigma does not drop enough anymore, check convergence
    if ( (sigRat<1.0f && sigRat>mStopScaleChange) ||
	 (zChange<mMinChangeZ && scaleSigma2>10*mSigma2Accept && sigRat>mStopScaleChange)) { // recheck 10*...
      if (zChange<mMinChangeZ && scaleSigma2<mSigma2Accept) { // converged, finalize the vertex	
	printf("Converged in dZ or dSigmaChange, stop iterations\n");
	finalize=true;
	break;
      }
      else if (vtx.mNTracks>mMinTracksPerVtx) { // stuck between 2 attractors?
	float zLR[2]={0.f};
	vector<int> trcIdsLR[2];
	//
	res = false;
	if (LRAttractors(trcIDs,vtx.mXYZ[2],trcIdsLR,zLR)) { // split to 2 subsample at left/right of current vtx
	  // side to be discarded has zLR == 2*kHugeF
	  //	for (int ilr=2;ilr--;) {
	  for (int ilr=0;ilr<2;ilr++) {
	    if (zLR[ilr]>=kHugeF || int(trcIdsLR[ilr].size())<mMinTracksPerVtx) DisableTracks(trcIdsLR[ilr]);
	    else res |= FindNextVertex(trcIdsLR[ilr],zLR[ilr],scaleSigma2);
	  }
	}
	else DisableTracks(trcIDs); // splitting failed - abandon tracks in current sample
	//
	break; // stop looping around current seed
      }
      else { // did not converge, disable contributors
	printf("Did not converge, ntrAcc=%d\n",vtx.mNTracks);
	res = false;
	break;
      }
    }
    //
  }
  if (finalize) {
    scaleSigma2 = 1.f;
    if ((res=FitVertex(trcIDs,vtx, scaleSigma2, true)) ) { // final fit with error extraction
      int vID = mVertices.size();
      mVertices.push_back(vtx);
      for (int itr=ntr;itr--;) {	
	vtxTrack &trc = mVtxTracks[trcIDs[itr]];              // don't use CanUse(zmn,zmx) version: z-range guaranteed 
	if (trc.CanUse() && trc.mWgh>0.f) {
	  trc.mVtxID = vID;   // when vector if track indices is used
	  mRemovedTracks++;
	}
      }
    }
  }
  else if (!res) { // disable tracks which failed to make a vertex
    DisableTracks(trcIDs);
  }
  //
  printf("Stopping level %d: res=%d\n",level,res);
  PrintVertices();
  level--;
 
  return res;
}

//______________________________________________
void MVertexFinder::LRAttractors(const std::vector<vtxTrack> &tracks,
				 float zmn,float zmx,float currZ,
				 float zLR[2]) const
{
  // find Z of attractors on the left and right from current Z
  float wghLR[2] = {0.f};
  zLR[0] = zLR[1] = 0.f;
  for (int itr=tracks.size();itr--;) {
    const vtxTrack &trc = tracks[itr];
    if (trc.CanUse(zmn,zmx) && trc.mWgh>0.f) { // track is free and has non-0 weight wrt vtx seed
      int ind = (trc.mZ<currZ) ? 0 : 1; // weight to the left and right of current vertex
      wghLR[ind] += trc.mWgh;
      zLR[ind] += trc.mWgh*trc.mZ;
    }    
  }
  for (int is=2;is--;) zLR[is] = wghLR[is]>kAlmost0F ? zLR[is]/wghLR[is] : 2.f*kHugeF;
  printf("Split to Zl=%+e  Zr=%+e\n",zLR[0],zLR[1]);
}

//______________________________________________
bool MVertexFinder::LRAttractors(const std::vector<int> &src, float currZ,
				 std::vector<int> tgt[2],float zLR[2]) const
{
  // find Z of attractors on the left and right from current Z
  float wghLR[2] = {0.f};
  zLR[0] = zLR[1] = 0.f;
  int ntr = src.size();
  tgt[0].reserve(ntr);
  tgt[1].reserve(ntr);  
  //
  for (int itr=ntr;itr--;) {
    int idx = src[itr];
    const vtxTrack &trc = mVtxTracks[idx];
    if (trc.mWgh>0.f) { // track has non-0 weight wrt vtx seed
      int side = (trc.mZ<currZ) ? 0 : 1; // to the left and right of current vertex ?
      tgt[side].push_back(idx);
      wghLR[side] += trc.mWgh;
      zLR[side] += trc.mWgh*trc.mZ;      
    }
  }
  bool splitOK = true;
  for (int is=2;is--;) {
    if (wghLR[is]>kAlmost0F) zLR[is] = zLR[is]/wghLR[is];
    else {
      zLR[is] = 2.f*kHugeF;
      splitOK = false;
    }
  }
  printf("Split to Zl=%+e (Ntr=%ld)  Zr=%+e (Ntr=%ld)\n",zLR[0],tgt[0].size(), zLR[1],tgt[1].size());
  return splitOK;
}

//______________________________________________
bool MVertexFinder::FindVertices()
{
  int ntr = mVtxTracks.size();
  if (ntr<mMinTracksPerVtx) return false;
  mRemovedTracks = 0;
  if (mUseZSorting) sort(mVtxTracks.begin(),mVtxTracks.end());
  
  const float zresChar = 200e-4; // characteristic Z resolution
  float sig2ini = (1.+2.*mZRange)/zresChar, zini=0.0f;
  sig2ini *= sig2ini;
  //
  //  while (FindNextVertex( mVtxTracks, zini, sig2ini,-mZRange,mZRange)) {}
  //
  vector<int> trcIDs;
  trcIDs.reserve(ntr);    
  do {
    trcIDs.clear();
    for (int itr=0;itr<ntr;itr++) if (mVtxTracks[itr].CanUse()) trcIDs.push_back(itr);
    PrintTracks();
  } while (FindNextVertex( trcIDs, zini, sig2ini) || (ntr-mRemovedTracks)>=mMinTracksPerVtx);
  
    
}

//______________________________________________
bool MVertexFinder::FitVertex(std::vector<MVertexFinder::vtxTrack> &tracks, MVertexFinder::vertex &vtx,
			      float &scaleSigma2, bool fillError, float zmin,float zmax)
{
  int ntAcc=0,ntr = tracks.size();
  if (ntr<2) return false;
  float* curVtx = vtx.mXYZ;
  //
  double wghSum=0,wghChi2=0; 
  double cxx=0,cxy=0,cxz=0,cx0=0,cyy=0,cyz=0,cy0=0,czz=0,cz0=0;
  float scaleSig2ITuk2I = mTukey2I/scaleSigma2;
  //
  for (int itr=ntr;itr--;) {
    vtxTrack &trc = tracks[itr];
    if (!trc.CanUse(zmin,zmax) || !CheckVertexTrackStamps(vtx,trc)) continue; // the track is invalidated or out of range, skip
    // determine weight of the track
    // current vertex in the track proper frame
    float vlocX =  curVtx[0]*trc.mCosAlp+curVtx[1]*trc.mSinAlp;
    float vlocY = -curVtx[0]*trc.mSinAlp+curVtx[1]*trc.mCosAlp;
    float vlocZ =  curVtx[2];
    // DCA from track to vertex. The track is straight-line defined by tangents
    float dy    = trc.mY + trc.mTgP*(vlocX-trc.mX) - vlocY;
    float dz    = trc.mZ + trc.mTgL*(vlocX-trc.mX) - vlocZ;
    // weighted distance to vertex
    float syyI(trc.mSig2YI),szzI(trc.mSig2ZI),syzI(trc.mSigYZI);
    float chi2T = 0.5f*(dy*dy*syyI + dz*dz*szzI) + dy*dz*syzI;
    float wghT = (1.f-chi2T*scaleSig2ITuk2I);
    if (wghT<kAlmost0F)  {
      trc.mWgh = 0.f;
      continue;
    }
    wghSum  += wghT;
    wghChi2 += wghT*chi2T;
    //
    double wdz = dz*wghT;
    syyI *= wghT;
    syzI *= wghT;
    szzI *= wghT;
    trc.mWgh = wghT;
    //
    // aux variables
    double tmpSP = trc.mSinAlp*trc.mTgP, tmpCP = trc.mCosAlp*trc.mTgP
      ,tmpSC = trc.mSinAlp+tmpCP, tmpCS = -trc.mCosAlp+tmpSP
      ,tmpCL = trc.mCosAlp*trc.mTgL, tmpSL = trc.mSinAlp*trc.mTgL
      ,tmpYXP = trc.mY-trc.mTgP*trc.mX, tmpZXL = trc.mZ-trc.mTgL*trc.mX
      ,tmpCLzz = tmpCL*szzI, tmpSLzz = tmpSL*szzI, tmpSCyz = tmpSC*syzI, tmpCSyz = tmpCS*syzI
      ,tmpCSyy = tmpCS*syyI, tmpSCyy = tmpSC*syyI, tmpSLyz = tmpSL*syzI, tmpCLyz = tmpCL*syzI;
    //
    // symmetric matrix equation 
    cxx += tmpCL*(tmpCLzz+tmpSCyz+tmpSCyz)+tmpSC*tmpSCyy;          // dchi^2/dx/dx
    cxy += tmpCL*(tmpSLzz+tmpCSyz)+tmpSL*tmpSCyz+tmpSC*tmpCSyy;    // dchi^2/dx/dy
    cxz += -trc.mSinAlp*syzI-tmpCLzz-tmpCP*syzI;                   // dchi^2/dx/dz
    cx0 += -(tmpCLyz+tmpSCyy)*tmpYXP-(tmpCLzz+tmpSCyz)*tmpZXL;     // RHS 
    //
    cyy += tmpSL*(tmpSLzz+tmpCSyz+tmpCSyz)+tmpCS*tmpCSyy;          // dchi^2/dy/dy
    cyz += -(tmpCSyz+tmpSLzz);                                     // dchi^2/dy/dz
    cy0 += -tmpYXP*(tmpCSyy+tmpSLyz)-tmpZXL*(tmpCSyz+tmpSLzz);     // RHS
    //
    czz += szzI;                                                    // dchi^2/dz/dz
    cz0 += tmpZXL*szzI+tmpYXP*syzI;                                 // RHS
    //
    ntAcc++;
  }
  //
  vtx.mNTracks = ntAcc;
  if (ntAcc<mMinTracksPerVtx) return false;
  //
  if (mUseConstraint) {
    // impose meanVertex constraint, i.e. account terms (V_i-Constrain_i)^2/sig2constr_i for i=X,Y 
    // in the fit chi2 definition
    cxx += mIRSig2I[0];
    cx0 += mIRSig2I[0]*mIRPos[0];
    cyy += mIRSig2I[1];
    cy0 += mIRSig2I[1]*mIRPos[1];
    //    czz += mIRSig2I[2];             // don't constraint Z
    //    cz0 += mIRSig2I[2]*mIRPos[2];
  }
  AliSymMatrix mat(3);
  double vec[3] = {cx0,cy0,cz0};
  mat(0,0) = cxx;
  mat(0,1) = cxy;
  mat(0,2) = cxz;
  mat(1,1) = cyy;
  mat(1,2) = cyz;
  mat(2,2) = czz;
  // 
  if (!mat.SolveChol(vec,true)) return false;
  for (int i=3;i--;) curVtx[i] = float(vec[i]);
  if (fillError) {
    vtx.mCov[0] = mat(0,0);
    vtx.mCov[1] = mat(0,1);
    vtx.mCov[2] = mat(1,1);
    vtx.mCov[3] = mat(0,2);
    vtx.mCov[4] = mat(1,2);
    vtx.mCov[5] = mat(2,2);
    vtx.mNTracks = ntAcc;
    vtx.mChi2 = 2.f*(ntAcc-wghSum)/scaleSig2ITuk2I;         // calculate chi^2
  }
  scaleSigma2 = wghChi2/wghSum;
  return true;
}

//______________________________________________
bool MVertexFinder::FitVertex(std::vector<int> &trcIDs, MVertexFinder::vertex &vtx,float &scaleSigma2, bool fillError)
{
  int ntAcc=0,ntr = trcIDs.size();
  if (ntr<2) return false;
  float* curVtx = vtx.mXYZ;
  //
  double wghSum=0,wghChi2=0; 
  double cxx=0,cxy=0,cxz=0,cx0=0,cyy=0,cyz=0,cy0=0,czz=0,cz0=0;
  float scaleSig2ITuk2I = mTukey2I/scaleSigma2;
  //
  for (int itr=ntr;itr--;) {
    vtxTrack &trc = mVtxTracks[trcIDs[itr]];
    if (!trc.CanUse() || !CheckVertexTrackStamps(vtx,trc)) continue; // the track is invalidated or out of range, skip
    // determine weight of the track
    // current vertex in the track proper frame
    float vlocX =  curVtx[0]*trc.mCosAlp+curVtx[1]*trc.mSinAlp;
    float vlocY = -curVtx[0]*trc.mSinAlp+curVtx[1]*trc.mCosAlp;
    float vlocZ =  curVtx[2];
    // DCA from track to vertex. The track is straight-line defined by tangents
    float dy    = trc.mY + trc.mTgP*(vlocX-trc.mX) - vlocY;
    float dz    = trc.mZ + trc.mTgL*(vlocX-trc.mX) - vlocZ;
    // weighted distance to vertex
    float syyI(trc.mSig2YI),szzI(trc.mSig2ZI),syzI(trc.mSigYZI);
    float chi2T = 0.5f*(dy*dy*syyI + dz*dz*szzI) + dy*dz*syzI;
    float wghT = (1.f-chi2T*scaleSig2ITuk2I);
    if (wghT<kAlmost0F)  {
      trc.mWgh = 0.f;
      continue;
    }
    wghSum  += wghT;
    wghChi2 += wghT*chi2T;
    //
    double wdz = dz*wghT;
    syyI *= wghT;
    syzI *= wghT;
    szzI *= wghT;
    trc.mWgh = wghT;
    //
    // aux variables
    double tmpSP = trc.mSinAlp*trc.mTgP, tmpCP = trc.mCosAlp*trc.mTgP
      ,tmpSC = trc.mSinAlp+tmpCP, tmpCS = -trc.mCosAlp+tmpSP
      ,tmpCL = trc.mCosAlp*trc.mTgL, tmpSL = trc.mSinAlp*trc.mTgL
      ,tmpYXP = trc.mY-trc.mTgP*trc.mX, tmpZXL = trc.mZ-trc.mTgL*trc.mX
      ,tmpCLzz = tmpCL*szzI, tmpSLzz = tmpSL*szzI, tmpSCyz = tmpSC*syzI, tmpCSyz = tmpCS*syzI
      ,tmpCSyy = tmpCS*syyI, tmpSCyy = tmpSC*syyI, tmpSLyz = tmpSL*syzI, tmpCLyz = tmpCL*syzI;
    //
    // symmetric matrix equation 
    cxx += tmpCL*(tmpCLzz+tmpSCyz+tmpSCyz)+tmpSC*tmpSCyy;          // dchi^2/dx/dx
    cxy += tmpCL*(tmpSLzz+tmpCSyz)+tmpSL*tmpSCyz+tmpSC*tmpCSyy;    // dchi^2/dx/dy
    cxz += -trc.mSinAlp*syzI-tmpCLzz-tmpCP*syzI;                   // dchi^2/dx/dz
    cx0 += -(tmpCLyz+tmpSCyy)*tmpYXP-(tmpCLzz+tmpSCyz)*tmpZXL;     // RHS 
    //
    cyy += tmpSL*(tmpSLzz+tmpCSyz+tmpCSyz)+tmpCS*tmpCSyy;          // dchi^2/dy/dy
    cyz += -(tmpCSyz+tmpSLzz);                                     // dchi^2/dy/dz
    cy0 += -tmpYXP*(tmpCSyy+tmpSLyz)-tmpZXL*(tmpCSyz+tmpSLzz);     // RHS
    //
    czz += szzI;                                                    // dchi^2/dz/dz
    cz0 += tmpZXL*szzI+tmpYXP*syzI;                                 // RHS
    //
    ntAcc++;
  }
  //
  vtx.mNTracks = ntAcc;
  if (ntAcc<mMinTracksPerVtx) return false;
  //
  if (mUseConstraint) {
    // impose meanVertex constraint, i.e. account terms (V_i-Constrain_i)^2/sig2constr_i for i=X,Y 
    // in the fit chi2 definition
    cxx += mIRSig2I[0];
    cx0 += mIRSig2I[0]*mIRPos[0];
    cyy += mIRSig2I[1];
    cy0 += mIRSig2I[1]*mIRPos[1];
    //    czz += mIRSig2I[2];             // don't constraint Z
    //    cz0 += mIRSig2I[2]*mIRPos[2];
  }
  AliSymMatrix mat(3);
  double vec[3] = {cx0,cy0,cz0};
  mat(0,0) = cxx;
  mat(0,1) = cxy;
  mat(0,2) = cxz;
  mat(1,1) = cyy;
  mat(1,2) = cyz;
  mat(2,2) = czz;
  // 
  if (!mat.SolveChol(vec,true)) return false;
  for (int i=3;i--;) curVtx[i] = float(vec[i]);
  if (fillError) {
    vtx.mCov[0] = mat(0,0);
    vtx.mCov[1] = mat(0,1);
    vtx.mCov[2] = mat(1,1);
    vtx.mCov[3] = mat(0,2);
    vtx.mCov[4] = mat(1,2);
    vtx.mCov[5] = mat(2,2);
    vtx.mNTracks = ntAcc;
    vtx.mChi2 = 2.f*(ntAcc-wghSum)/scaleSig2ITuk2I;        // calculate chi^2
  }
  scaleSigma2 = wghChi2/wghSum;
  return true;
}

//______________________________________________
void MVertexFinder::Reset()
{
  mVtxTracks.clear();
  mVertices.clear();
}

//______________________________________________
void MVertexFinder::AddTrack(float x,float y,float z,float sy2,float sz2, float syz, float snp, float tgl, float alp,
			     unsigned int stamp)
{
  // add new track to the pool
  float cs = (1.f-snp)*(1.f+snp);
  if (cs<kAlmost0F) return;
  double detI = sy2*sz2 - syz*syz;
  if (fabs(detI)<kAlmost0D) return;
  detI = 1./detI;
  //
  vtxTrack trc;
  //
  trc.mStamp = stamp;
  trc.mTgP = snp/sqrtf(cs);
  trc.mTgL = tgl;
  trc.mX = x;
  trc.mY = y;
  trc.mZ = z;
  trc.mSig2YI = sz2*detI;
  trc.mSig2ZI = sy2*detI;
  trc.mSigYZI =-syz*detI;
  trc.mCosAlp = cosf(alp);
  trc.mSinAlp = sinf(alp);
  trc.mWgh = 0;
  trc.mVtxID = vtxTrack::kNoVtx;
  mVtxTracks.push_back(trc);
  return;
}

//______________________________________________
void MVertexFinder::PrintVertices() const
{
  ///< print current vertices
  int nv = mVertices.size();
  for (int iv=0;iv<nv;iv++) {
    const vertex& vtx = mVertices[iv];
    printf("#%2d %+e %+e %+e | Ntracks = %4d | Chi2/ntr = %7.2f | Stamp: %d\n",
	   iv,vtx.mXYZ[0],vtx.mXYZ[1],vtx.mXYZ[2],vtx.mNTracks,vtx.mChi2/vtx.mNTracks,vtx.mStamp);
    int id = 0;
    for (int ie=0;ie<3;ie++) {printf("    "); for (int je=ie;je<3;je++) printf("%+e ",vtx.mCov[id++]); printf("\n");}
  }
}

//______________________________________________
void MVertexFinder::PrintTracks() const
{
  ///< print tracks
  int nt = mVtxTracks.size();
  printf("Ntracks = %5d Nremoved = %5d\n",nt,mRemovedTracks);
  for (int it=0;it<nt;it++) {
    const vtxTrack &trc = mVtxTracks[it];
    printf("#%4d Z:%+e W:%+e tgL:%+5.2f | Vtx: = %3d | Stamp: %d\n",it,trc.mZ,trc.mWgh,trc.mTgL,trc.mVtxID,trc.mStamp);
  }
}

//______________________________________________
float MVertexFinder::GetTukey() const
{
  // convert 1/tukey^2 to tukey
  return sqrtf(1./mTukey2I);
}

//______________________________________________
void MVertexFinder::SelectFreeTracks(const std::vector<int> &src,
				     std::vector<int> &tgt,
				     float zmn,float zmx) const
{
  // select free tracks in z-range from src and store their pointers in tgt
  int ntr = src.size();
  tgt.clear();
  tgt.reserve(ntr);
  //
  if (mUseZSorting) {
    for (int it=0;it<ntr;++it) {
      int tID = src[it];
      const vtxTrack& trc = mVtxTracks[tID];
      if (trc.mZ<zmn) continue;  // use binary search or lower_bound ?
      if (trc.mZ>zmx) break;
      if (trc.CanUse()) tgt.push_back(tID);
    }
  }
  else {
    for (int it=0;it<ntr;++it) {
      int tID = src[it];
      if (mVtxTracks[tID].CanUse(zmn,zmx)) tgt.push_back(tID);
    }
  }
  
}

//______________________________________________
bool MVertexFinder::CheckVertexTrackStamps(vertex& vtx, const vtxTrack& trc)
{
  // check if vertex and track track-stamps are compatible. Undefined time-stamp is
  // compatible with both defined and undefined ones, but only defined track stamp
  // may enforce the vertex time stamp
  bool res = true;
  //  return res;
  //
  if (IsStampDefined(trc.mStamp)) {
    if (IsStampDefined(vtx.mStamp)) res = trc.mStamp==vtx.mStamp; // at the moment require compatibility as strict equality
    else vtx.mStamp = trc.mStamp; //  track will set its stamp to the vertex
  }
  // undefined track-stamp is compatible with any vertex stamp...
  return  res;
  //
}
