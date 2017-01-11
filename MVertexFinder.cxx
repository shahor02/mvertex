#include "MVertexFinder.h"
#include "AliSymMatrix.h"
#include <math.h>
//#include <array>

ClassImp(MVertexFinder)

const float kAlmost0F = 1e-12;
const float kAlmost0  = 1e-16;

//______________________________________________
MVertexFinder::MVertexFinder() :
  mVtxTracks()
  ,mVertices()
  ,mMaxVtxIter(100)
  ,mMinTracksPerVtx(2)
  ,mScalSigma2Start(1e5)
  ,mMinChangeZ(10e-4f)
  ,mStopScaleChange(0.95)
  ,mSigma2Accept(3.0)
  ,mSigma2Push(10.)
  ,mTukey2(6)
  ,mZRange(30.0f)
{
  mVtxConstraint[0]=mVtxConstraint[1]=mVtxConstraint[2]=0.f;
  mVtxConstrErr[0]=mVtxConstrErr[1]=mVtxConstrErr[2]=0.f;
}

//______________________________________________
bool MVertexFinder::FindNextVertex(float zseed)
{
  //
  static int level = 0;
  level++;
  printf("\nStarting at level %d\n",level);
  //
  vertex vtx;
  vtx.mXYZ[0] = mVtxConstraint[0];
  vtx.mXYZ[1] = mVtxConstraint[1];
  vtx.mXYZ[2] = zseed;
  //
  int nIter = 0;
  int ntr = mVtxTracks.size();
  if (ntr<mMinTracksPerVtx) {
    printf("Sopping level %d: ntr=%d\n",level,ntr);
    level--;
    return false;
  }
  //
  float scaleSigma2 = mScalSigma2Start;
  bool vtxOK = false;

  bool res = false;
  while(nIter++<mMaxVtxIter) {

    double cxx=0,cxy=0,cxz=0,cx0=0,cyy=0,cyz=0,cy0=0,czz=0,cz0=0;
    //
    printf(">>#%3d Ntr:%4d Vtx: %+e %+e %+e Sig: %f\n",nIter,ntr, vtx.mXYZ[0],vtx.mXYZ[1],vtx.mXYZ[2], scaleSigma2);
    //
    float scaleSigma2Old = scaleSigma2, oldZ = vtx.mXYZ[2];
    res = FitVertex(mVtxTracks,vtx, scaleSigma2, false); // fit with current seed and scaling sigma
    if (!res) break;

    float zChange = vtx.mXYZ[2] - oldZ;
    //
    float sigRat = scaleSigma2/scaleSigma2Old;
    printf("<<#%3d Ntacc:%4d Vtx: %+e %+e %+e Sig:%f -> %f %f (Dz:%+.4f)\n",nIter,vtx.mNTracks, vtx.mXYZ[0],vtx.mXYZ[1],vtx.mXYZ[2],
	   scaleSigma2Old,scaleSigma2,sigRat,zChange);
    //
    if (sigRat>mStopScaleChange) { // sigma does not drop enough anymore, check convergence
      if ((fabs(zChange)<mMinChangeZ && scaleSigma2<mSigma2Push) || scaleSigma2<mSigma2Accept) { // converged, finalize the vertex
	break;
      }
      else if (scaleSigma2<mSigma2Push) { // decrease sigme to get rid of outliers
	scaleSigma2 = (mSigma2Accept+scaleSigma2)*0.5;
	if (scaleSigma2<1.f) scaleSigma2 = 1.f; 
	printf("-->pushing sigma to %f\n",scaleSigma2);
	continue;
      }
      else {
	scaleSigma2 *= 0.5;
	if (scaleSigma2<1.f) scaleSigma2 = 1.f; 
	printf("split sigma to %f\n",scaleSigma2);
	continue;
      }
    }
    //
  }
  scaleSigma2 = 1.f;
  if (res && (res=FitVertex(mVtxTracks,vtx, scaleSigma2, true)) ) { // final fit with error extraction
    int vID = mVertices.size();
    mVertices.push_back(vtx);
    for (int itr=ntr;itr--;) {
      vtxTrack &trc = mVtxTracks[itr];
      if (trc.IsUsed() || trc.mWgh==0.f) continue; // the track is invalidated, skip
      trc.SetUsed();
      trc.mVtxID = vID;
      trc.mWgh = 0.f;
    }
  }
  //
  printf("Stopping level %d: res=%d\n",level,res);
  level--;
  return res;
}

//______________________________________________
bool MVertexFinder::FindNextVertex(std::vector<MVertexFinder::vtxTrack> &tracks,float zseed,float sigScale2Ini, float zmin,float zmax)
{
  //
  static int level = 0;
  level++;
  printf("\nStarting at level %d : Zseed=%+e sg2=%e | %+e<Z<%e\n",level,zseed,sigScale2Ini,zmin,zmax);
  
  vertex vtx;
  vtx.mXYZ[0] = mVtxConstraint[0];
  vtx.mXYZ[1] = mVtxConstraint[1];
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
    res = FitVertex(mVtxTracks,vtx, scaleSigma2, false, zmin, zmax); // fit with current seed and scaling sigma
    if (!res) break;

    float zChange = vtx.mXYZ[2] - oldZ;
    //
    float sigRat = scaleSigma2/scaleSigma2Old;
    printf("<<#%3d Ntacc:%4d Vtx: %+e %+e %+e Sig:%f -> %f %f (Dz:%+.4f)\n",nIter,vtx.mNTracks, vtx.mXYZ[0],vtx.mXYZ[1],vtx.mXYZ[2],
	   scaleSigma2Old,scaleSigma2,sigRat,zChange);
    //
    if (scaleSigma2<0.5f) {
      printf("sigmaScale went below min., stop iterations\n");
      finalize=true;
      break;
    }
    //
    
    if (sigRat<1.0f && sigRat>mStopScaleChange) { // sigma does not drop enough anymore, check convergence
      if ((fabs(zChange)<mMinChangeZ && scaleSigma2<mSigma2Push) || scaleSigma2<mSigma2Accept) { // converged, finalize the vertex
	finalize=true;
	break;
      }
      /*
      else if (scaleSigma2<mSigma2Push) { // decrease sigma to get rid of outliers ???
	scaleSigma2 = (mSigma2Accept+scaleSigma2)*0.5;
	if (scaleSigma2<1.f) scaleSigma2 = 1.f; 
	printf("-->pushing sigma to %f\n",scaleSigma2);
	continue;
      }
      */
      else if (vtx.mNTracks>(mMinTracksPerVtx+mMinTracksPerVtx)) { // stuck between 2 attractors?
	double wghSide[2]={0},wghZSide[2]={0}; //
	int ntSide[2] = {0};
	for (int itr=ntr;itr--;) {
	  vtxTrack &trc = tracks[itr];
	  if (trc.IsUsed() || trc.mZ>zmax || trc.mZ<zmin) continue; // the track is invalidated, skip
	  float wghT = trc.mWgh;
	  if (wghT>0) {
	    int ind = (trc.mZ<vtx.mXYZ[2]) ? 0 : 1; // weight to the left and right of current vertex
	    wghSide[ind] += wghT;
	    wghZSide[ind] += wghT*trc.mZ;
	    ntSide[ind]++;
	  }
	}
	printf("Split to L[Z=%+e N=%3d W:%e] R[Z=%+e N=%3d W=%e]\n",
	       wghSide[0]>1e-6 ? wghZSide[0]/wghSide[0] : -999,ntSide[0], wghSide[0],
	       wghSide[1]>1e-6 ? wghZSide[1]/wghSide[1] : -999,ntSide[1], wghSide[1]);
	//
	{
	  int iside=0;
	  if (wghSide[iside]>1e-6) {  // FIXME constant
	    float zseedSide = wghZSide[iside]/wghSide[iside];
	    printf("Doing %s with Z=%f Sig2Scl=%f\n",iside ? "Right":"Left",zseedSide, scaleSigma2);
	    FindNextVertex(tracks,zseedSide,scaleSigma2/2.f, zmin,vtx.mXYZ[2]);
	  }
	  //
	  iside=1;
	  if (wghSide[iside]>1e-6) {  // FIXME constant
	    float zseedSide = wghZSide[iside]/wghSide[iside];
	    printf("Doing %s with Z=%f Sig2Scl=%f\n",iside ? "Right":"Left",zseedSide, scaleSigma2);
	    FindNextVertex(tracks,zseedSide,scaleSigma2/2.f, vtx.mXYZ[2],zmax);
	  }
	  printf("Break from Left/Right, res=%d\n",res);
	  break; // stop looping around current seed, res must be 0 at this stage
	  //
	}
	/*
	scaleSigma2 *= 0.5;
	if (scaleSigma2<1.f) scaleSigma2 = 1.f; 
	printf("split sigma to %f\n",scaleSigma2);
	continue;
	*/
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
    if ((res=FitVertex(mVtxTracks,vtx, scaleSigma2, true, zmin,zmax)) ) { // final fit with error extraction
      int vID = mVertices.size();
      mVertices.push_back(vtx);
      for (int itr=ntr;itr--;) {
	vtxTrack &trc = mVtxTracks[itr];
	if (trc.IsUsed() || trc.mWgh==0.f || trc.mZ>zmax || trc.mZ<zmin) continue; // the track is invalidated, skip
	trc.SetUsed();
	trc.mVtxID = vID;
	trc.mWgh = 0.f;
      }
    }
  }
  else if (!res) {
    int ntrAcc = 0;
    for (int itr=ntr;itr--;) {
      vtxTrack &trc = mVtxTracks[itr];
      if (trc.IsUsed() || trc.mWgh==0.f || trc.mZ>zmax || trc.mZ<zmin) continue; // the track is invalidated, skip
      trc.SetUsed(); // actually, should disable, TODO
      trc.mVtxID = vtxTrack::kDiscarded; // flag as discarded
      trc.mWgh = 0.f;
      ntrAcc++;
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
bool MVertexFinder::FindVertices()
{
  int ntr = mVtxTracks.size();
  if (ntr<2) return false;
  
  const float zresChar = 200e-4; // characteristic Z resolution
  float sig2ini = (1.+2.*mZRange)/zresChar, zini=0.0f;
  sig2ini *= sig2ini;
  //
  while (FindNextVertex( mVtxTracks, zini, sig2ini,-mZRange,mZRange)) {}
    
}

/*
//______________________________________________
bool MVertexFinder::FindVertices()
{
  int ntr = mVtxTracks.size();
  if (ntr<2) return false;
  
  const float zresChar = 200e-4; // characterisic Z resolution
  float iterCov = sqrtf(mScalSigma2Start)*zresChar;
  int nbin = int(1.+2.*mZRange/iterCov);
  float binw = 2.*mZRange/nbin;
  //std::array<int, nbin> occZ;
  int occZ[nbin] = {0};  // rough histogram of occupancy in Z to generate the seeds
  for (int itr=ntr;itr--;) {
    vtxTrack &trc = mVtxTracks[itr];
    if (trc.IsUsed()) continue; // the track is invalidated, skip
    int iz = (trc.mZ+mZRange)/binw;
    if (iz<0||iz>=nbin) continue;
    occZ[iz]++;
    //
  }
  //
  while(1) {
    int maxBin = -1, maxOcc = 0;
    for (int ibn=nbin;ibn--;) {
      if (occZ[ibn]>maxOcc) {
	maxBin = ibn;
	maxOcc = occZ[ibn];
      }
    }
    if (maxBin<0) break;
    float zseed = (maxBin+0.5f)*binw - mZRange;
    printf("\nStart from max ZBin %d at z=%f\n",maxBin,zseed);
    for (int i=0;i<nbin;i++) printf("(%d) %d ",i,occZ[i]); printf("\n");

    if (!FindNextVertex(zseed)) {
      printf("disabling bin %d\n",maxBin);
      occZ[maxBin] = -1; // ignore this seed in future
    }
    else { // refill tracks
      for (int ib=nbin;ib--;) if (occZ[ib]>0) occZ[ib]=0; // reset all bins except disabled ones
      for (int itr=ntr;itr--;) {
	vtxTrack &trc = mVtxTracks[itr];
	if (trc.IsUsed()) continue; // the track is invalidated, skip
	int iz = (trc.mZ+mZRange)/binw;
	if (iz<0||iz>=nbin || occZ[iz]<0) continue;
	occZ[iz]++;
	//
      }
      //
    }
  }
    
}
*/
//______________________________________________
bool MVertexFinder::FitVertex(std::vector<MVertexFinder::vtxTrack> &tracks, MVertexFinder::vertex &vtx,
			      float &scaleSigma2, bool fillError, float zmin,float zmax)
{
  const double kTiny = 1e-6;
  const double kTukey2 = 6;
  int ntAcc=0,ntr = tracks.size();
  if (ntr<2) return false;
  float* curVtx = vtx.mXYZ;
  //
  double wghSum=0,wghChi2=0; 
  double cxx=0,cxy=0,cxz=0,cx0=0,cyy=0,cyz=0,cy0=0,czz=0,cz0=0;
  //
  for (int itr=ntr;itr--;) {
    vtxTrack &trc = tracks[itr];
    if (trc.IsUsed() || trc.mZ>zmax || trc.mZ<zmin) continue; // the track is invalidated, skip
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
    float wghT = (1.f-chi2T/mTukey2/scaleSigma2);
    if (wghT<kTiny)  {
      trc.mWgh = 0;
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
  if (ntAcc<2) return false;
  //
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
void MVertexFinder::AddTrack(float x,float y,float z,float sy2,float sz2, float syz, float snp, float tgl, float alp)
{
  // add new track to the pool
  float cs = (1.f-snp)*(1.f+snp);
  if (cs<kAlmost0F) return;
  double detI = sy2*sz2 - syz*syz;
  if (fabs(detI)<kAlmost0) return;
  detI = 1./detI;
  //
  vtxTrack trc;
  //
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
  trc.mFlags = 0;
  trc.mVtxID = vtxTrack::kNoVtx;
  mVtxTracks.push_back(trc);
  return;
}

//______________________________________________
void MVertexFinder::SetConstraint(float x,float y,float z)
{
  mVtxConstraint[0] = x;
  mVtxConstraint[1] = y;
  mVtxConstraint[2] = z;
}

//______________________________________________
void MVertexFinder::SetConstraintError(float sy2,float sz2,float syz)
{
  mVtxConstrErr[0] = sy2;
  mVtxConstrErr[1] = sz2;
  mVtxConstrErr[2] = syz;
}

//______________________________________________
void MVertexFinder::PrintVertices() const
{
  ///< print current vertices
  int nv = mVertices.size();
  for (int iv=0;iv<nv;iv++) {
    const vertex& vtx = mVertices[iv];
    printf("#%2d %+e %+e %+e | Ntracks = %d\n",iv,vtx.mXYZ[0],vtx.mXYZ[1],vtx.mXYZ[2],vtx.mNTracks);
    int id = 0;
    for (int ie=0;ie<3;ie++) {printf("    "); for (int je=ie;je<3;je++) printf("%+e ",vtx.mCov[id++]); printf("\n");}
  }
}

//______________________________________________
void MVertexFinder::PrintTracks() const
{
  ///< print tracks
  int nt = mVtxTracks.size();
  for (int it=0;it<nt;it++) {
    const vtxTrack &trc = mVtxTracks[it];
    printf("#%4d Z:%+e W:%+e tgL:%+5.2f | Vtx: = %d\n",it,trc.mZ, trc.mWgh, trc.mTgL, trc.mVtxID);
  }
}
