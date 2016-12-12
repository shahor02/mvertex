#include "MVertexFinder.h"
#include "AliSymMatrix.h"
#include <math.h>

ClassImp(MVertexFinder)

const float kAlmost0F = 1e-12;
const float kAlmost0  = 1e-16;

//______________________________________________
MVertexFinder::MVertexFinder() :
  mVtxTracks()
  ,mVertices()
  ,mMaxVtxIter(30)
  ,mStopScaleChange(0.8)
  ,mSigma2Accept(3.0)
  ,mSigma2Push(10.)
  ,mTukey2(6)
{
  mVtxConstraint[0]=mVtxConstraint[1]=mVtxConstraint[2]=0.f;
  mVtxConstrErr[0]=mVtxConstrErr[1]=mVtxConstrErr[2]=0.f;
}

//______________________________________________
bool MVertexFinder::FindNextVertex()
{
  //
  const double kTiny = 1e-6;
  const double kTukey2 = 6;
  //
  vertex vtx;
  float* curVtx = vtx.mXYZ;
  for (int i=3;i--;) curVtx[i] = mVtxConstraint[i];
  float curVtxErr[6] = {0.f};
  float scaleSigma2=1e6;       // initial sigma scaling
  //
  int nIter = 0;
  int ntAcc=0,ntr = mVtxTracks.size();
  if (ntr<2) return false;
  AliSymMatrix mat(3);
  //
  bool vtxOK = false;
  
  while(nIter++<mMaxVtxIter) {
    ntAcc = 0;
    double mean=0,rms=0,wghL=0,wghR=0;
    double wghSum=0,wghChi2=0; 
    double cxx=0,cxy=0,cxz=0,cx0=0,cyy=0,cyz=0,cy0=0,czz=0,cz0=0;
    //
    printf(">>#%3d Ntr:%4d Vtx: %+e %+e %+e Sig: %f\n",nIter,ntr, curVtx[0],curVtx[1],curVtx[2], scaleSigma2);
    int ntAv = 0;
    for (int itr=ntr;itr--;) {
      //
      vtxTrack &trc = mVtxTracks[itr];
      if (trc.IsUsed()) continue; // the track is invalidated, skip
      //
      ntAv++;
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
	trc.mWgh = 0.f;
	continue;
      }
      wghSum  += wghT;
      wghChi2 += wghT*chi2T;
      //
      double wdz = dz*wghT;
      mean += wdz;
      if (dz<0) wghL += wghT; // sum of weights to the left and right of current vertex
      else      wghR += wghT; // to push it in case it is stuck
      //
      rms  += wdz*wdz;
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
    if (ntAcc<2) break;   // failed
    //
    double vec[3] = {cx0,cy0,cz0};
    mat(0,0) = cxx;
    mat(0,1) = cxy;
    mat(0,2) = cxz;
    mat(1,1) = cyy;
    mat(1,2) = cyz;
    mat(2,2) = czz;
    // 
    if (!mat.SolveChol(vec,true)) return false;
    //
    float scaleSigma2New = wghChi2/wghSum;
    for (int i=3;i--;) curVtx[i] = float(vec[i]);
    //
    float sigRat = scaleSigma2New/scaleSigma2;
    printf("<<#%3d Ntacc:%4d (%4d) Vtx: %+e %+e %+e Sig:%f -> %f %f\n",nIter,ntAcc,ntAv, curVtx[0],curVtx[1],curVtx[2],
	   scaleSigma2,scaleSigma2New,sigRat);
    //
    if (sigRat>mStopScaleChange) { // sigma does not drop enough anymore, check convergence
      if (scaleSigma2<mSigma2Accept) { // converged, finalize the vertex
	if (scaleSigma2>1) { // do one more iteration with nominal errors
	  printf("Forcing sigma to %f\n",scaleSigma2);
	  scaleSigma2 = 1.;
	  continue;
	}
	vtxOK = true;
	break;
      }
      else if (scaleSigma2<mSigma2Push) { // decrease sigme to get rid of outliers
	scaleSigma2 = (mSigma2Accept+scaleSigma2)*0.5;
	printf("-->pushing sigma to %f\n",scaleSigma2);
      }
      else { // no convergence, push toward stronger side
	rms /= wghSum;
	mean /= wghSum;
	rms -= mean*mean;
	rms = rms>0 ? sqrt(rms) : 0.0f;
	printf("-->pushing vertex by %f (%.3f : %.3f)\n",wghL>wghR ? -rms : rms,wghL,wghR);
	curVtx[2] += wghL>wghR ? -rms : rms;
	scaleSigma2 *= 0.5f;
      }
    }
    else scaleSigma2 = scaleSigma2New;
    //
  }  
  if (vtxOK) {
    // flag tracks attached to the vertex
    vtx.mCov[0] = mat(0,0);
    vtx.mCov[1] = mat(0,1);
    vtx.mCov[2] = mat(1,1);
    vtx.mCov[3] = mat(0,2);
    vtx.mCov[4] = mat(1,2);
    vtx.mCov[5] = mat(2,2);
    vtx.mNTracks = ntAcc;
    mVertices.push_back(vtx);
    for (int itr=ntr;itr--;) {
      vtxTrack &trc = mVtxTracks[itr];
      if (trc.IsUsed() || trc.mWgh<kTiny) continue; // the track is invalidated, skip
      trc.SetUsed();
    }
  }
  //
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
  trc.flags = 0;
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
