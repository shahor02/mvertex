#ifndef MVERTEXFINDER_H
#define MVERTEXFINDER_H

#include <TObject.h>
#include <vector>

class MVertexFinder : public TObject
{
 public:
  
  struct vtxTrack {
    /** Straight track parameterization in the frame defined by alpha angle.
	Assumed to be defined in the proximity to vertex, so that the straight-line
	extrapolation Y=mY+mTgP*(x-mX) and Z=mZ+mTgL*(x-mX) is precise
    */
    enum {kUsed, kNoVtx=-1,kDiscarded=kNoVtx-1};
    unsigned int mStamp;///< track stamp (time, event id, etc)
    float mX;           ///< reference X
    float mY;           ///< Y at X
    float mZ;           ///< Z at X
    float mSig2YI;      ///< YY component of inverse cov.matrix
    float mSig2ZI;      ///< ZZ component of inverse cov.matrix
    float mSigYZI;      ///< YZ component of inverse cov.matrix
    float mTgP;         ///< tangent(phi) in tracking frame
    float mTgL;         ///< tangent(lambda)
    float mCosAlp;      ///< cos of alpha frame
    float mSinAlp;      ///< sin of alpha frame
    float mWgh;         ///< track weight wrt current vertex seed
    short mVtxID;       ///< assigned vertex
    //
    bool CanUse()                      const {return mVtxID==kNoVtx;}
    bool CanUse(float zmin,float zmax) const {return mVtxID==kNoVtx && mZ>zmin && mZ<zmax;}
    bool operator < (const vtxTrack& trc) const {return mZ<trc.mZ;}
  };

  struct vertex {
    float mXYZ[3];      ///< vertex position
    float mCov[6];      ///< sym.cov. matrix in lower triangle representation
    float mChi2;        ///< total chi^2
    int   mNTracks;     ///< number of tracks associated
    unsigned int mStamp;///< vertex stamp (time, event id, etc)
  vertex() : /*mXYZ{0.f},mCov{0.f},*/mChi2(0.f),mNTracks(0),mStamp(0) { memset(this,0,sizeof(vertex));}
  };
  
  MVertexFinder();
  void Reset();
  bool FindVertices();
  bool FindNextVertex(std::vector<vtxTrack> &tracks,float zseed,float sigScale2Ini, float zmin,float zmax);
  bool FindNextVertex(std::vector<int> &trcID,float zseed,float sigScale2Ini);
  
  bool FitVertex(std::vector<vtxTrack> &tracks, vertex &vtx, float &scaleSigma2, bool fillError, float zmin=-999.f,float zmax=999.f);
  bool FitVertex(std::vector<int> &trcITs, vertex &vtx, float &scaleSigma2, bool fillError);
  void AddTrack(float x,float y,float z,float sy2,float sz2, float syz, float snp, float tgl, float alp, unsigned int stamp=0);

  bool CheckVertexTrackStamps(vertex& vtx, const vtxTrack& trc);
  bool IsStampDefined(UInt_t stamp) const {return stamp;}
  
  /// Interaction region settings
  void SetIRPos(float x=0.f,float y=0.f, float z=0.f);
  void SetIRSig2(float sigY2=1.0f,float sigZ2=1.0f, float sigYZ=0.f);
  //
  bool  GetUseConstraint() const                   {return mUseConstraint;}
  void  SetUseConstraint(bool v=true)              {mUseConstraint = v;}
  //
  void  SetMinChangeZ(float v)                     {mMinChangeZ = v;}
  float GetMinChangeZ() const                      {return mMinChangeZ;}
  //
  void  SetZRange(float z)                         {mZRange = z;}
  float GetZRange()     const                      {return mZRange;}
  //
  void  SetTukey(float t)                          {mTukey2I = t>0.f ? 1.f/(t*t) : 1.f/(kDefTukey*kDefTukey);}
  float GetTukey()      const;
  //
  void  PrintVertices() const;
  void  PrintTracks()   const;
  //
 protected:

  void LRAttractors(const std::vector<vtxTrack> &tracks,float zmn,float zmx,float currZ,float zLR[2]) const;
  bool LRAttractors(const std::vector<int> &trcID,float currZ, std::vector<int> tgt[2],float zLR[2]) const;
  void SelectFreeTracks(const std::vector<int> &src,std::vector<int> &tgt,float zmn,float zmx) const;
  void DisableTracks(const std::vector<int> &src);  
  
  std::vector<vtxTrack> mVtxTracks;           ///< container for input tracks
  std::vector<vertex> mVertices;              ///< container for found vertices
  bool  mUseZSorting;                         ///< optionally presort tracks in Z
  bool  mUseConstraint;                       ///< impose mean-vertex constraint

  int   mRemovedTracks;                       ///< number of tracks removed from the pool
  int   mMaxVtxIter;                          ///< max number of iterations per vertex
  int   mMinTracksPerVtx;                     ///< min number of tracks per vertex
  float mMinChangeZ;                          ///< stop if Z changes by less than this amount
  float mStopScaleChange;                     ///< stopping condition: max sigma2New/sigma2
  float mSigma2Accept;                        ///< stopping condition: acceptable Sigma2
  float mTukey2I;                             ///< 1./[Tukey parameter]^2
  float mZRange;                              ///< ZRange to accept
  float mIRPos[3];                            ///< nominal vertex constraint
  float mIRSig2I[3];                          ///< nominal vertex constraint inverted errors^2
  //
  static const float kDefTukey;               ///< def.value for tukey constant
  static const float kHugeF;                  ///< very large float
  static const float kAlmost0F;               ///< tiny float
  static const double kAlmost0D;              ///< tiny double
  //
  ClassDef(MVertexFinder,1)
};

//______________________________________________
inline void MVertexFinder::SetIRPos(float x,float y, float z)
{
  /// set mean vertex position
  mIRPos[0] = x;
  mIRPos[1] = y;
  mIRPos[2] = z;
}

//______________________________________________
inline void MVertexFinder::SetIRSig2(float sgx2,float sgy2, float sgz2)
{
  mIRSig2I[0] = sgx2>kAlmost0F ? 1.f/sgx2 : 0.f;
  mIRSig2I[1] = sgy2>kAlmost0F ? 1.f/sgy2 : 0.f;
  mIRSig2I[2] = sgz2>kAlmost0F ? 1.f/sgz2 : 0.f;  
}

//______________________________________________
inline void MVertexFinder::DisableTracks(const std::vector<int> &src)
{
  // flag tracks as disabled
  int ntrAcc = 0;
  for (int it=src.size();it--;) {
    vtxTrack &trc = mVtxTracks[src[it]];
    if (trc.CanUse() && trc.mWgh>0.f) {
      trc.mVtxID = vtxTrack::kDiscarded; // flag as discarded
      trc.mWgh = 0.f;
      ntrAcc++; // tmp
      mRemovedTracks++;
    }
  }
  printf("Disabling %d tracks after failure\n",ntrAcc);
}


#endif
