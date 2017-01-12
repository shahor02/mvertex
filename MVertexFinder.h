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
    bool CanUse(float zmin,float zmax) const {return mVtxID==kNoVtx && mZ>zmin && mZ<zmax;}
  };

  struct vertex {
    float mXYZ[3];      ///< vertex position
    float mCov[6];      ///< sym.cov. matrix in lower triangle representation
    int   mNTracks;     ///< number of tracks associated
  vertex() : /*mXYZ{0.f},mCov{0.f},*/mNTracks(0) { memset(this,0,sizeof(vertex));}
  };
  
  MVertexFinder();
  void Reset();
  bool FindVertices();
  bool FindNextVertex(std::vector<vtxTrack> &tracks,float zseed,float sigScale2Ini, float zmin,float zmax);
  bool FitVertex(std::vector<vtxTrack> &tracks, vertex &vtx, float &scaleSigma2, bool fillError, float zmin=-999.f,float zmax=999.f);
  void AddTrack(float x,float y,float z,float sy2,float sz2, float syz, float snp, float tgl, float alp);
  void SetConstraint(float x,float y,float z);
  void SetConstraintError(float sy2,float sz2,float syz);
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
  
  std::vector<vtxTrack> mVtxTracks;           ///< container for input tracks
  std::vector<vertex> mVertices;              ///< container for found vertices
  int   mMaxVtxIter;                          ///< max number of iterations per vertex
  int   mMinTracksPerVtx;                     ///< min number of tracks per vertex
  float mMinChangeZ;                          ///< stop if Z changes by less than this amount
  float mStopScaleChange;                     ///< stopping condition: max sigma2New/sigma2
  float mSigma2Accept;                        ///< stopping condition: acceptable Sigma2
  float mSigma2Push;                          ///< trigger algo push if iterations stuck at sigma2 above this
  float mTukey2I;                             ///< 1./[Tukey parameter]^2
  float mZRange;                              ///< ZRange to accept
  float mVtxConstraint[3];                    ///< nominal vertex constraint
  float mVtxConstrErr[3];                     ///< nominal vertex constraint errors
  //
  static const float kDefTukey;               ///< def.value for tukey constant
  static const float kHugeF;                  ///< very large float
  static const float kAlmost0F;               ///< tiny float
  static const double kAlmost0D;              ///< tiny double
  //
  ClassDef(MVertexFinder,1)
};

#endif
