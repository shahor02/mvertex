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
    enum {kUsed};
    void SetBit(int i)        {flags |= 0x1<<i;}
    bool TestBit(int i) const {return flags&(0x1<<i);}
    //
    bool IsUsed()       const {return TestBit(kUsed);}
    void SetUsed()            {SetBit(kUsed);}
    
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
    char  flags;        ///< status bits
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
  bool FindNextVertex(float zseed=0);
  bool FindNextVertex(std::vector<vtxTrack> &tracks,float zseed,float sigScale2Ini);
  bool FitVertex(std::vector<vtxTrack> &tracks, vertex &vtx, float &scaleSigma2, bool fillError);
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
  void  PrintVertices() const;
  //
 protected:
  std::vector<vtxTrack> mVtxTracks;           ///< container for input tracks
  std::vector<vertex> mVertices;              ///< container for found vertices
  int   mMaxVtxIter;                          ///< max number of iterations per vertex
  float mScalSigma2Start;                     ///< initial value of scaling sigma^2
  float mMinChangeZ;                          ///< stop if Z changes by less than this amount
  float mStopScaleChange;                     ///< stopping condition: max sigma2New/sigma2
  float mSigma2Accept;                        ///< stopping condition: acceptable Sigma2
  float mSigma2Push;                          ///< trigger algo push if iterations stuck at sigma2 above this
  float mTukey2;                              ///< Tukey parameter^2
  float mZRange;                              ///< ZRange to accept
  float mVtxConstraint[3];                    ///< nominal vertex constraint
  float mVtxConstrErr[3];                     ///< nominal vertex constraint errors
  //
  ClassDef(MVertexFinder,1)
};

#endif
