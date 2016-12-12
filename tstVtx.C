
TObjArray genVtx;
MVertexFinder* vtf = 0;

const float sigVZ = 5;
const float sigVR = 60e-4;
const float meanV[3] = {0.03,0.3,0.2};

const float sigTZ = 200e-4;
const float sigTR = 100e-4;


void tstVtx(int nv=1, int ntrMean=20, int nNoise=20)
{
  if (!vtf) {
    vtf = new MVertexFinder();
    vtf->SetConstraint(meanV[0],meanV[1],meanV[2]);
    vtf->SetConstraintError(sigVR*sigVR,sigVR*sigVR,sigVZ*sigVZ);
  }
  vtf->Reset();
  
  genVtx.Delete();
  for (int iv=0;iv<nv;iv++) {
    double vXYZ[3] = {
      vXYZ[0] = gRandom->Gaus(meanV[0],sigVR),
      vXYZ[1] = gRandom->Gaus(meanV[1],sigVR),
      vXYZ[2] = gRandom->Gaus(meanV[2],sigVZ)
    };
    AliESDVertex* vtx = new AliESDVertex();
    vtx->SetXYZ(vXYZ);
    genVtx.Add(vtx);
    //
    int ntr = 0;
    while (ntr<2) ntr = gRandom->Poisson(ntrMean);
    for (int itr=0;itr<ntrMean;itr++) {
      double alp = gRandom->Rndm()*TMath::Pi()*2;
      double sna = TMath::Sin(alp);
      double csa = TMath::Cos(alp);
      float vlocX =  vXYZ[0]*csa+vXYZ[1]*sna;
      float vlocY = -vXYZ[0]*sna+vXYZ[1]*csa;
      //
      double xtr = gRandom->Gaus(vlocX,sigTR);
      double ytr = gRandom->Gaus(vlocY,sigTR);
      double ztr = gRandom->Gaus(vXYZ[2],sigTZ);
      //
      double tgl = (gRandom->Rndm()-0.5)*2;
      //
      vtf->AddTrack(xtr,ytr,ztr, sigTR*sigTR, sigTR*sigTR, 0, 0., tgl,alp);
    }
    vtx->SetNContributors(ntr);
  }
  //
  for (int i=0;i<nNoise;i++) {
    double alp = gRandom->Rndm()*TMath::Pi()*2;
    double sna = TMath::Sin(alp);
    double csa = TMath::Cos(alp);
    float vlocX =  vXYZ[0]*csa+vXYZ[1]*sna;
    float vlocY = -vXYZ[0]*sna+vXYZ[1]*csa;
    //
    double xtr = gRandom->Gaus(vlocX,sigTR);
    double ytr = gRandom->Gaus(vlocY,sigTR);
    double ztr = gRandom->Gaus(meanV[2],sigVZ);
    //
    double tgl = (gRandom->Rndm()-0.5)*2;
    //
    vtf->AddTrack(xtr,ytr,ztr, sigTR*sigTR, sigTR*sigTR, 0, 0., tgl,alp);
  }

  genVtx.Print();
  vtf->FindNextVertex();
}


void PrintGen()
{
  int nv = genVtx.GetEntriesFast();
  for (int i=0;i<nv;i++) {
    AliESDVertex* vt = (AliESDVertex*) genVtx[i];
    printf("#%2d %+e %+e %+e | %d\n",i,vt->GetX(),vt->GetY(),vt->GetZ(),vt->GetNContributors());
  }
}
