#ifndef JunParticle_HH
#define JunParticle_HH
#include <iostream>
#include <string>
#include <TROOT.h>
#include <TVector3.h>
#include <TMath.h>

using namespace std;
class JunParticle
{
public:
  JunParticle();
  JunParticle(string pName,Double_t pEnergy,TVector3 pDirection,Double_t pTime=-99999,int a=0,int z=0,int tele=0);
  JunParticle(string pName,Double_t pEnergy,Double_t pTheta,Double_t pPhi,Double_t pTime=-99999,int a=0,int z=0,int tele=0);
  void SetNote(string tNote);
  virtual ~JunParticle();
  string name;
  TVector3 direction;
  Double_t time;//channel of tdc
  Double_t energy;//*MeV
  Double_t theta;//*deg
  Double_t phi;//*deg
  Double_t des[3];//energy on different Detectors *MeV
  int A;
  int Z;
  int tflag;//l1:1000 l0:100 r0:10 r1:1
  string note;// He4 : t0ssd/t0bb7/t1ssd ; Be9 : t0recoil/t0break/t1more
  ClassDef(JunParticle,1);
};
  ClassImp(JunParticle);
#endif