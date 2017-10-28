#include "JunParticle.hh"


JunParticle::JunParticle()
{
  name      = "null";
  direction.SetXYZ(0,0,-1);
  theta     = TMath::QuietNaN();
  phi       = TMath::QuietNaN();
  energy    = TMath::QuietNaN();
  time      = TMath::QuietNaN();
  des[0]    = TMath::QuietNaN();
  des[1]    = TMath::QuietNaN();
  des[2]    = TMath::QuietNaN();
  A         = 0;
  Z         = 0;
  tflag     = 0;
  note      = "null";
}

JunParticle::JunParticle(string pName,Double_t pEnergy,TVector3 pDirection,Double_t pTime,int a,int z,int tele)
{
  name      = pName;
  direction = pDirection.Unit();
  energy    = pEnergy;
  theta     = direction.Theta()*TMath::RadToDeg();
  phi       = direction.Phi()*TMath::RadToDeg();
  time      = pTime;
  A         = a;
  Z         = z;
  tflag     = tele;
  note      = "null";
}

JunParticle::JunParticle(string pName,Double_t pEnergy,Double_t pTheta,Double_t pPhi,Double_t pTime,int a,int z,int tele)
{
  name      = pName;
  energy    = pEnergy;
  direction.SetMagThetaPhi(1,pTheta,pPhi);
  theta     = direction.Theta()*TMath::RadToDeg();
  phi       = direction.Phi()*TMath::RadToDeg();
  time      = pTime;
  A         = a;
  Z         = z;
  tflag     = tele;
  note      = "null";
}

void JunParticle::SetNote(string tNote)
{
  note = tNote;
}

JunParticle::~JunParticle()
{}