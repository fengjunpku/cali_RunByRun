#include "JunLossCorrection.hh"

JunLossCorrection::JunLossCorrection()
{}

void JunLossCorrection::addDataFile(string filename,string discription)
{
  ifstream fin;
  string datadir = "../rebuild_13C/rebuild_v0/srim_data/";
  fin.open(datadir+filename,ios::in);
  if(!fin) MiaoError("Can not Open : "+filename);
  string buff;
  TGraph *gr = new TGraph();
  TGraph *rg = new TGraph();
  int n=0;
  while(getline(fin,buff))
  {
    stringstream buffstream;
    buffstream.str(buff);
    if(buff[0]=='#') continue;
    double energy,tmp1,tmp2,range;
    string unit1,unit2;
    buffstream>>energy>>unit1>>tmp1>>tmp2>>range>>unit2;
    if(unit1=="keV") energy=energy/1000.;//keV->MeV
    if(unit2=="A") range = range/10000.;//A->um
    if(unit2=="mm") range = range*1000.;//mm->mm
    gr->SetPoint(n,energy,range);
    rg->SetPoint(n,range,energy);
    n++;
  }
  mRvE[discription] = new TSpline2(discription.c_str(),gr);
  mEvR[discription] = new TSpline2(discription.c_str(),rg);
}

double JunLossCorrection::correctEnergy(const double range,const double energy,string discription)
{
  if(!mRvE.count(discription)||!mEvR.count(discription))
    MiaoError("Undefined "+discription+" as index.");
  double r0 = mRvE[discription]->Eval(energy);
  return mEvR[discription]->Eval(r0+range);
}

double JunLossCorrection::GetE(const double *range,const double *energys,const int ne,string discription,double theta)
{
  double re = 0;
  double costh = TMath::Cos(theta);
  for(int ine=0;ine<ne;ine++)
    re = correctEnergy(range[ne-ine-1]/costh,re+energys[ne-ine-1],discription);
  return re;
}

double JunLossCorrection::calAngle(double th,double ph,string tname)
{
  double deg = TMath::DegToRad();
  const double *angles = NULL;
  if("l0" == tname) angles = angles_l0;
  if("l1" == tname) angles = angles_l1;
  if("r0" == tname) angles = angles_r0;
  if("r1" == tname) angles = angles_r1;
  TVector3 vd,vp;
  vp.SetMagThetaPhi(1,th,ph);
  vd.SetMagThetaPhi(1,angles[0]*deg,angles[1]*deg);
  return vp.Angle(vd);
}

JunLossCorrection::~JunLossCorrection()
{}