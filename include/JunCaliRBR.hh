#ifndef JunCaliRBR_HH
#define JunCaliRBR_HH 1

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TGraph.h>
#include <TF1.h>
#include <vector>
#include <map>

#include "JunLossCorrection.hh"
#include "JunArrayOfParticle.hh"

//
extern double PP_Be[3];
//thickness of dead layer [um]
extern double DL_l0[3];//w1 bb7 ssd
extern double DL_r0[3];//w1 bb7 ssd
extern double DL_l1[2];//w1 ssd
extern double DL_r1[2];//w1 ssd


struct sInfo
{
  int flag;
  double e[3];
  double th;
  double ph;
  double angle;
};

class JunCaliRBR
{
public:
  JunCaliRBR(int runnum);
  virtual ~JunCaliRBR();
  void Load();
  void Save();

private:
  void InitTree();
  void Loop();
  void Update();
  double Calculate(sInfo _si);
  double getLoss(TH1D *hd,TGraph *gd);
  JunParticle getMM(double _et,TVector3 _dir,int _flag);
  JunArrayOfParticle *ps;
  int runno;
  JunLossCorrection *ploss;
  TFile *pfile;
  TFile *pfout;
  TTree *ptree;
  TH1D *hql;
  TH1D *hqr;
  TGraph *gql;
  TGraph *gqr;
  vector<sInfo> vsi;
  string *par_name;
  int par_num;
  map<string,double> m_Par;
  map<string,double> m_dPar;
  map<string,double> m_step;
  map<int,string> m_tname;
};

#endif