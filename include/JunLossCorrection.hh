#ifndef JunLossCorrection_HH
#define JunLossCorrection_HH

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <TMVA/TSpline2.h>
#include <TMath.h>
#include <TVector3.h>

#include "JunErrors.hh"
#include "reDefine.hh"

using namespace std;
using namespace reDef;
using namespace TMVA;


class JunLossCorrection
{
public:
  JunLossCorrection();
  virtual ~JunLossCorrection();
  void addDataFile(string filename,string discription);
  double correctEnergy(const double range,const double energy,string discription);
  double GetE(const double *range,const double *energys,const int ne,string discription,double theta=0);
  //unit: um, MeV,rad
  double calAngle(double th,double ph,string tname);
private:
  map<string,TSpline2*> mRvE;//range vs. energy
  map<string,TSpline2*> mEvR;//energy vs. range
};
#endif