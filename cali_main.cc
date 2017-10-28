#include <iostream>
#include <string>
#include "stdio.h"
#include "stdlib.h"

#include <TStopwatch.h>
#include <TString.h>

#include "reDefine.hh"
#include "JunErrors.hh"
#include "JunParMan.hh"

#include "JunParticle.hh"
#include "JunLossCorrection.hh"
#include "JunCaliRBR.hh"
using namespace std;
using namespace reDef;

//double three peaks
double PP_Be[3] = {0,3.68,7.55};
//thickness of dead layer [um]
double DL_l0[3] = {1.25,1.1,1.3};//w1 bb7 ssd
double DL_r0[3] = {0.75,0.8,0.7};//w1 bb7 ssd
double DL_l1[2] = {1,1};//w1 ssd
double DL_r1[2] = {1,1};//w1 ssd



int main(int argc,char** argv)
{
  //---------------------------------------------------------------------------
  TStopwatch watch;
  //---------------------------------------------------------------------------
  if(argc != 2)
    MiaoError("Please give a run num !");
  int runnum = atoi(argv[1]);
  //=======================================================
  JunCaliRBR *conFront = new JunCaliRBR(runnum);
  conFront->Load();
  conFront->Save();
  //---------------------------------------------------------------------------
  printf("\n ======End of  %04d ====\n",runnum);
  printf("   ==> %04d CPU_Time: %s%f%s s, RealTime: %s%f%s s\n",runnum,RED,watch.CpuTime(),COLOR_END,RED,watch.RealTime(),COLOR_END);
}