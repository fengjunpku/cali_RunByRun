#ifndef JunParMan_HH
#define JunParMan_HH 1

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <cstdlib>
using namespace std;

struct Par_DSSD
{
  string name;
  double heigh;
  double width;
  double thick;
  int vstrips;
  int hstrips;
  double distance;
  double angle;
  double xoff;
  double yoff;
};

class JunParMan
{
private:
  JunParMan();
  static JunParMan* thisParMan;
public:
  virtual ~JunParMan();
  static JunParMan* Instance();
  void ReadParFile(string filename);
  double GetPar(string parname);
  Par_DSSD* GetDSSDPar(string dssdname);
  void PrintAll();
  vector<string> GetAllDSSDName();
private:
  map<string,double> mPar;
  vector<Par_DSSD> vDSSDPar;
};
#endif
