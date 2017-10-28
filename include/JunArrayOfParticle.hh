#ifndef JunArrayOfParticle_HH
#define JunArrayOfParticle_HH 1

#include "JunParticle.hh"
#include "JunParMan.hh"

class JunArrayOfParticle
{
public:
  JunArrayOfParticle();
  virtual ~JunArrayOfParticle();
  JunParticle ap[10];
  int GetNum(int Z,int A,int tele);
  int GetNum(int Z,int A);
  int GetNum(int Z,int A,const string note);
  void Add(JunParticle tp);
  JunParticle* GetParticle(int Z,int A,int tele);
  JunParticle* GetParticle(int Z,int A);
  JunParticle* GetParticle(int Z,int A,const string note);
  int _num;
  int _num_he4;
  int _num_be9;
  int _code;
};

#endif