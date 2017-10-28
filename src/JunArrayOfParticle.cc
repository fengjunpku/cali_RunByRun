#include "JunArrayOfParticle.hh"

JunArrayOfParticle::JunArrayOfParticle()
{
  _num = 0;
  _num_he4 = 0;
  _num_be9 = 0;
  _code = 0;
}

int JunArrayOfParticle::GetNum(int Z,int A,int tele)
{
  int cnt = 0;
  for(int i=0;i<_num;i++)
  {
    if(ap[i].Z==Z && ap[i].A==A && ap[i].tflag==tele)
      cnt++;
  }
  return cnt;
}

int JunArrayOfParticle::GetNum(int Z,int A)
{
  int cnt = 0;
  for(int i=0;i<_num;i++)
  {
    if(ap[i].Z==Z && ap[i].A==A)
      cnt++;
  }
  return cnt;
}

int JunArrayOfParticle::GetNum(int Z,int A,const string note)
{
  int cnt = 0;
  for(int i=0;i<_num;i++)
  {
    if(ap[i].Z==Z && ap[i].A==A && ap[i].note == note)
      cnt++;
  }
  return cnt;
}

JunParticle *JunArrayOfParticle::GetParticle(int Z,int A,int tele)
{
  int numP = GetNum(Z,A,tele);
  if(0 == numP) return NULL;
  JunParticle *stp = new JunParticle[numP];
  int k = 0;
  for(int i=0;i<_num;i++)
  {
    if(ap[i].Z==Z && ap[i].A==A && ap[i].tflag==tele)
      stp[k++] = ap[i];
  }
  return stp;
}

JunParticle *JunArrayOfParticle::GetParticle(int Z,int A)
{
  int numP = GetNum(Z,A);
  if(0 == numP) return NULL;
  JunParticle *stp = new JunParticle[numP];
  int k = 0;
  for(int i=0;i<_num;i++)
  {
    if(ap[i].Z==Z && ap[i].A==A)
      stp[k++] = ap[i];
  }
  return stp;
}

JunParticle *JunArrayOfParticle::GetParticle(int Z,int A,const string note)
{
  int numP = GetNum(Z,A,note);
  if(0 == numP) return NULL;
  JunParticle *stp = new JunParticle[numP];
  int k = 0;
  for(int i =0;i<_num;i++)
  {
    if(ap[i].Z==Z && ap[i].A==A && ap[i].note==note)
      stp[k++] = ap[i];
  }
  return stp;
}

void JunArrayOfParticle::Add(JunParticle tp)
{
  ap[_num%10] = tp;
  _num++;
  _code += tp.tflag;
  if(tp.Z == 2 && tp.A == 4) _num_he4++;
  if(tp.Z == 4 && tp.A == 9 && tp.note != "t0more") _num_be9++;
}

JunArrayOfParticle::~JunArrayOfParticle()
{}