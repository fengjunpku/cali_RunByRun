#include "JunParMan.hh"

JunParMan* JunParMan::thisParMan=NULL;

JunParMan::JunParMan()
{;}

JunParMan::~JunParMan()
{
  cout<<"End of my ParMan!"<<endl;
}

JunParMan* JunParMan::Instance()
{
  if(thisParMan==NULL)
    thisParMan = new JunParMan();
  return thisParMan;
}

void JunParMan::ReadParFile(string filename)
{
  ifstream fin;
  fin.open(filename,ios::in);
  if(!fin)
  {
    cout<<"Miao! Error: Can not open "<<filename<<endl;
    exit(0);
  }
  string buff;
  while(getline(fin,buff))
  {
    stringstream buffstream;
    buffstream.str(buff);
    int type;
    switch(buff[0])
    {
      case '0':
      {
        //cout<<"raw line: "<<buff<<endl;
        string name;
        double value;
        buffstream>>type>>name>>value;
        //cout<<name<<" "<<value<<endl;
        mPar[name]=value;
        //cout<<name<<" "<<mPar[name]<<endl;
      }
        break;
      case '1':
        {
          //cout<<"raw line: "<<buff<<endl;
          Par_DSSD tmp;
          buffstream>>type>>tmp.name>>tmp.heigh>>tmp.width>>tmp.thick>>tmp.vstrips>>tmp.hstrips>>tmp.distance>>tmp.angle>>tmp.xoff>>tmp.yoff;
          vDSSDPar.push_back(tmp);
        }
        break;
      default:
        break;
    }

  }
  fin.close();
}

double JunParMan::GetPar(string parname)
{
  if(!mPar.count(parname))
  {
    cout<<"Miao! Error: Key: "<<parname<<" is not in par Map!"<<endl;
    exit(0);
  }
    return mPar[parname];
}

vector<string> JunParMan::GetAllDSSDName()
{
  vector<string> alldssdname;
  for(auto iter = vDSSDPar.begin();iter != vDSSDPar.end();++iter)
  {
    alldssdname.push_back((*iter).name);
  }
  return alldssdname;
}

Par_DSSD* JunParMan::GetDSSDPar(string dssdname)
{
  for(auto iter = vDSSDPar.begin();iter != vDSSDPar.end();++iter)
  {
    if((*iter).name==dssdname) return &(*iter);
  }
  return NULL;
}

void JunParMan::PrintAll()
{
  cout<<"# Num of standard Parameters : ";
  cout<<mPar.size()<<endl;
  for(auto iter = mPar.begin();iter != mPar.end();++iter)
  {
    cout<<"# "<<iter->first<<" "<<iter->second<<endl;
  }
  cout<<"# Num of par_dssd Parameters : ";
  cout<<vDSSDPar.size()<<endl;
  for(auto iter = vDSSDPar.begin();iter != vDSSDPar.end();++iter)
  {
    cout<<"# "<<(*iter).name;
    cout<<" "<<(*iter).heigh;
    cout<<" "<<(*iter).width;
    cout<<" "<<(*iter).thick;
    cout<<" "<<(*iter).vstrips;
    cout<<" "<<(*iter).hstrips;
    cout<<" "<<(*iter).distance;
    cout<<" "<<(*iter).angle;
    cout<<" "<<(*iter).xoff;
    cout<<" "<<(*iter).yoff;
    cout<<endl;
  }
}
