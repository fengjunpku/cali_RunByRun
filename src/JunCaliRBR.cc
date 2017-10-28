#include "JunCaliRBR.hh"

JunCaliRBR::JunCaliRBR(int runnum)
{
  ps = NULL;//very important !!!!
  runno = runnum;
  //Load EnergyLoss Correction
  ploss = new JunLossCorrection();
  ploss->addDataFile("He4_in_Al.txt","He4InAl");
  ploss->addDataFile("Be9_in_Al.txt","Be9InAl");
  ploss->addDataFile("He4_in_Be9.txt","He4InBe");
  ploss->addDataFile("Be9_in_Be9.txt","Be9InBe");
  //readfile
  char _filename[30];
  sprintf(_filename,"extr%04d.root",runno);
  TString _infilepath(in_dir+_filename);
  pfile = TFile::Open(_infilepath,"READ");
  if(!pfile|| !pfile->IsOpen())
  {
    char err[30];
    sprintf(err," * JunCaliRBR::JunCaliRBR() :Can not open %s !",_infilepath.Data())
;
    MiaoError(err);
  }
  pfile->GetObject("tree",ptree);
  InitTree();
  //initial pars
  m_tname[100] = "l0";
  m_tname[10] ="r0";
  //
  string name_list[] = {"e_beam","k_l0w","b_l0w","k_l0b","b_l0b","k_r0w","b_r0w","k_r0b","b_r0b"};
  par_num = sizeof(name_list)/sizeof(name_list[0]);
  par_name = name_list;
  for(int i=0;i<par_num;i++)
  {
    if(par_name[i][0]=='e') m_Par[par_name[i]] = 65;
    if(par_name[i][0]=='k') m_Par[par_name[i]] = 1;
    if(par_name[i][0]=='b') m_Par[par_name[i]] = 0;
  }
  //outputfile
  sprintf(_filename,"chek%04d.root",runno);
  TString _oufilepath(ou_dir+_filename);
  pfout = new TFile(_oufilepath,"RECREATE");
  hql = new TH1D("hql","hql",100,-5,15);
  gql = new TGraph();
  gql->SetNameTitle("gql","gql");
}

JunCaliRBR::~JunCaliRBR()
{
}

void JunCaliRBR::InitTree()
{
  if(!ptree) MiaoError(" * JunCaliRBR::InitTree() : ptree is NULL !");
  ptree->SetBranchAddress("ps",&ps);
}

void JunCaliRBR::Load()
{
  long numOfEntries = ptree->GetEntries();
  printf(" # cali extr%04d #\n",runno);
  for(long ie=0;ie<numOfEntries;ie++)
  {
    ptree->GetEntry(ie);
    if(ps->_num != 1) continue;
    JunParticle be = ps->ap[0];
    if(be.tflag != 10 && be.tflag != 100) continue;
    //
    
    sInfo si;
    si.flag = be.tflag;
    si.e[0] = be.des[0];
    si.e[1] = be.des[1];
    si.e[2] = 0;
    
    si.th = be.direction.Theta();
    si.ph = be.direction.Phi();
    //
    si.angle = ploss->calAngle(si.th,si.ph,m_tname[si.flag]);
    vsi.push_back(si);
    double et = Calculate(si);
    TVector3 sdir(0,0,1);
    sdir.SetMagThetaPhi(1,si.th,si.ph);
    JunParticle mm = getMM(et,sdir,si.flag);
    if(si.flag==100)
    {
      hql->Fill(mm.energy);
      gql->SetPoint(gql->GetN(),mm.theta,mm.energy);
    }
  }
}

double JunCaliRBR::Loss()
{
  //for hist
  
  //
  return 0; 
}

double JunCaliRBR::Calculate(sInfo _si)
{
  double *dl = NULL;//deadlayer
  double e[2];
  if(_si.flag == 100) 
  {
    dl = DL_l0;
    e[0] = _si.e[0]*m_Par["k_l0w"]+m_Par["b_l0w"];
    e[1] = _si.e[1]*m_Par["k_l0b"]+m_Par["b_l0b"];
  }
  if(_si.flag == 10) 
  {
    dl = DL_r0;
    e[0] = _si.e[0]*m_Par["k_r0w"]+m_Par["b_r0w"];
    e[1] = _si.e[1]*m_Par["k_r0b"]+m_Par["b_r0b"];
  }
  double eb = ploss->GetE(dl,_si.e,2,"Be9InAl",_si.angle);//dead layer loss
  eb = ploss->correctEnergy(2*halfTT/TMath::Cos(_si.th),eb,"Be9InBe");//target loss
  return eb;
}

void JunCaliRBR::Save()
{
  pfout->cd();
  hql->Write();
  gql->Write();
  pfout->Close();
}

JunParticle JunCaliRBR::getMM(double _et,TVector3 _dir,int _flag)
{
  double bEn = m_Par["e_beam"];//*MeV
  double epr = _et;
  TVector3 dirR = TMath::Sqrt(2*Mass_Be9*epr)*_dir;
  TVector3 dir0(0,0,1);
  dir0 = TMath::Sqrt(2*Mass_C13*bEn)*dir0;
  TVector3 dir_recon = dir0 - dirR;
  double ene_recon = bEn - epr - dir_recon*dir_recon/Mass_C13/2.;
  JunParticle MM("mm",ene_recon,_dir);
  MM.tflag = _flag;
  return MM;
}