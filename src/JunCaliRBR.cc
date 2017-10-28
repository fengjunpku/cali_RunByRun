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
  par_name = new string[par_num];
  for(int i=0;i<par_num;i++)
  {
    par_name[i] = name_list[i];
    if(par_name[i][0]=='e') 
    {
      m_Par[par_name[i]] = 65;
      m_step[par_name[i]] = 1e-2;
    }
    if(par_name[i][0]=='k') 
    {
      m_Par[par_name[i]] = 1;
      m_step[par_name[i]] = 1e-3;
    }
    if(par_name[i][0]=='b') 
    {
      m_Par[par_name[i]] = 0;
      m_step[par_name[i]] = 1e-3;
    }
    m_dPar[par_name[i]] = 0;
  }
  //outputfile
  sprintf(_filename,"chek%04d.root",runno);
  TString _oufilepath(ou_dir+_filename);
  pfout = new TFile(_oufilepath,"RECREATE");
  //
  hql = new TH1D("hql","hql",100,-5,15);
  gql = new TGraph();gql->SetNameTitle("gql","gql");
  hqr = new TH1D("hqr","hqr",100,-5,15);
  gqr = new TGraph();gqr->SetNameTitle("gqr","gqr");
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
  }
  //Update();
  Loop();
}

double JunCaliRBR::getLoss(TH1D *hd,TGraph *gd)
{
  //for hist
  double p[9];
  TF1 g1("g1","gaus",-1,1);
  TF1 g2("g2","gaus",2,5);
  TF1 g3("g3","gaus",5.5,9);
  TF1 gt("gt","gaus(0)+gaus(3)+gaus(6)",-2,10);
  hd->Fit(&g1,"QNR");
  hd->Fit(&g2,"QNR+");
  hd->Fit(&g3,"QNR+");
  g1.GetParameters(&p[0]);
  g2.GetParameters(&p[3]);
  g3.GetParameters(&p[6]);
  gt.SetParameters(p);
  hd->Fit(&gt,"QNR+");
  gt.GetParameters(&p[0]);
  double cut1 = p[1]+p[2]*(p[4]-p[1])/(p[5]+p[2]);
  double cut2 = p[4]+p[5]*(p[7]-p[4])/(p[8]+p[5]);
  //cout<<cut1<<" "<<cut2<<endl;
  //cout<<p[0]<<" "<<p[3]<<" "<<p[6]<<endl;
  //cout<<p[1]<<" "<<p[4]<<" "<<p[7]<<endl;
  //cout<<p[2]<<" "<<p[5]<<" "<<p[8]<<endl;
  //
  TGraph p1;
  TGraph p2;
  TGraph p3;
  int count = gd->GetN();
  for(int i=0;i<count;i++)
  {
    double x,y;
    gd->GetPoint(i,x,y);
    if(y>cut2)
    {
      p3.SetPoint(p3.GetN(),x,y);
    }
    else if(y>cut1)
    {
      p2.SetPoint(p2.GetN(),x,y);
    }
    else
    {
      p1.SetPoint(p1.GetN(),x,y);
    }
  }
  double k[3];
  TF1 line("line","pol1");
  p1.Fit(&line,"Qrob=0.75");
  k[0] = line.GetParameter(1);
  p2.Fit(&line,"Qrob=0.75");
  k[1] = line.GetParameter(1);
  p3.Fit(&line,"Qrob=0.75");
  k[2] = line.GetParameter(1);
  //cout<<k[0]<<" "<<k[1]<<" "<<k[2]<<endl;
  //
  double loss = 0;
  for(int i=0;i<3;i++)
  {
    double cst = 100;//adjust loss function
    double part = cst*(p[3*i+1]-PP_Be[i])*k[i]/p[3*i+2];
    loss += (part*part);
  }
  //cout<<loss<<endl;
  return loss; 
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
  double eb = ploss->GetE(dl,e,2,"Be9InAl",_si.angle);//dead layer loss
  eb = ploss->correctEnergy(2*halfTT/TMath::Cos(_si.th),eb,"Be9InBe");//target loss
  return eb;
}

void JunCaliRBR::Loop()
{
  double eps = 100;
  double loss_lst = 0;
  double loss_now = 0;
  while(eps>1e-2)
  {
    //update pars
    cout<<"==========="<<endl;
    for(int i=0;i<par_num;i++)
    {
      m_Par[par_name[i]] += m_dPar[par_name[i]];
      cout<<m_Par[par_name[i]]<<" ";
    }
    cout<<endl;
    Update();
    loss_now = getLoss(hql,gql)+getLoss(hqr,gqr);
    //calculate d_loss
    for(int i=0;i<par_num;i++)
    {
      double tmp = m_Par[par_name[i]];
      m_Par[par_name[i]] += m_step[par_name[i]];
      Update();
      double loss = getLoss(hql,gql)+getLoss(hqr,gqr);
      //cout<<par_name[i]<<" : "<<loss_now-loss<<endl;
      m_Par[par_name[i]] = tmp;
      m_dPar[par_name[i]] = (loss_now-loss)*0.1*m_step[par_name[i]];
      //cout<<par_name[i]<<" : "<<m_dPar[par_name[i]]<<endl;
    }
    eps = TMath::Abs(loss_now - loss_lst);
    loss_lst = loss_now;
    cout<<eps<<" "<<loss_now<<endl;
  }
}

void JunCaliRBR::Update()
{
  //update hist & graph
  if(gql->GetN()>0) gql->Set(0);
  if(gqr->GetN()>0) gqr->Set(0);
  if(hql->GetSum()>0) hql->Reset();
  if(hqr->GetSum()>0) hqr->Reset();
  int count = vsi.size();
  for(int i=0;i<count;i++)
  {
    double et = Calculate(vsi[i]);
    TVector3 sdir(0,0,1);
    sdir.SetMagThetaPhi(1,vsi[i].th,vsi[i].ph);
    JunParticle mm = getMM(et,sdir,vsi[i].flag);
    if(vsi[i].flag==100)
    {
      hql->Fill(mm.energy);
      gql->SetPoint(gql->GetN(),mm.theta,mm.energy);
    }
    if(vsi[i].flag==10)
    {
      hqr->Fill(mm.energy);
      gqr->SetPoint(gqr->GetN(),mm.theta,mm.energy);
    }
  }
}

void JunCaliRBR::Save()
{
  pfout->cd();
  hql->Write();
  hqr->Write();
  gql->Write();
  gqr->Write();
  pfout->Close();
  //
  FILE *txt;
  char buff[20];
  sprintf(buff,"out/cprbr%04d.txt",runno);
  txt = fopen(buff,"w+");
  if(!txt) MiaoError("JunCaliRBR::Save(), Can't open output file !");
  for(int i=0;i<par_num;i++)
  {
    fprintf(txt," %f ",m_Par[par_name[i]]);
  }
  fprintf(txt,"\n");
  fclose(txt);
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