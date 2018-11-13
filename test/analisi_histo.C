void analisi_histo(string directory, string channel, double l){
 
  //gStyle->SetOptStat(0000);
  
  string name_file_data = (directory + "/treesTTJets.root").c_str();
  string name_file_histo = ("histo_" + channel + ".root").c_str();
  
  TFile *f = new TFile(name_file_data.c_str());
  TTree * fChain = new TTree;
  fChain = (TTree*) (f->Get("DMTreesDumper/ttDM__noSyst"));
  
  string pre_sel = "(muonsLoose_size==0 && electronsVeto_size==0 && jetsAK8CHS_size>=1 && ((jetsAK8CHS_tau3CHS[0]/jetsAK8CHS_tau2CHS[0]<0.81 &&jetsAK8CHS_prunedMassCHS[0]>105)))";

  std::vector<TString> var,name,sel,presel;
  std::vector<double> nbins, Min, Max;
  
  var.push_back("jetsAK4CHS_size");name.push_back("jetsAK4CHS_size"); nbins.push_back(60); Min.push_back(0); Max.push_back(60); sel.push_back(""); presel.push_back(pre_sel); 
  var.push_back("jetsAK4CHSTight_size");name.push_back("jetsAK4CHSTight_size"); nbins.push_back(32); Min.push_back(0); Max.push_back(16); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK4CHS_E[0]");name.push_back("jetsAK4CHS_E"); nbins.push_back(15); Min.push_back(0); Max.push_back(3000); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK4CHSTight_E[0]");name.push_back("jetsAK4CHSTight_E"); nbins.push_back(15); Min.push_back(0); Max.push_back(3000); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK4CHS_Pt[0]");name.push_back("jetsAK4CHS_Pt"); nbins.push_back(25); Min.push_back(0); Max.push_back(1500); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK4CHSTight_Pt[0]");name.push_back("jetsAK4CHSTight_Pt"); nbins.push_back(25); Min.push_back(0); Max.push_back(1500); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK4CHS_Eta[0]");name.push_back("jetsAK4CHS_Eta"); nbins.push_back(24); Min.push_back(-6); Max.push_back(6); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK4CHSTight_Eta[0]");name.push_back("jetsAK4CHSTight_Eta"); nbins.push_back(50); Min.push_back(-5); Max.push_back(5); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK4CHS_Phi[0]");name.push_back("jetsAK4CHS_Phi"); nbins.push_back(40); Min.push_back(-4); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK4CHSTight_Phi[0]");name.push_back("jetsAK4CHSTight_Phi"); nbins.push_back(40); Min.push_back(-4); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);
  
  var.push_back("jetsAK8CHS_size");name.push_back("jetsAK8CHS_size"); nbins.push_back(35); Min.push_back(0); Max.push_back(7); sel.push_back(""); presel.push_back(pre_sel); 
  var.push_back("jetsAK8CHS_E[0]");name.push_back("jetsAK8CHS_E"); nbins.push_back(40); Min.push_back(0); Max.push_back(4000); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK8CHS_Pt[0]");name.push_back("jetsAK8CHS_Pt"); nbins.push_back(24); Min.push_back(0); Max.push_back(1200); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK8CHS_Eta[0]");name.push_back("jetsAK8CHS_Eta"); nbins.push_back(24); Min.push_back(-4); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK8CHS_Phi[0]");name.push_back("jetsAK8CHS_Phi"); nbins.push_back(40); Min.push_back(-6); Max.push_back(6); sel.push_back(""); presel.push_back(pre_sel);
    
  var.push_back("jetsAK8CHS_Eta[0]-jetsAK4CHSTight_Eta[0]");name.push_back("deltaEta"); nbins.push_back(40); Min.push_back(-10); Max.push_back(10); sel.push_back("(sqrt(pow(jetsAK8CHS_Phi[0]-jetsAK4CHSTight_Phi[0],2)+pow(jetsAK8CHS_Eta[0]-jetsAK4CHSTight_Eta[0],2))>1.2 && jetsAK4CHSTight_IsCSVM[0]>0)"); presel.push_back((" && "+ pre_sel));
  var.push_back("jetsAK8CHS_Phi[0]-jetsAK4CHSTight_Phi[0]");name.push_back("deltaPhi"); nbins.push_back(40); Min.push_back(-10); Max.push_back(10); sel.push_back("sqrt(pow(jetsAK8CHS_Phi[0]-jetsAK4CHSTight_Phi[0],2)+pow(jetsAK8CHS_Eta[0]-jetsAK4CHSTight_Eta[0],2))>1.2  && jetsAK4CHSTight_IsCSVM[0]>0"); presel.push_back((" && "+ pre_sel));

  var.push_back("sqrt(pow(jetsAK8CHS_E[0]+jetsAK4CHSTight_E[0],2)-(pow((jetsAK8CHS_Pt[0]*cos(jetsAK8CHS_Phi[0])+jetsAK4CHSTight_Pt[0]*cos(jetsAK4CHSTight_Phi[0])),2)-pow((jetsAK8CHS_Pt[0]*sin(jetsAK8CHS_Phi[0])+jetsAK4CHSTight_Pt[0]*sin(jetsAK4CHSTight_Phi[0])),2)-pow((jetsAK8CHS_Pt[0]*sinh(jetsAK8CHS_Eta[0])+jetsAK4CHSTight_Pt[0]*sinh(jetsAK4CHSTight_Eta[0])),2)))");name.push_back("massa_invariante"); nbins.push_back(50); Min.push_back(0); Max.push_back(4000); sel.push_back("(sqrt(pow(jetsAK8CHS_Phi[0]-jetsAK4CHSTight_Phi[0],2)+pow(jetsAK8CHS_Eta[0]-jetsAK4CHSTight_Eta[0],2))>1.2 && jetsAK4CHSTight_IsCSVM[0]>0)"); presel.push_back(" && "+pre_sel);

  //var.push_back("jetsAK8Puppi_size");name.push_back("jetsAK8Puppi_size"); nbins.push_back(35); Min.push_back(0); Max.push_back(7); sel.push_back("");
  //var.push_back("jetsAK8Puppi_E");name.push_back("jetsAK8Puppi_E"); nbins.push_back(30); Min.push_back(0); Max.push_back(3000); sel.push_back("");
  //var.push_back("jetsAK8Puppi_Pt[0]");name.push_back("jetsAK8Puppi_Pt"); nbins.push_back(30); Min.push_back(0); Max.push_back(1200); sel.push_back("");
  //var.push_back("jetsAK8Puppi_Eta[0]");name.push_back("jetsAK8Puppi_Eta"); nbins.push_back(30); Min.push_back(-3); Max.push_back(3); sel.push_back("");
  //var.push_back("jetsAK8Puppi_Phi[0]");name.push_back("jetsAK8Puppi_Phi"); nbins.push_back(40); Min.push_back(-4); Max.push_back(4); sel.push_back("");
  
  var.push_back("jetsAK8CHS_tau1CHS[0]");name.push_back("jetsAK8CHS_tau1CHS"); nbins.push_back(50); Min.push_back(0); Max.push_back(1); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK8CHS_tau2CHS[0]");name.push_back("jetsAK8CHS_tau2CHS"); nbins.push_back(50); Min.push_back(0); Max.push_back(1); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK8CHS_tau3CHS[0]");name.push_back("jetsAK8CHS_tau3CHS"); nbins.push_back(50); Min.push_back(0); Max.push_back(1); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK8CHS_tau3CHS[0]/jetsAK8CHS_tau2CHS[0]");name.push_back("tau3_2"); nbins.push_back(50); Min.push_back(0); Max.push_back(1); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("jetsAK8CHS_tau2CHS[0]/jetsAK8CHS_tau1CHS[0]");name.push_back("tau2_1"); nbins.push_back(50); Min.push_back(0); Max.push_back(1); sel.push_back(""); presel.push_back(pre_sel);

  var.push_back("subjetsAK8CHS_size");name.push_back("subjetsAK8CHS_size"); nbins.push_back(20); Min.push_back(0); Max.push_back(10); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("subjetsAK8CHS_E[0]");name.push_back("subjetsAK8CHS_E"); nbins.push_back(20); Min.push_back(0); Max.push_back(2000); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("subjetsAK8CHS_Pt[0]");name.push_back("subjetsAK8CHS_Pt"); nbins.push_back(20); Min.push_back(0); Max.push_back(1000); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("subjetsAK8CHS_Eta[0]");name.push_back("subjetsAK8CHS_Eta"); nbins.push_back(40); Min.push_back(-4); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("subjetsAK8CHS_Phi[0]");name.push_back("subjetsAK8CHS_Phi"); nbins.push_back(40); Min.push_back(-4); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);

  var.push_back("subjetsAK8Puppi_size");name.push_back("subjetsAK8Puppi_size"); nbins.push_back(20); Min.push_back(0); Max.push_back(10); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("subjetsAK8Puppi_E[0]");name.push_back("subjetsAK8Puppi_E"); nbins.push_back(20); Min.push_back(0); Max.push_back(2000); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("subjetsAK8Puppi_Pt[0]");name.push_back("subjetsAK8Puppi_Pt"); nbins.push_back(20); Min.push_back(0); Max.push_back(1000); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("subjetsAK8Puppi_Eta[0]");name.push_back("subjetsAK8Puppi_Eta"); nbins.push_back(40); Min.push_back(-4); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("subjetsAK8Puppi_Phi[0]");name.push_back("subjetsAK8Puppi_Phi"); nbins.push_back(40); Min.push_back(-4); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);

  var.push_back("metFull_size");name.push_back("metFull_size"); nbins.push_back(30); Min.push_back(0); Max.push_back(3); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("metFull_E[0]");name.push_back("metFull_E"); nbins.push_back(20); Min.push_back(-1); Max.push_back(1); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("metFull_Pt[0]");name.push_back("metFull_Pt"); nbins.push_back(16); Min.push_back(0); Max.push_back(800); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("metFull_Px[0]");name.push_back("metFull_Px"); nbins.push_back(32); Min.push_back(-800); Max.push_back(800); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("metFull_Py[0]");name.push_back("metFull_Py"); nbins.push_back(32); Min.push_back(-800); Max.push_back(800); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("metFull_Eta[0]");name.push_back("metFull_Eta"); nbins.push_back(20); Min.push_back(-1); Max.push_back(1); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("metFull_Phi[0]");name.push_back("metFull_Phi"); nbins.push_back(40); Min.push_back(-4); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("metFull_uncorPt[0]");name.push_back("metFull_uncorPt"); nbins.push_back(16); Min.push_back(0); Max.push_back(800); sel.push_back(""); presel.push_back(pre_sel);
  var.push_back("metFull_uncorPhi[0]");name.push_back("metFull_uncorPhi"); nbins.push_back(20); Min.push_back(-4); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);
  
  var.push_back("muonsLoose_size");name.push_back("muonsLoose_size"); nbins.push_back(40); Min.push_back(0); Max.push_back(4); sel.push_back(""); presel.push_back(pre_sel);

  
  TFile * output = TFile::Open(name_file_histo.c_str(),"RECREATE");
  //TFile * output = TFile::Open(name_file_histo.c_str(),"UPDATE");
  output->Close();
  
  for (int i=0; i<var.size();++i){
    
    TH1D a ("a", name[i], nbins[i],Min[i],Max[i]);
    if(l==0) fChain->Project("a",var[i],sel[i]);
    else fChain->Project("a",var[i],sel[i]+presel[i]);
    output = TFile::Open(name_file_histo.c_str(),"UPDATE");
    a.Write(name[i]);
    a.Scale(1/a.Integral());
    a.Write(name[i]+"_normalizzato");
    output->Close();
    
  }
  
}
