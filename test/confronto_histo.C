{  
  TFile *f_fondo = TFile::Open("histo_fondo.root"); 
  TFile *f_segnale = TFile::Open("histo_segnale.root");
  TFile *f_QCD = TFile::Open("histo_QCD.root");

  TFile * output = TFile::Open("confronto_histo.root","RECREATE");
  //TFile * output = TFile::Open("confronto_histo.root","UPDATE");   
  output->Close();

  std::vector<TString> name_var,name_histo;
  
  name_var.push_back("jetsAK4CHS_size");name_histo.push_back("jetsAK4CHS_size_normalizzato");
  name_var.push_back("jetsAK4CHSTight_size");name_histo.push_back("jetsAK4CHSTight_size_normalizzato");
  name_var.push_back("jetsAK4CHS_E");name_histo.push_back("jetsAK4CHS_E_normalizzato");
  name_var.push_back("jetsAK4CHSTight_E");name_histo.push_back("jetsAK4CHSTight_E_normalizzato");
  name_var.push_back("jetsAK4CHS_Pt");name_histo.push_back("jetsAK4CHS_Pt_normalizzato");
  name_var.push_back("jetsAK4CHSTight_Pt");name_histo.push_back("jetsAK4CHSTight_Pt_normalizzato"); 
  name_var.push_back("jetsAK4CHS_Eta");name_histo.push_back("jetsAK4CHS_Eta_normalizzato");
  name_var.push_back("jetsAK4CHSTight_Eta");name_histo.push_back("jetsAK4CHSTight_Eta_normalizzato");
  name_var.push_back("jetsAK4CHS_Phi");name_histo.push_back("jetsAK4CHS_Phi_normalizzato"); 
  name_var.push_back("jetsAK4CHSTight_Phi");name_histo.push_back("jetsAK4CHSTight_Phi_normalizzato");
  
  name_var.push_back("jetsAK8CHS_size");name_histo.push_back("jetsAK8CHS_size_normalizzato");
  name_var.push_back("jetsAK8CHS_E");name_histo.push_back("jetsAK8CHS_E_normalizzato");
  name_var.push_back("jetsAK8CHS_Pt");name_histo.push_back("jetsAK8CHS_Pt_normalizzato");
  name_var.push_back("jetsAK8CHS_Eta");name_histo.push_back("jetsAK8CHS_Eta_normalizzato"); 
  name_var.push_back("jetsAK8CHS_Phi");name_histo.push_back("jetsAK8CHS_Phi_normalizzato"); 

  name_var.push_back("deltaEta");name_histo.push_back("deltaEta_normalizzato"); 
  name_var.push_back("deltaPhi");name_histo.push_back("deltaPhi_normalizzato");
  
  name_var.push_back("massa_invariante");name_histo.push_back("massa_invariante_normalizzato");
  
  //name_var.push_back("jetsAK8Puppi_size");name_histo.push_back("jetsAK8Puppi_size_normalizzato");
  //name_var.push_back("jetsAK8Puppi_E");name_histo.push_back("jetsAK8Puppi_E_normalizzato"); 
  //name_var.push_back("jetsAK8Puppi_Pt");name_histo.push_back("jetsAK8Puppi_Pt_normalizzato"); 
  //name_var.push_back("jetsAK8Puppi_Eta");name_histo.push_back("jetsAK8Puppi_Eta_normalizzato"); 
  //name_var.push_back("jetsAK8Puppi_Phi");name_histo.push_back("jetsAK8Puppi_Phi_normalizzato");
  
  name_var.push_back("jetsAK8CHS_tau1CHS");name_histo.push_back("jetsAK8CHS_tau1CHS_normalizzato");
  name_var.push_back("jetsAK8CHS_tau2CHS");name_histo.push_back("jetsAK8CHS_tau2CHS_normalizzato");
  name_var.push_back("jetsAK8CHS_tau3CHS");name_histo.push_back("jetsAK8CHS_tau3CHS_normalizzato");
  name_var.push_back("tau3_2");name_histo.push_back("tau3_2_normalizzato");
  name_var.push_back("tau2_1");name_histo.push_back("tau2_1_normalizzato");
  
  name_var.push_back("subjetsAK8CHS_size");name_histo.push_back("subjetsAK8CHS_size_normalizzato");
  name_var.push_back("subjetsAK8CHS_E");name_histo.push_back("subjetsAK8CHS_E_normalizzato"); 
  name_var.push_back("subjetsAK8CHS_Pt");name_histo.push_back("subjetsAK8CHS_Pt_normalizzato"); 
  name_var.push_back("subjetsAK8CHS_Eta");name_histo.push_back("subjetsAK8CHS_Eta_normalizzato");
  name_var.push_back("subjetsAK8CHS_Phi");name_histo.push_back("subjetsAK8CHS_Phi_normalizzato");

  name_var.push_back("subjetsAK8Puppi_size");name_histo.push_back("subjetsAK8Puppi_size_normalizzato"); 
  name_var.push_back("subjetsAK8Puppi_E");name_histo.push_back("subjetsAK8Puppi_E_normalizzato"); 
  name_var.push_back("subjetsAK8Puppi_Pt");name_histo.push_back("subjetsAK8Puppi_Pt_normalizzato");
  name_var.push_back("subjetsAK8Puppi_Eta");name_histo.push_back("subjetsAK8Puppi_Eta_normalizzato");
  name_var.push_back("subjetsAK8Puppi_Phi");name_histo.push_back("subjetsAK8Puppi_Phi_normalizzato");

  name_var.push_back("metFull_size");name_histo.push_back("metFull_size_normalizzato");
  name_var.push_back("metFull_E");name_histo.push_back("metFull_E_normalizzato");
  name_var.push_back("metFull_Pt");name_histo.push_back("metFull_Pt_normalizzato");
  name_var.push_back("metFull_Px");name_histo.push_back("metFull_Px_normalizzato");
  name_var.push_back("metFull_Py");name_histo.push_back("metFull_Py_normalizzato");
  name_var.push_back("metFull_Eta");name_histo.push_back("metFull_Eta_normalizzato");
  name_var.push_back("metFull_Phi");name_histo.push_back("metFull_Phi_normalizzato");
  name_var.push_back("metFull_uncorPt");name_histo.push_back("metFull_uncorPt_normalizzato");
  name_var.push_back("metFull_uncorPhi");name_histo.push_back("metFull_uncorPhi_normalizzato");
  
  name_var.push_back("muonsLoose_size");name_histo.push_back("muonsLoose_size_normalizzato");  
  
  for (int i=0; i<name_histo.size();++i){

    TFile * output = TFile::Open("confronto_histo.root","UPDATE"); 
    TCanvas *c1 = new TCanvas();
    
    TH1D *h_segnale = (TH1D*)f_segnale->Get(name_histo[i]);
    double max_sig = h_segnale->GetMaximum();   
    TH1D *h_fondo = (TH1D*)f_fondo->Get(name_histo[i]);
    double max_bkg = h_fondo->GetMaximum();
    TH1D *h_qcd = (TH1D*)f_QCD->Get(name_histo[i]);
    double max_qcd = h_qcd->GetMaximum();

    if(max_sig>=max_bkg && max_sig>=max_qcd){
      
      h_segnale->SetTitle(name_var[i]); 
      h_segnale->SetLineColor(kRed);
      h_segnale->Draw();
      h_fondo->SetLineColor(kBlue);
      h_fondo->Draw("same");
      h_qcd->SetLineColor(kGreen);
      h_qcd->Draw("same");

    }

    if(max_bkg>=max_sig && max_bkg>=max_qcd){

      h_fondo->SetTitle(name_var[i]); 
      h_fondo->SetLineColor(kBlue);
      h_fondo->Draw();
      h_segnale->SetLineColor(kRed);
      h_segnale->Draw("same");      
      h_qcd->SetLineColor(kGreen);
      h_qcd->Draw("same");

    }
    
    if(max_qcd>=max_sig && max_qcd>=max_bkg){

      h_qcd->SetTitle(name_var[i]); 
      h_qcd->SetLineColor(kGreen);
      h_qcd->Draw();
      h_segnale->SetLineColor(kRed);
      h_segnale->Draw("same");
      h_fondo->SetLineColor(kBlue);
      h_fondo->Draw("same");

    }
  
    TLegend *legend = new TLegend(0.75,0.65,0.9,0.75);
    legend->AddEntry(h_segnale,"segnale","l");
    legend->AddEntry(h_fondo,"fondo","l");
    legend->AddEntry(h_qcd,"QCD","l");
    legend->Draw();

    c1->Write(name_histo[i]);
    c1->Close();
   
    output->Close();
    
  }

}
