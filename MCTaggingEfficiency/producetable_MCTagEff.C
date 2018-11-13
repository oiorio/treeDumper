#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TFile.h"
#include "TStyle.h"
#include "TDirectory.h"
#include "TString.h"
#include "TGaxis.h"
#include "TLine.h"
#include <vector>
#include <iostream>
using namespace std;

void h12ascii (TH2D* h, ofstream &myfile ){

  Int_t nx = h->GetNbinsX();
  Int_t ny = h->GetNbinsY();
  
  for (Int_t i=1; i<=nx; i++) {
    for (Int_t j=1; j<=ny; j++) {
      myfile << "\t\tif(((eta>=" <<  h->GetYaxis()->GetBinLowEdge(i) << ") && (eta <" << h->GetYaxis()->GetBinUpEdge(i) << ")) && ((pt >=  " << h->GetXaxis()->GetBinLowEdge(j) << ") && (pt < " << h->GetXaxis()->GetBinUpEdge(j)  << "))) return "<< h->GetBinContent(i,j) << "; \n"; 
    }
  }

}

void producetable_MCTagEff()
{

  ofstream myfile;
  myfile.open ("MCTagEff_AK8SoftDropSubj.txt", std::ofstream::app);
  cout << "started creating file txt..." << endl;
  
  TFile *file_t = new TFile("Maps/bTaggingEfficiency_AK8SoftDropSubj_CSVv2t_AK8Subj_CSVT_bTaggingEfficiencyMap.root");
  TFile *file_m = new TFile("Maps/bTaggingEfficiency_AK8SoftDropSubj_CSVv2m_AK8Subj_CSVM_bTaggingEfficiencyMap.root");
  TFile *file_l = new TFile("Maps/bTaggingEfficiency_AK8SoftDropSubj_CSVv2l_AK8Subj_CSVL_bTaggingEfficiencyMap.root");
  
  cout << "opened files..." << endl;

  myfile << "if ((abs(flavor) ==5)){ " << endl;
  myfile << "\t if (algo==\"csvt\") {" << endl;
  h12ascii((TH2D*)file_t->Get("efficiency_b"),myfile);
  myfile << "\t }" << endl;
  myfile << "\t if (algo==\"csvm\") {" << endl;
  h12ascii((TH2D*)file_m->Get("efficiency_b"),myfile);
  myfile << "\t }" << endl;
  myfile << "\t if (algo==\"csvl\") {" << endl;
  h12ascii((TH2D*)file_l->Get("efficiency_b"),myfile);
  myfile << "\t }" << endl;
  myfile << "}" << endl;
  cout << "flavor b done!" << endl;

  myfile << "if ((abs(flavor) ==4)){ " << endl;
  myfile << "\t if (algo==\"csvt\") {" << endl;
  h12ascii((TH2D*)file_t->Get("efficiency_c"),myfile);
  myfile << "\t }" << endl;
  myfile << "\t if (algo==\"csvm\") {" << endl;
  h12ascii((TH2D*)file_m->Get("efficiency_c"),myfile);
  myfile << "\t}" << endl;
  myfile << "\t if (algo==\"csvl\") {" << endl; 
  h12ascii((TH2D*)file_l->Get("efficiency_c"),myfile);
  myfile << "\t}" << endl;
  myfile << "}" << endl;
  cout << "flavor c done!" << endl;

  myfile << "if ((abs(flavor)!= 5 && (abs(flavor)!= 4 )){ " << endl;
  myfile << "\t if (algo==\"csvt\") {" << endl;
  h12ascii((TH2D*)file_t->Get("efficiency_udsg"),myfile);
  myfile << "\t }" << endl;
  myfile << "\t if (algo==\"csvm\") {" << endl;
  h12ascii((TH2D*)file_m->Get("efficiency_udsg"),myfile);
  myfile << "\t }" << endl;
  myfile << "\t if (algo==\"csvl\") {" << endl; 
  h12ascii((TH2D*)file_l->Get("efficiency_udsg"),myfile);
  myfile << "\t }" << endl;
  myfile << "}" << endl;
  cout << "flavor light done!" << endl;

  myfile.close();
}
