void compare_shapes_DR() {
 
  // Braden Allmond -- Nov 7th, 2023
  //
  // comparing DRsr and DRar shapes in one plot
  //
  // for i in 0j 1j 2j
  //   get new canvas to draw on
  //   for graph in 1 2 3 final
  //     open file, get hist, draw it, close file
  //   save canvas
  
  string directory = "inputs_2p5";
  cout << " is " << directory << " the directory you want to be reading?" << endl;
  string channel = "ditau";
  cout << " is " << channel << " the channel you want to be using?" << endl;
  map<string, vector<string>> channelGraphs{
                                {"ditau", {"mvis", "tau1pt", "tau2pt"}}
                              };
  
  vector<string> possibleNumJets = {"0j", "1j", "2j"};
  vector<string> possibleGraphs  = channelGraphs[channel];

  for(string numJets : possibleNumJets){
    for (string graphVar : possibleGraphs){
      TCanvas *can = new TCanvas(("Compare DRsr and DRar for "+numJets).c_str(), ("DR shape and ratio for "+numJets).c_str(), 600, 600);
      can->cd();
      can->SetLeftMargin(0.14); // always do this so your axis tick marks aren't clipped
      can->SetRightMargin(0.06);
      TLegend *legend = new TLegend(0.6, 0.6, 0.94, 0.9);
      legend->SetFillStyle(0);

      string data_DRar_var = "data_"+graphVar+"_DRar_"+numJets;
      TFile* DRar_file = TFile::Open((directory+"/"+data_DRar_var+"_1.root").c_str(), "READ");
      TH1F* DRar_hist = (TH1F*)DRar_file->Get((data_DRar_var).c_str());
      DRar_hist->SetDirectory(0);
      //DRar_hist->SetMarkerStyle(3);
      DRar_hist->SetStats(false);
      stringstream DRar_integ;
      DRar_integ << fixed << setprecision(2) << DRar_hist->Integral();
      legend->AddEntry(DRar_hist,("DRar ("+DRar_integ.str()+")").c_str());
      //legend->AddEntry(DRar_hist,"DRar");
      DRar_hist->Draw();
      DRar_hist->GetXaxis()->SetTitle(graphVar.c_str());
      DRar_hist->GetYaxis()->SetTitle("Events");
      DRar_file->Close(); 

      cout << "working on file for " << numJets << " jets and " << graphVar << " graph" << endl;
      // "bkg_tau1pt_DRar_0j_1.root"
      // should be all bkgs, and data
      // but just do data for now...
      // "data_tau1pt_DRar_0j_1.root"
      string data_DRsr_var = "data_"+graphVar+"_DRsr_"+numJets;
      TFile* DRsr_file = TFile::Open((directory+"/"+data_DRsr_var+"_1.root").c_str(), "READ");
      TH1F* DRsr_hist = (TH1F*)DRsr_file->Get((data_DRsr_var).c_str());
      DRsr_hist->SetDirectory(0);
      DRsr_hist->SetLineColor(2);
      //DRsr_hist->SetMarkerStyle(14);
      DRsr_hist->SetStats(false);
      stringstream DRsr_integ;
      DRsr_integ << fixed << setprecision(2) << DRsr_hist->Integral();
      legend->AddEntry(DRsr_hist,("DRsr ("+DRsr_integ.str()+")").c_str());
      //legend->AddEntry(DRsr_hist,"DRsr");
      DRsr_hist->Draw("SAME");
      DRsr_file->Close(); 

      legend->Draw();
      can->SaveAs(("compare_DRs_"+graphVar+"_"+numJets+".png").c_str());
      delete can;
    }
  }
}
