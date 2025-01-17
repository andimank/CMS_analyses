void plot_ratio_effcor() {

  const int nBin=12;
  double n[nBin];

  double vn_Cfit_reco[nBin];
  double vne_Cfit_reco[nBin];
  double vn_Cfit_recocor[nBin];
  double vne_Cfit_recocor[nBin];
  double tmp;

  double vn_avgcos_reco[nBin];
  double vn_avgcos_recocor[nBin];

  // --- Set Centrality --- //
  // -----------------------//
  int centmin = 5;          //
  int centmax = 10;         //
  // -----------------------//

  // Read in the vn results from correlation function fit files

  // no eff. cor.
  ifstream in_vn_Cfit_reco;
  in_vn_Cfit_reco.open(Form("merged_all_pt_vn_%d_%dcent_2k_noeffcor.txt", centmin, centmax));
  if(!in_vn_Cfit_reco.good())    cout << " Corr. function input fail - no eff. cor." << endl;
  for(int i=0;i<nBin;i++){
    in_vn_Cfit_reco>>tmp;
    in_vn_Cfit_reco>>tmp;
    //in_vn_Cfit_reco>>n[i];
    in_vn_Cfit_reco>>tmp;
    in_vn_Cfit_reco>>tmp;
    in_vn_Cfit_reco>>tmp;
    in_vn_Cfit_reco>>vn_Cfit_reco[i];
    in_vn_Cfit_reco>>vne_Cfit_reco[i];
  }

  // eff. cor.
  ifstream in_vn_Cfit_recocor;
  in_vn_Cfit_recocor.open(Form("merged_all_pt_vn_%d_%dcent_2k.txt", centmin, centmax));
  if(!in_vn_Cfit_recocor.good())    cout << " Corr. function input fail - eff. cor" << endl;
  for(int i=0;i<nBin;i++){
    in_vn_Cfit_recocor>>tmp;
    in_vn_Cfit_recocor>>tmp;
    in_vn_Cfit_recocor>>n[i];
    in_vn_Cfit_recocor>>tmp;
    in_vn_Cfit_recocor>>tmp;
    in_vn_Cfit_recocor>>vn_Cfit_recocor[i];
    in_vn_Cfit_recocor>>vne_Cfit_recocor[i];
  }

  // change the type casting for correlation function results for easier plotting
  for (int i = 0; i < nBin; i++) {
      vn_Cfit_reco[i] = static_cast<float_t>(vn_Cfit_reco[i]);
      vn_Cfit_recocor[i] = static_cast<float_t>(vn_Cfit_recocor[i]);
  }



  ////////// Analyze the avgcos method output root file to get the vn

  TFile*f = new TFile ("../test/interactive/higherHarmonics_PbPb_2pc_final_test_run_2k.root"); //higherHarmonics_PbPb_2pc_test_run_3k_effcor.root

  TDirectory* dirG1d= (TDirectory*)f->Get(Form("defaultAnalysis_%d%d",centmin, centmax));

  TDirectory* dir1= (TDirectory*)dirG1d->Get("Signal");
  TDirectory* dir2= (TDirectory*)dirG1d->Get("Background");

  // no eff. cor.
  TH1D* hSignal_all_reco_evtavgcos_n_[nBin];
  TH1D* hBackground_all_reco_evtavgcos_n_[nBin];
  TH1D* hBackground_all_reco_evtavgcos_2n_[nBin];

  // eff. cor
  TH1D* hSignal_all_recocor_evtavgcos_n_[nBin];
  TH1D* hBackground_all_recocor_evtavgcos_n_[nBin];
  TH1D* hBackground_all_recocor_evtavgcos_2n_[nBin];

  // Read (individual harmonics) histograms from file
  for(unsigned int n_i = 2; n_i<nBin+1; n_i++){

    // no eff. cor.
    hSignal_all_reco_evtavgcos_n_[n_i] = (TH1D*)dir1->Get(Form("signal_all_reco_evtavgcos_n_%d", n_i));
    //hBackground_all_reco_evtavgcos_n_[n_i] = (TH1D*)dir2->Get(Form("background_all_reco_evtavgcos_n_%d", n_i));
    //hBackground_all_reco_evtavgcos_2n_[n_i] = (TH1D*)dir2->Get(Form("background_all_reco_evtavgcos_2n_%d", n_i));
    hBackground_all_reco_evtavgcos_n_[n_i] = (TH1D*)f->Get(Form("defaultAnalysis_510/Background/background_all_reco_evtavgcos_n_%d", n_i));
    hBackground_all_reco_evtavgcos_2n_[n_i] = (TH1D*)f->Get(Form("defaultAnalysis_510/Background/background_all_reco_evtavgcos_2n_%d", n_i));

    // eff. cor
    hSignal_all_recocor_evtavgcos_n_[n_i] = (TH1D*)dir1->Get(Form("signal_all_recocor_evtavgcos_n_%d", n_i));
    //hBackground_all_recocor_evtavgcos_n_[n_i] = (TH1D*)dir2->Get(Form("background_all_recocor_evtavgcos_n_%d", n_i));
    //hBackground_all_recocor_evtavgcos_2n_[n_i] = (TH1D*)dir2->Get(Form("background_all_recocor_evtavgcos_2n_%d", n_i));
    hBackground_all_recocor_evtavgcos_n_[n_i] = (TH1D*)f->Get(Form("defaultAnalysis_510/Background/background_all_recocor_evtavgcos_n_%d", n_i));
    hBackground_all_recocor_evtavgcos_2n_[n_i] = (TH1D*)f->Get(Form("defaultAnalysis_510/Background/background_all_recocor_evtavgcos_2n_%d", n_i));

  }


 
  // Claculate vn
  for(unsigned int n_i = 2; n_i<nBin+1; n_i++){

    // no eff. cor.
    Float_t mean_signal_reco = hSignal_all_reco_evtavgcos_n_[n_i]->GetMean();
    Float_t mean_background_reco = hBackground_all_reco_evtavgcos_n_[n_i]->GetMean();
    Float_t mean_denom_reco = hBackground_all_reco_evtavgcos_2n_[n_i]->GetMean();

    Float_t vn_delta_reco = (mean_signal_reco - mean_background_reco)/(1 + mean_denom_reco);
    Float_t vn_reco_tmp = sqrt(std::abs(vn_delta_reco));
    vn_avgcos_reco[n_i-1] = vn_reco_tmp;

    // eff. cor.
    Float_t mean_signal_recocor = hSignal_all_recocor_evtavgcos_n_[n_i]->GetMean();
    Float_t mean_background_recocor = hBackground_all_recocor_evtavgcos_n_[n_i]->GetMean();
    Float_t mean_denom_recocor = hBackground_all_recocor_evtavgcos_2n_[n_i]->GetMean();

    Float_t vn_delta_recocor = (mean_signal_recocor - mean_background_recocor)/(1 + mean_denom_recocor);
    Float_t vn_recocor_tmp = sqrt(std::abs(vn_delta_recocor));
    vn_avgcos_recocor[n_i-1] = vn_recocor_tmp;


  }


  // Create canvas and configure
  TCanvas *c1 = new TCanvas("c1","c_ratios",700,900);

  // Create the two panels and configure
  auto *panel2 = new TPad("panel2","panel2",0.,0.,1.,0.35); panel2->Draw();
  panel2->SetLeftMargin(0.12);
  panel2->SetTopMargin(0.001);
  panel2->SetBottomMargin(0.3);
  panel2->SetGrid();
  panel2->SetTicks();
  auto *panel1 = new TPad("panel1","panel1",0.,0.35,1.,1.);  panel1->Draw();
  panel1->SetLeftMargin(0.12);
  panel1->SetBottomMargin(0.001);
  panel1->cd();
  panel1->SetGrid();
  panel1->SetTicks();
  panel1->SetLogy();

  // Create top-panel graphs and configure
  TGraph*graph_recocor = new TGraph(nBin,n,vn_Cfit_recocor);
  TGraph*graph_reco = new TGraph(nBin,n,vn_Cfit_reco);
  TGraph*graph_avgcos_recocor = new TGraph(nBin,n,vn_avgcos_recocor);
  TGraph*graph_avgcos_reco = new TGraph(nBin,n,vn_avgcos_reco);

  graph_recocor->SetTitle("");
  graph_recocor->SetMarkerColor(1);
  graph_recocor->SetMarkerStyle(20); // full circle
  graph_recocor->SetMarkerSize(1.3);
  graph_recocor->GetYaxis()->SetTitle("v_{n}{2, |#Delta#eta|>2}");
  graph_recocor->GetYaxis()->CenterTitle();
  graph_recocor->GetYaxis()->SetTitleOffset(1.0);
  graph_recocor->GetYaxis()->SetTitleSize(0.05);
  graph_recocor->GetXaxis()->SetRangeUser(1.5,13);
  graph_recocor->GetYaxis()->SetRangeUser(0.00004,0.95);
  graph_recocor->GetXaxis()->SetTickSize(0.);
  graph_recocor->Draw("AP");

  graph_reco->SetMarkerColor(1);
  graph_reco->SetMarkerStyle(24); // open circle
  graph_reco->SetMarkerSize(1.3);
  graph_reco->GetXaxis()->SetRangeUser(1.5,13);
  graph_reco->Draw("P");

  graph_avgcos_recocor->SetMarkerColor(2);
  graph_avgcos_recocor->SetMarkerStyle(20); // full circle
  graph_avgcos_recocor->SetMarkerSize(1.3);
  graph_avgcos_recocor->SetLineColor(2);
  graph_avgcos_recocor->GetXaxis()->SetRangeUser(1.5,13);
  graph_avgcos_recocor->Draw("P");

  graph_avgcos_reco->SetMarkerColor(2);
  graph_avgcos_reco->SetMarkerStyle(24); // open circle
  graph_avgcos_reco->SetMarkerSize(1.3);
  graph_avgcos_reco->SetLineColor(2);
  graph_avgcos_reco->GetXaxis()->SetRangeUser(1.5,13);
  graph_avgcos_reco->Draw("P");

  // Create legend and configure
  TLegend *leg = new TLegend(0.14,0.05,0.48,0.25);
  leg->AddEntry(graph_recocor,"v_{n} C-fit, (eff. cor.)","lp");
  leg->AddEntry(graph_reco,"v_{n} C-fit, (no eff. cor.)","lp");
  leg->AddEntry(graph_avgcos_recocor,"v_{n} avg-cos, (eff. cor.)","lp");
  leg->AddEntry(graph_avgcos_reco,"v_{n} avg-cos, (no eff. cor.)","lp");
  leg->SetTextSize(0.03);
  leg->Draw();

  // Create text annotations and configure
  TLatex *tex1 = new TLatex(0,0,"CMS");
  tex1->SetTextColor(1);
  tex1->SetTextSize(0.065);
  tex1->SetTextFont(42);
  //tex1->Draw();
  tex1->DrawLatexNDC(0.165,0.80,"CMS                  PbPb 5.02 TeV");

  TLatex *tex2 = new TLatex(0,0,Form("%d-%d%s cent.", centmin, centmax, "%"));
  tex2->SetTextColor(1);
  tex2->SetTextSize(0.048);
  tex2->SetTextFont(42);
  //tex2->Draw();
  tex2->DrawLatexNDC(0.58,0.16,Form("%d-%d%s cent.", centmin, centmax, "%"));

  TLatex *tex3 = new TLatex(0,0,"|#eta| < 2.4");
  tex3->SetTextColor(1);
  tex3->SetTextSize(0.045);
  tex3->SetTextFont(42);
  //tex3->Draw();
  tex3->DrawLatexNDC(0.58,0.10,"|#eta| < 2.4");

  TLatex *tex4 = new TLatex(0,0,"0.5 < p_{T} < 3.0 GeV");
  tex4->SetTextColor(1);
  tex4->SetTextSize(0.045);
  tex4->SetTextFont(42);
  //tex4->Draw();
  tex4->DrawLatexNDC(0.58,0.04,"0.5 < p_{T} < 3.0 GeV");


  // Create bottom-panel graph and configure
  panel2->cd();
  TGraph*ratio = new TGraph(nBin); ratio->SetTitle("");
  for (int i=1; i<nBin; i++) ratio->SetPoint(i, n[i], vn_Cfit_recocor[i]/vn_Cfit_reco[i]);
  TGraph*ratio_avgcos = new TGraph(nBin); ratio_avgcos->SetTitle("");
  for (int i=1; i<nBin; i++) ratio_avgcos->SetPoint(i, n[i], vn_avgcos_recocor[i]/vn_avgcos_reco[i]);
  // for (int i=1; i<nBin; i++){
  //  std::cout << "   " << i << "   " << n[i] << "   " << vn_avgcos_recocor[i]/vn_avgcos_reco[i] << std::endl;
  // }

  ratio->SetMarkerSize(1.5);
  ratio->SetMarkerColor(1);

  ratio_avgcos->SetMarkerSize(1.5);
  ratio_avgcos->SetMarkerColor(2);

  TMultiGraph *combined_ratio = new TMultiGraph();
  combined_ratio->Add(ratio);
  combined_ratio->Add(ratio_avgcos);
  combined_ratio->GetYaxis()->SetTitle("eff. cor. / no eff. cor.");
  combined_ratio->GetYaxis()->CenterTitle();
  combined_ratio->GetYaxis()->SetTitleOffset(0.7);
  combined_ratio->GetYaxis()->SetTitleSize(0.06);
  combined_ratio->GetXaxis()->SetTitle("n");
  combined_ratio->GetXaxis()->CenterTitle();
  combined_ratio->GetXaxis()->SetTitleOffset(1.0);
  combined_ratio->GetXaxis()->SetTitleSize(0.1);
  combined_ratio->GetXaxis()->SetLabelSize(0.075);
  combined_ratio->GetYaxis()->SetLabelSize(0.075);
  combined_ratio->GetXaxis()->SetLimits(1.44, 13.1);
  combined_ratio->GetYaxis()->SetRangeUser(-0.5, 2.9);

  // Create a horizontal line
  TGraph* hline = new TGraph(2);
  hline->SetPoint(0, 0, 1);
  hline->SetPoint(1, 14, 1);
  hline->SetLineColor(34);
  combined_ratio->Add(hline);
  // Get the individual graphs in the TMultiGraph and set the line width to 0 (leaving the g_i=2 hline utouched)
  TList* graphs = combined_ratio->GetListOfGraphs();
  for (int g_i = 0; g_i < graphs->GetSize(); g_i++) {
      TGraph* g = (TGraph*)graphs->At(g_i);
      if (g_i != 2) {
          g->SetLineWidth(0);
      }
  }

  combined_ratio->Draw("AL*");

  // save figures
  c1->Print(Form("./plot/plot_vn_%d_%dcen_effcor_ratios_2k_events.png", centmin, centmax));
  c1->Print(Form("./plot/plot_vn_%d_%dcen_effcor_ratios_2k_events.pdf", centmin, centmax));

}
