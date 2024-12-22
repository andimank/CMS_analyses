#include <iostream>
#include "TROOT.h"
#include "TObject.h"
#include "TProfile.h"
#include "TFile.h"
#include "TTree.h"
#include "math.h"
#include "TH1.h"

int centmin = 5;
int centmax = 10;

//const int mergepTbin = 12;
const int npTbin = 15;

double phimax = 1.5;
double etamax = 1.0;
double etamin_flow = 2.0;
double etamax_flow = 4.0;
double MaxEtaRange = 4.0;
double flow_cal_min = -TMath::Pi()*0.5;
double flow_cal_max = TMath::Pi()*0.5*3.;

double pTAvg[npTbin+1] = {0.5, 0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0};

void charged_Hadron_flow_Cent()
{
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadBottomMargin(0.10);
  gStyle->SetPadTopMargin(0.10);
  gStyle->SetOptStat(0);  // hide statistics box
  gStyle->SetHistLineWidth(2);
  gStyle->SetLabelFont(132,"xyz");
  gStyle->SetLabelSize(0.04,"xyz"); // size of axis labels (rel. to pad size)
  gStyle->SetTitleFont(132,"xyz");  // roman font
  gStyle->SetTitleSize(0.04,"xyz"); // size of axis titles
  gStyle->SetLabelOffset(0.01,"xyz"); // label offset (to fit bigger font
  gStyle->SetNdivisions(505,"xyz");
  gStyle->SetPalette(1);


  //TFile*file3 = new TFile ("MC_HINtestRun3_0p5pttrg10p0_1p0ptass3p0_eta2p4_cent070_Nov21_2022.root");
  //TFile*file3 = new TFile ("MC_HINtestRun3_0p5pttrg3p0_0p5ptass3p0_eta2p4_cent070_Nov23_2022.root");
  //TFile*file3 = new TFile ("MC_HINtestRun3_PixelTracks_0p3pttrg3p0_0p3ptass3p0_eta2p4_cent070_Nov23_2022.root");
  //TFile*file3 = new TFile ("2022_PbPbData_WoMBTrg_0p5pttrg3p0_0p5ptass3p0_GeneralTrack.root");
  TFile*file3 = new TFile ("./merged_files_combined/merged_cent0510_combined_3k_159_ev.root");
  //TFile*file3 = new TFile ("2022_PbPbData_MBTrg_WALGOMVA_vnCent_0p5pttrg3p0_0p5ptass3p0_GeneralTrack.root");

  TDirectory*dirG1d= (TDirectory*)file3->Get(Form("defaultAnalysis_%d%d",centmin, centmax));
  TDirectory*dir1d= (TDirectory*)dirG1d->Get("Signal");
  TDirectory*dir2d= (TDirectory*)dirG1d->Get("Background");
  TDirectory*dir3d= (TDirectory*)dirG1d->Get("Global");

  TH1D* hNtrig[npTbin];
  TH2D* hSignal[npTbin];
  TH2D* hMixed[npTbin];
  TH2D* hCorr[npTbin];

  TH1D* hNtrig_M[npTbin];
  TH2D* hSignal_M[npTbin];
  TH2D* hMixed_M[npTbin];

  TH1D* hNtrig_merge;
  TH2D* hSignal_merge;
  TH2D* hMixed_merge;
  TH2D* hCorr_merge;


  TH1D* heta_proj[npTbin];
  TH1D* hphi_proj[npTbin];
  TH1D* hphi_proj_flow[npTbin];

  TF1 * fit_dphi_flow[npTbin];

  TH1D* hphi_proj_merge_flow;

  TF1 * fit_dphi_merge_flow;

  TCanvas* c_signal[npTbin];
  TCanvas* c_mixed[npTbin];
  TCanvas* c_corr[npTbin];
  TCanvas* c_deta[npTbin];
  TCanvas* c_dphi[npTbin];
  TCanvas* c_dphi_flow[npTbin];
  TCanvas* c_corr_merge;
  TCanvas* c_dphi_merge_flow;


  hSignal_merge = (TH2D*)dir1d->Get("signal_all_recocor");

  hMixed_merge = (TH2D*)dir2d->Get("background_all_recocor");

  hNtrig_merge = (TH1D*)dir3d->Get("hMult_trg");

  //~~~~~~~~~~~Merge correlation~~~~~~~~~~~~~~~~
  c_corr_merge = new TCanvas("c_corr_merge", "c_corr_merge", 600, 530);
  c_corr_merge->cd();

  //hSignal_merge->Scale(1./hNtrig_merge->GetBinContent(2));
  hSignal_merge->Scale(1./((hSignal_merge->GetXaxis()->GetBinWidth(1))*(hSignal_merge->GetYaxis()->GetBinWidth(1))));

  //hMixed_merge->Scale(1./hNtrig_merge->GetBinContent(2));
  hMixed_merge->Scale(1./((hMixed_merge->GetXaxis()->GetBinWidth(1))*(hMixed_merge->GetYaxis()->GetBinWidth(1))));
  //hMixed_merge->Scale(1./hMixed_merge->GetBinContent(hMixed_merge->FindBin(0.,0.)));

  hCorr_merge = (TH2D*)hSignal_merge->Clone("hCorr_merge");
  hCorr_merge->Reset("ICSEM");
  hCorr_merge->Add(hSignal_merge);
  hCorr_merge->Divide(hMixed_merge);
  hCorr_merge->GetXaxis()->SetTitle("#Delta#eta");
  hCorr_merge->GetYaxis()->SetTitle("#Delta#varphi (rad)");
  hCorr_merge->GetZaxis()->SetTitle("C_{2}^{Corr.}(#Delta#eta, #Delta#varphi)");

  hCorr_merge->GetXaxis()->CenterTitle(1);
  hCorr_merge->GetYaxis()->CenterTitle(1);
  hCorr_merge->GetZaxis()->CenterTitle(1);

  hCorr_merge->GetYaxis()->SetTitleSize(0.05);
  hCorr_merge->GetYaxis()->SetLabelSize(0.05);
  hCorr_merge->GetYaxis()->SetTitleOffset(1.9);

  hCorr_merge->GetXaxis()->SetTitleSize(0.05);
  hCorr_merge->GetXaxis()->SetLabelSize(0.05);
  hCorr_merge->GetXaxis()->SetTitleOffset(1.9);

  hCorr_merge->GetZaxis()->SetTitleSize(0.05);
  hCorr_merge->GetZaxis()->SetLabelSize(0.05);
  hCorr_merge->GetZaxis()->SetTitleOffset(1.7);

  hCorr_merge->GetXaxis()->SetNdivisions(505);
  hCorr_merge->GetYaxis()->SetNdivisions(505);
  hCorr_merge->GetZaxis()->SetNdivisions(505);

  hCorr_merge->GetXaxis()->SetRangeUser(-MaxEtaRange, MaxEtaRange);
  hCorr_merge->Draw("surf1fbb");

  //~~~~~~~~~~~phi merge flow projection~~~~~~~~~~~~~~~~~~~~~
  c_dphi_merge_flow = new TCanvas("c_dphi_merge_flow", "c_dphi_merge_flow", 600, 530);
  c_dphi_merge_flow->cd();

  Double_t etaMin_merge_flow = hCorr_merge->GetXaxis()->FindBin(etamin_flow);
  Double_t etaMax_merge_flow = hCorr_merge->GetXaxis()->FindBin(etamax_flow);

  hphi_proj_merge_flow = hCorr_merge->ProjectionY("hphi_proj_merge_flow", etaMin_merge_flow, etaMax_merge_flow);

  int nBins_eta_merge_flow = etaMax_merge_flow - etaMin_merge_flow + 1;
  hphi_proj_merge_flow->Scale(1./(nBins_eta_merge_flow));

  hphi_proj_merge_flow->SetMarkerStyle(kFullCircle);
  hphi_proj_merge_flow->SetMarkerColor(kBlue);
  hphi_proj_merge_flow->SetMarkerSize(1.5);

  hphi_proj_merge_flow->GetXaxis()->SetTitle("#Delta#varphi");
  hphi_proj_merge_flow->GetYaxis()->SetTitle("C_{2}^{Corr.}(#Delta#varphi)");
  hphi_proj_merge_flow->GetXaxis()->CenterTitle(1);
  hphi_proj_merge_flow->GetYaxis()->CenterTitle(1);

  hphi_proj_merge_flow->GetYaxis()->SetTitleFont(42);
  hphi_proj_merge_flow->GetYaxis()->SetTitleSize(0.05);
  hphi_proj_merge_flow->GetYaxis()->SetLabelSize(0.05);
  hphi_proj_merge_flow->GetYaxis()->SetTitleOffset(1.8);

  hphi_proj_merge_flow->GetXaxis()->SetTitleFont(42);
  hphi_proj_merge_flow->GetXaxis()->SetTitleSize(0.05);
  hphi_proj_merge_flow->GetXaxis()->SetLabelSize(0.05);
  hphi_proj_merge_flow->GetXaxis()->SetTitleOffset(0.9);

  hphi_proj_merge_flow->Draw("EP");

  //fit_dphi_merge_flow = new TF1("fit_dphi_merge_flow", "[0]*(1.+2.0*([1]*cos(1.0*x) + [2]*cos(2.0*x) + [3]*cos(3.0*x) + [4]*cos(4.0*x) + [5]*cos(5.0*x) + [6]*cos(6.0*x)))", flow_cal_min , flow_cal_max);
  //fit_dphi_merge_flow = new TF1("fit_dphi_merge_flow", "[0]*(1.+2.0*([1]*cos(1.0*x) + [2]*cos(2.0*x) + [3]*cos(3.0*x) + [4]*cos(4.0*x) + [5]*cos(5.0*x) + [6]*cos(6.0*x) + [7]*cos(7.0*x)+ [8]*cos(8.0*x)+ [9]*cos(9.0*x)+ [10]*cos(10.0*x)+ [11]*cos(11.0*x)+ [12]*cos(12.0*x)+ [13]*cos(13.0*x)+ [14]*cos(14.0*x)+ [15]*cos(15.0*x)+ [16]*cos(16.0*x)+ [17]*cos(17.0*x)+ [18]*cos(18.0*x)+ [19]*cos(19.0*x)+ [20]*cos(20.0*x) ))", flow_cal_min , flow_cal_max);
  fit_dphi_merge_flow = new TF1("fit_dphi_merge_flow", "[0]*(1.+2.0*([1]*cos(1.0*x) + [2]*cos(2.0*x) + [3]*cos(3.0*x) + [4]*cos(4.0*x) + [5]*cos(5.0*x) + [6]*cos(6.0*x) + [7]*cos(7.0*x)+ [8]*cos(8.0*x)+ [9]*cos(9.0*x)+ [10]*cos(10.0*x)+ [11]*cos(11.0*x)+ [12]*cos(12.0*x) ))", flow_cal_min , flow_cal_max);

  fit_dphi_merge_flow->SetLineColor(kRed);
  fit_dphi_merge_flow->SetParameter(0, 10.);
  fit_dphi_merge_flow->SetParameter(1, 1.);
  fit_dphi_merge_flow->SetParameter(2, 1.);
  fit_dphi_merge_flow->SetParameter(3, 1.);
  fit_dphi_merge_flow->SetParameter(4, 0.);
  fit_dphi_merge_flow->SetParameter(5, 0.);
  fit_dphi_merge_flow->SetParameter(6, 0.);

  fit_dphi_merge_flow->SetParameter(7, 0.);
  fit_dphi_merge_flow->SetParameter(8, 0.);
  fit_dphi_merge_flow->SetParameter(9, 0.);
  fit_dphi_merge_flow->SetParameter(10, 0.);
  fit_dphi_merge_flow->SetParameter(11, 0.);
  fit_dphi_merge_flow->SetParameter(12, 0.);
  //fit_dphi_merge_flow->SetParameter(13, 0.);
  //fit_dphi_merge_flow->SetParameter(14, 0.);
  //fit_dphi_merge_flow->SetParameter(15, 0.);
  //fit_dphi_merge_flow->SetParameter(16, 0.);
  //fit_dphi_merge_flow->SetParameter(17, 0.);
  //fit_dphi_merge_flow->SetParameter(18, 0.);
  //fit_dphi_merge_flow->SetParameter(19, 0.);
  //fit_dphi_merge_flow->SetParameter(20, 0.);

  hphi_proj_merge_flow->Fit(fit_dphi_merge_flow, "QSER", "", flow_cal_min, flow_cal_max);
  hphi_proj_merge_flow->Fit(fit_dphi_merge_flow, "QSER", "", flow_cal_min, flow_cal_max);
  hphi_proj_merge_flow->Fit(fit_dphi_merge_flow, "SER", "", flow_cal_min, flow_cal_max);

  double v2_merge_coff = fit_dphi_merge_flow->GetParameter(2);
  double v2_merge_coff_error = fit_dphi_merge_flow->GetParError(2);
  double v3_merge_coff = fit_dphi_merge_flow->GetParameter(3);
  double v3_merge_coff_error = fit_dphi_merge_flow->GetParError(3);

  double v2_merge = TMath::Sqrt(v2_merge_coff);
  double v2_merge_error = v2_merge_coff_error/(2.*v2_merge);

  double v3_merge = TMath::Sqrt(v3_merge_coff);
  double v3_merge_error = v3_merge_coff_error/(2.*v3_merge);

  std::cout<<"v2 and v3 merge : "<<(centmin+centmax)/2.<<"   "<<v2_merge<<"    "<<v2_merge_error<<"   "<<v3_merge<<"   "<<v3_merge_error<<std::endl;

  ofstream outmerged;
  outmerged.open(Form("merged_all_pt_vn_%d_%dcent_3k.txt", centmin, centmax));
  double vn_merge_coff[13] = {0};
  double vn_merge_coff_error[13] = {0};
  double vn_merge[13] = {0};
  double vn_merge_error[13] = {0};
  for(int ni=0; ni<13; ni++){
    vn_merge_coff[ni] = fit_dphi_merge_flow->GetParameter(ni);
    vn_merge_coff_error[ni] = fit_dphi_merge_flow->GetParError(ni);
    vn_merge[ni] = TMath::Sqrt(fabs(vn_merge_coff[ni]));
    vn_merge_error[ni] = vn_merge_coff_error[ni]/(2.*vn_merge[ni]);
    cout<<centmin<<"  "<<centmax<<"  "<<ni<<"  "<<vn_merge_coff[ni]<<"  "<<vn_merge_coff_error[ni]<<"  "<<vn_merge[ni]<<"  "<<vn_merge_error[ni]<<endl;
    if(ni>0){
      outmerged<<centmin<<"  "<<centmax<<"  "<<ni<<"  "<<vn_merge_coff[ni]<<"  "<<vn_merge_coff_error[ni]<<"  "<<vn_merge[ni]<<"  "<<vn_merge_error[ni]<<endl;
    }
  }

  c_corr_merge->SaveAs(Form("./plot/vnCent_Corr_GeneralTrack_effCorr_cent%d_%d.pdf", centmin, centmax));
  c_dphi_merge_flow->SaveAs(Form("./plot/vnCent_DphiFit_GeneralTrack_effCorr_cent%d_%d.pdf", centmin, centmax));
}
