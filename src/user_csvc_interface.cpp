#include "csvc_interface.h"
#include <iostream>
#include <TRandom3.h>
#include <cmath>
#include <TH1D.h>
#include <TH2D.h>
#include "fom.h"
#include "TMath.h"
#include <vector>
#include <unistd.h>
#include "TString.h"
#include "timer.h"
//#include "csvc_constants.h"
using namespace csvc_constants;
using std::cout;
using std::endl;
void csvc_interface::obtain_probabilities(const double c_p , const double g_p, const int eval_bkg, std::vector<int> * output ){
  cout << "\n Classifying the evaluation sample... " << endl;
  cout << " input variables - C: "                   << c_p               << 
    "  gamma: "                                      << g_p               << 
    "  highest significance cut on probability: "    << highest_accur_cut << 
    "  # of bkg events in eval sample: "             << eval_bkg          << endl;
  timer deltat;
  disc_S              = new TH1D ("svm_disc_signal"         , "SVM probability signal"     , disc_nbin  , 0., 1.);
  disc_B              = new TH1D ("svm_disc_background"     , "SVM probability background" , disc_nbin  , 0., 1.);
  cuteff              = new TH1D ("svm_disc_cuteff"         , "SVM Cut Efficiency"         , disc_nbin  , 0., 1.);
  TH1D*  disc_S_roc   = new TH1D ("svm_disc_signal_roc"     , "SVM probability signal"     , 1000       , 0., 1.);
  TH1D*  disc_B_roc   = new TH1D ("svm_disc_background_roc" , "SVM probability background" , 1000       , 0., 1.);
  TH2D*  dist_B       = new TH2D ("svm_dist_background"     , "SVM background"             , 202, 0., 1.01, 202, 0., 1.01);
  TH2D*  dist_S       = new TH2D ("svm_dist_signal"         , "SVM signal"                 , 202, 0., 1.01, 202, 0., 1.01);
  TH2D*  dist_SV      = new TH2D ("svm_dist_SV"             , "SVM SV"                     , 202, 0., 1.01, 202, 0., 1.01);
  double bkg_yield = 0., sig_yield = 0.;
  /* to have proper unc calculation from the weighted histograms */
  disc_S->Sumw2();  
  disc_B->Sumw2();  
  /* svm prob container */
  double prob[2];
  bool pre_doprobabilitycalc = doprobabilitycalc ;
  doprobabilitycalc = true;
  set_parameters();
  set_indexes();
  svm_model* csvc_svm_model= (svm_model*)malloc(sizeof(svm_model));
  cout << " Number of events to be trained:   " << sample_train.l                    << 
    "  with the total weight of : "                  << sum_of_weights(sample_train) <<  endl;
  cout << " Number of events to be evaluated: " << sample_eval.l                     << 
    "  with the total weight of : "                  << sum_of_weights(sample_eval)  <<  endl;
  if(min_features.empty() && max_features.empty()){
    get_extrema_features(sample_train, min_features, max_features);
  }
  scale(sample_train, min_features, max_features, true);
  scale(sample_eval,  min_features, max_features);
  //  adjust_weights(sample_train, kBkg, (double)nsig_train/(double)nbkg_train);
  deltat.start();
  /* train the method */
  csvc_svm_model =  svm_train(&(sample_train), get_parameters(c_p,g_p));
  deltat.stop("Training time");
  sig_yield = 0.;
  bkg_yield = 0.;
  int * svm_sv = (int*)calloc(svm_get_nr_sv(csvc_svm_model),sizeof(int));
  svm_get_sv_indices(csvc_svm_model, svm_sv);
  int sv_ind = 0; 
  for(int comp = 0; comp < nsamp_train; comp++){
    if(svm_sv[sv_ind] == comp){
      dist_SV->Fill(sample_train.x[comp][0].value,sample_train.x[comp][1].value,sample_train.W[comp]);
      sv_ind ++;
    }
    if(nbkg_train > comp){
      dist_B->Fill(sample_train.x[comp][0].value,sample_train.x[comp][1].value,sample_train.W[comp]);
    } else {
      dist_S->Fill(sample_train.x[comp][0].value,sample_train.x[comp][1].value,sample_train.W[comp]);
    }
  }
  for(int comp = 0; comp < nsamp_eval; comp++){
    svm_predict_probability(csvc_svm_model, sample_eval.x[comp], prob);
    if(prob[1] > highest_accur_cut) { 
      if(output) output    -> push_back(1);
    } else{ 
      if(output) output    -> push_back(0);
    }
    if(eval_bkg  > comp){
      disc_B    -> Fill(prob[1], sample_eval.W[comp]);
      disc_B_roc-> Fill(prob[1], sample_eval.W[comp]);
      if(prob[1] > highest_accur_cut) {bkg_yield += sample_eval.W[comp];}
    } else{
      disc_S    -> Fill(prob[1], sample_eval.W[comp]);
      disc_S_roc-> Fill(prob[1], sample_eval.W[comp]);
      if(prob[1] > highest_accur_cut) {sig_yield += sample_eval.W[comp];}
    }
  }
  TFile * tut = new TFile("tutorial_debug2.root","RECREATE");
  dist_B->Write();
  dist_S->Write();
  dist_SV->Write();
  tut->Close();
  /* assume 25% uncertainty */
  fom FOM (systematical_unc);
  if(!pre_doprobabilitycalc){
    FOM.maxSignificance(disc_S, disc_B, true, 0., cuteff);
  } else{
    FOM.maxSignificance(disc_S, disc_B, false, 0., cuteff);
  } 
  roc = FOM.ROC(disc_S_roc, disc_B_roc);
  FOM.setSignal    (sig_yield);
  FOM.setBackground(bkg_yield);
  if(pre_doprobabilitycalc){
    cout << " Significance obtained from the given cut: "    << FOM.getSignificance(fom::asimov)                     << 
    "  signal yield: "                                     << sig_yield                                            << 
    "  background yield: "                                 << bkg_yield                                            << endl; 
  }
}
