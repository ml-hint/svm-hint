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
#define  fEpsilon 1.e-2
using std::cout;
using std::endl;
void csvc_interface::obtain_probabilities(double c_p , double g_p, int eval_bkg){
  cout << "Entering the obtain probabilities function... "   << endl;
  cout << " input variables - C: "  << c_p               << 
    "  gamma: "                          << g_p               << 
    "  highest accuracy cut: "           << highest_accur_cut << 
    "  # of bkg events in eval sample: " << eval_bkg          << endl;
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
  scale(sample_train, min_features, max_features);
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
  FOM.maxSignificance(disc_S, disc_B, true, 0., cuteff);
  roc = FOM.ROC(disc_S_roc, disc_B_roc);
  FOM.setSignal    (sig_yield);
  FOM.setBackground(bkg_yield);
  cout << " Significance obtained from the given cut: " << FOM.getSignificance(fom::asimov)                          << 
    "  signal yield: "                                     << sig_yield                                            << 
    "  background yield: "                                 << bkg_yield                                            << endl; 
}
void csvc_interface::set_gamma_array(std::vector<double> & gamma_arr, unsigned iter){
  int mid = (gamma_arr.size() - 1)/2;
  static int f = 0;
  /* after every fourth iteration the gamma array focuses on the gamma value which gave the highest significance in the previous iterations */
  if((iter%4)==0){
    f++;
    gamma_arr.at(mid) = highest_accur_gamma;
    for(int ind = 0; ind < (gamma_arr.size() - 1)/2; ind++){
      gamma_arr.at(mid-1-ind) = gamma_arr.at(mid-ind)* 2./(log(exp(1)*(double)f/2.));
      gamma_arr.at(mid+1+ind) = gamma_arr.at(mid+ind)/(2./ log(exp(1)*(double)f/2.));
    }
    cout << "gamma array set to: {";
    for(int ind = 0;ind < gamma_arr.size()-1; ind++)
      cout << gamma_arr[ind] << ", ";
    cout << gamma_arr[gamma_arr.size()-1] << "}" << endl;
  }
}
void csvc_interface::parameter_scan(double C, double pre_accur){ 
  static int iteration = 0;
  cout << "Entering the parameter scan function... Iteration: " << iteration << endl; 
  if(!iteration){
    deep_copy_svm_pro(sample_test , nsamp_test , svm_node_max, 0          , para_scan_problem_test );
    deep_copy_svm_pro(sample_train, nsamp_train, svm_node_max, nsamp_train, para_scan_problem_train);
    get_extrema_features(sample_train, min_features, max_features);
    scale(para_scan_problem_test , min_features, max_features);
    scale(para_scan_problem_train, min_features, max_features);
  }
   /*dynamically setting gamma array */
  set_gamma_array(gamma,iteration);
  int maxbin[gamma_array_size];
  int highest = 0;
  std::vector<svm_parameter> scan_parameters;
  for(std::vector<double>::iterator g = gamma.begin(); g != gamma.end(); g++){
    scan_parameters.push_back(*get_parameters(C,*g));
  }
  double sig_error = 0.,dummy_error = 0.;
  static std::vector<std::vector<std::vector<double> > > vSig_test (gamma_array_size, std::vector<std::vector<double> >(nsig_test,  std::vector<double>(2)));
  static std::vector<std::vector<std::vector<double> > > vBkg_test (gamma_array_size, std::vector<std::vector<double> >(nbkg_test,  std::vector<double>(2)));
  static std::vector<std::vector<std::vector<double> > > vSig_train(gamma_array_size, std::vector<std::vector<double> >(nsig_train, std::vector<double>(2)));
  static std::vector<std::vector<std::vector<double> > > vBkg_train(gamma_array_size, std::vector<std::vector<double> >(nbkg_train, std::vector<double>(2)));
  svm_model * csvc_svm_model[gamma_array_size];
  for(int indg = 0; indg < gamma_array_size; indg++){
      csvc_svm_model[indg] = svm_train(&(para_scan_problem_train), &(scan_parameters.at(indg)));
  }
  int indg ;
  //#pragma omp parallel for shared(vSig,vBkg) private (indg) 
  for( indg = 0; indg < gamma_array_size; indg++){
    classify(csvc_svm_model[indg], para_scan_problem_test,  vSig_test[indg],  vBkg_test[indg],  nbkg_test );
    classify(csvc_svm_model[indg], para_scan_problem_train, vSig_train[indg], vBkg_train[indg], nbkg_train);
    svm_free_and_destroy_model(&(csvc_svm_model[indg]));    
  }
  iteration++;
  pre_accur = parameter_comparison(vSig_test, vBkg_test, vSig_train, vBkg_train, iteration, pre_accur, C);
}
void csvc_interface::classify(const svm_model* given_model, const svm_problem& given_problem, vvector<double>& vSig, vvector<double>& vBkg, int nbkg){
  double prob[] = {0.,0.} ;   
  double dec[]  = {0.,0.} ;   
  for(int comp = 0; comp < given_problem.l; comp++){
    if(doprobabilitycalc) {
      svm_predict_probability(given_model,given_problem.x[comp],prob);
    } else{
      svm_predict_values(given_model,given_problem.x[comp],dec);
      prob[1] = -dec[0];                              // works as well - do not miss the sign!
      prob[0] = 1-prob[1];
    }
    if(nbkg > comp){ 
      vBkg.at(comp).at(0) = prob[1]; 
      vBkg.at(comp).at(1) = given_problem.W[comp];
    } else { 
      vSig.at(comp-nbkg).at(0) = prob[1]; 
      vSig.at(comp-nbkg).at(1) = given_problem.W[comp];
    }
  }
}
double csvc_interface::parameter_comparison(const vvvector<double>& vSig, const vvvector<double>& vBkg, const vvvector<double>& VarvSig, const vvvector<double>& VarvBkg, const int iteration, const double pre_accur, double C){
  cout << " Entering parameter comparison " << endl;
  static int    plato_counter = 0;
  accuracy    = new std::vector<double>(gamma_array_size);
  max_cut_bin = new std::vector<int>   (gamma_array_size); 
  std::vector<double> high_accuracy_settings (4); /* accuracy, gamma, C, bin */ 
  double sig_error = 0., dummy_error = 0., train_accur = 0.;
  int    highest = 0 ,     dummy_bin = 0;
  //assume 25% unc
  fom    FOM(systematical_unc);
  /*obtaining significance for each value in the gamma array. 
    This can be included inside of the above loop but due to OMP thread exclusivity, it is better to run it outside */
  for(int ind = 0; ind < gamma_array_size; ind++){
    cout << "Gamma " << gamma.at(ind) << endl;
    if(!doprobabilitycalc){
      cout << " Test sample: " ;
      accuracy->at(ind) = FOM.unbinned_maxSignificance(vSig[ind], vBkg[ind], max_cut_bin->at(ind), sig_error, false, min_signal);  
      cout << " Train sample: " ;
      train_accur = FOM.unbinned_maxSignificance(VarvSig[ind], VarvBkg[ind], dummy_bin , dummy_error, false, min_signal);
    } else {
      cout << " Test sample: " ;
      accuracy->at(ind) = FOM.maxSignificance(vSig[ind], vBkg[ind], max_cut_bin->at(ind), sig_error, false, min_signal);  
      cout << " Train sample: " ;
      train_accur = FOM.getSignificance(VarvSig[ind], VarvBkg[ind], max_cut_bin->at(ind));
    }
    accuracy->at(ind) *= (1. - fabs(accuracy->at(ind) - train_accur)/(accuracy->at(ind) + train_accur));
  }
  sig_error;
  for(int ind = 0; ind<gamma_array_size-1; ind++){
    if( accuracy->at(ind) > accuracy->at(ind+1) || ( accuracy->at(ind) == accuracy->at(ind+1) && max_cut_bin->at(ind) < max_cut_bin->at(ind+1) ) ) {
      accuracy->at(ind+1)    = accuracy->at(ind); 
      max_cut_bin->at(ind+1) = max_cut_bin->at(ind);
    } else { 
      highest   = ind + 1;
      sig_error = dummy_error; //propagate only the error of maximum significance
    }
  }
  //this will ensure the iteration algorithm does not stuck in a plato
  if(accuracy->at(highest) == pre_accur){
    plato_counter++;
  } else{
    plato_counter = 0;
  }
  /* for the low C values keep the fluctions allowance higher (0.5) for the rest 1 sigma band */
  if( ( ( accuracy->at(highest) > ( pre_accur * 0.68 ) ) || ( C < 0.6 && accuracy->at(highest) > ( pre_accur * 0.5 )) ) && iteration < max_iter && plato_counter < 5 ){
    high_accuracy_settings[0] = accuracy->at(highest);
    high_accuracy_settings[1] = (double)(max_cut_bin->at(highest)+1.)/(double)disc_nbin;
    high_accuracy_settings[2] = C;
    high_accuracy_settings[3] = gamma[highest];
    highest_accur_gamma       = gamma[highest];
    cout <<" For given C value: "      << std::setprecision(16) << C                        << 
      "  highest accuracy: "           << std::setprecision(6)  << accuracy->at(highest)    << 
      "  obtained with gamma: "        << std::setprecision(16) << gamma [highest]          << 
      "  on the bin "                  << std::setprecision(6)  << max_cut_bin->at(highest) << 
      "  in the iteration "            << iteration             << endl;
    overall_accuracy.push_back(high_accuracy_settings);
    if(iteration == max_iter - 5 && precise_scan) { 
      parameter_decision();
    } else { 
      if      (C < 0.5) {parameter_scan(below_cut_inc * C, accuracy->at(highest));}
      else              {parameter_scan(above_cut_inc * C, accuracy->at(highest));}
    }
  } else {
    if(precise_scan) max_iter = iteration + 5;
    parameter_decision();
    return 1.;
  }
}
bool csvc_interface::parameter_decision()   {
  int max_ind = 0; 
  auto it_max = overall_accuracy.begin();
  for(auto it = overall_accuracy.begin(); it != overall_accuracy.end() -1;it++ ){
    if(it_max->at(0) < (it+1)->at(0)) {it_max = it+1;} 
  }
  highest_accur_cut   = it_max->at(1);
  highest_accur_C     = it_max->at(2);
  highest_accur_gamma = it_max->at(3);
  if(precise_scan) {
    do_precise_scan(false);
    if (highest_accur_C > C_cut && highest_accur_C / above_cut_inc > C_cut){
      above_cut_inc -= (above_cut_inc - 1.) * 4. / 5.;
      below_cut_inc -= (below_cut_inc - 1.) * 4. / 5.;
      cout << " entering precise scan with initial C " << highest_accur_C / (above_cut_inc * above_cut_inc) << " and previous accuracy " << it_max->at(0) << endl;
      this->parameter_scan(highest_accur_C - 2. * above_cut_inc , it_max->at(0));
    } else {
      above_cut_inc -= (above_cut_inc - 1.) * 4. / 5.;
      below_cut_inc -= (below_cut_inc - 1.) * 4. / 5.;
      cout << " entering precise scan with initial C " << highest_accur_C / (below_cut_inc * below_cut_inc) << " and previous accuracy " << it_max->at(0) << endl;
      this->parameter_scan(highest_accur_C / (below_cut_inc * below_cut_inc) , it_max->at(0));
    }
  } else {
    cout << " Optimized C "      << std::setprecision(16) << std::fixed << highest_accur_C      << 
      " optimized gamma "             << std::setprecision(16) << std::fixed << highest_accur_gamma  << 
      " optimized discriminator cut " << std::setprecision(4)  << std::fixed << highest_accur_cut    << 
      " with the accuracy of "        << std::setprecision(4)  << std::fixed << it_max->at(0)        << 
      " in the test sample "          << std::setprecision(12) << endl;
    clean(para_scan_problem_test);
    clean(para_scan_problem_train);
    return true;
  }
} 
void csvc_interface::clean(svm_problem& tbc){
  free(tbc.W);
  free(tbc.x);
}
void csvc_interface::deep_copy_svm_pro(const svm_problem& tbc, int x1_max, int x2_max, int y_max, svm_problem& copy){
  copy.x = (svm_node **) calloc(x1_max,sizeof(svm_node*)); 
  for (unsigned indx = 0; indx < x1_max; indx++){ 
    copy.x[indx] = (svm_node*) calloc(x2_max,sizeof(svm_node)); 
    for(unsigned indy = 0; indy < x2_max; indy++) {
      copy.x[indx][indy] = tbc.x[indx][indy]; 
    } 
  }
  copy.y = (double*)calloc(y_max,sizeof(double));
  for(unsigned ind =0; ind < y_max; ind++) {
    copy.y[ind] = tbc.y[ind];
  }
  copy.W = (double*)calloc(x1_max,sizeof(double));
  for(unsigned ind =0; ind < x1_max; ind++) {
    copy.W[ind] = tbc.W[ind];
  }
  copy.l = tbc.l;
};
double csvc_interface::sum_of_weights (const svm_problem& tbs){
  double sum_of_w = 0.;
  for(unsigned ind = 0; ind < tbs.l; ind++){
    sum_of_w += tbs.W[ind];
  }
  return sum_of_w;
};
double csvc_interface::sum_of_weights (const svm_problem& tbs, double t){
  double sum_of_w = 0.;
  for(unsigned ind = 0; ind < tbs.l; ind++){
    if(tbs.y[ind] == t ) sum_of_w += tbs.W[ind];
  }
  return sum_of_w;
};
double csvc_interface::sum_of_value   (const svm_problem& tbs, int index){
  if(index > nfeature+1) return -1.;
  double sum_of_i = 0.;
  for(unsigned ind = 0; ind < tbs.l; ind++){
    sum_of_i += tbs.x[ind][index-1].value;
  }
  return sum_of_i;
};
int csvc_interface::sum_of_index      (const svm_problem& tbs, int index){
  if(index > nfeature+1) return -1.;
  int sum_of_i = 0.;
  for(unsigned ind = 0; ind < tbs.l; ind++){
    sum_of_i += tbs.x[ind][index-1].index;
  }
  return sum_of_i;
};
double csvc_interface::sum_of_y_index (const svm_problem& tbs, int index){
  if(index > nfeature) return -1.;
  double sum_of_i = 0.;
  for(unsigned ind = 0; ind < tbs.l; ind++){
    sum_of_i += tbs.x[ind][index-1].value * (double)tbs.y[ind] * tbs.W[ind];
  }
  return sum_of_i;
};
void csvc_interface::adjust_weights(svm_problem& tbs, double t, double adj){
  cout << " Weights of the label " << t << " are adjusted with factor " << adj << endl;
  for(unsigned ind = 0; ind < tbs.l; ind++){
    if(tbs.y[ind] == t ) tbs.W[ind] *= adj;
  }
};
void csvc_interface::get_extrema_features (const svm_problem &tbs, std::vector<double>& min,  std::vector<double>& max){
  for(unsigned ind = 0; ind < tbs.l; ind++){
    for( unsigned feat = 0; feat < nfeature; feat++){
      if(0 == ind){
	min.push_back( tbs.x[ind][feat].value);
	max.push_back( tbs.x[ind][feat].value);
      } else {
	min.at(feat) = min.at(feat) < tbs.x[ind][feat].value ? min.at(feat) : tbs.x[ind][feat].value;
	max.at(feat) = max.at(feat) > tbs.x[ind][feat].value ? max.at(feat) : tbs.x[ind][feat].value;
      }
    }
  }
}
void csvc_interface::scale                (svm_problem &tbs, const std::vector<double> & min, const std::vector<double> & max){
  cout << " Entering scale function: " << endl;
  for (unsigned feat = 0; feat < nfeature; feat++) {
    cout << " feature: " << feat + 1 << " min: " << min.at(feat) << "  max: " << max.at(feat) << endl; 
  }
  for(unsigned ind = 0; ind < tbs.l; ind++) {
    for( unsigned feat = 0; feat < nfeature; feat++) {
      tbs.x[ind][feat].value = (tbs.x[ind][feat].value - min.at(feat)) / (max.at(feat) - min.at(feat));
    }
  }
}
bool csvc_interface::split_set_sample (const svm_container & svm){
  //c++ std vectors should be converted to the plain C arrays - allocating memory
  unsigned indtra  = 0;
  unsigned indtest = 0;
  cout << " Entering the sample splitter... "                              << endl;
  cout << " Number of events to be seperated: "               << nsamp_tot << 
    "  background events: "                                   << nbkg_tot  << 
    "  signal events: "                                       << nsig_tot  << endl;
  nbkg_test_w = 0.;
  nsig_test_w = 0.;
  int nrand   = 0;
  for(unsigned indx= 0 ; indx < nsamp_tot; indx++){
    //allocating memory for the node array
    //    nrand = sep.Uniform(1.,3.);
    if(indx < nbkg_tot){
      if(1 == indx%2){
	nbkg_test_w += svm.weights->at(indx);
	sample_test.x[indtest] = (svm_node*) calloc(svm.svm_cont->at(0).size(), sizeof(svm_node));
	sample_test.W[indtest] = svm.weights->at(indx);
	for(unsigned indy = 0; indy <svm.svm_cont->at(indx).size() ; indy ++) {
	  sample_test.x[indtest][indy] = (svm.svm_cont->at(indx).at(indy));
	}  
	indtest++;
      } else {
	sample_train.x[indtra] = (svm_node*) calloc(svm.svm_cont->at(0).size(), sizeof(svm_node));
	sample_train.W[indtra] = svm.weights->at(indx);
	for(unsigned indy = 0; indy < svm.svm_cont->at(indx).size() ; indy ++) {
	  sample_train.x[indtra][indy] = (svm.svm_cont->at(indx).at(indy)); 
	} 
	nbkg_train++; 
	indtra++;
      }
    } else {
      if(1 == indx%2) {
	nsig_test_w += svm.weights->at(indx);
	sample_test.x[indtest] = (svm_node*) calloc(svm.svm_cont->at(0).size(), sizeof(svm_node));
	sample_test.W[indtest] = svm.weights->at(indx);
	for(unsigned indy = 0; indy <svm.svm_cont->at(indx).size(); indy ++) {
	  sample_test.x[indtest][indy] = (svm.svm_cont->at(indx).at(indy));
	}  
	indtest++;
      } else {
	sample_train.x[indtra] = (svm_node*) calloc(svm.svm_cont->at(0).size(), sizeof(svm_node));
	sample_train.W[indtra] = svm.weights->at(indx);
	for(unsigned indy = 0; indy < svm.svm_cont->at(indx).size(); indy ++) {
	  sample_train.x[indtra][indy] = (svm.svm_cont->at(indx).at(indy)); 
	}
	nsig_train++; 
	indtra++;
      }
    }
  } 
  cout << " Sum of event weights for the test samples:"      << 
    " Signal "      << nsig_test_w                           << 
    " Background "  << nbkg_test_w                           << endl;
  nsamp_train  = nsig_train+nbkg_train; 
  svm_node_max = svm.svm_cont->at(0).size();
  nfeature     = svm_node_max - 1; 
  cout << " Number of features: "                << nfeature << endl;
  /* take minimum signal efficiency of 1%, if it is less than 1, take one signal event as minimum instead */
  min_signal   = nsig_test_w * 0.01; 
  if(min_signal < 1.) min_signal = 1.;  
  nsamp_test   = nsamp_tot - nsamp_train;
  nsig_test    = nsig_tot  -  nsig_train;
  nbkg_test    = nbkg_tot  -  nbkg_train;
  svm_pro  ->  set_test_samp (sample_test);
  svm_pro  ->  set_training_samp(sample_train);
  return true;
}
void csvc_interface::set_sample (const svm_container & svm, samp_type type) {
  svm_problem sample;
  sample.x = (svm_node **) calloc(svm.svm_cont->size(), sizeof(svm_node *));
  sample.W = (double   * ) calloc(svm.svm_cont->size(), sizeof(double    ));
  unsigned indx= 0;
  for(indx = 0; indx < svm.svm_cont->size(); indx++){ 
    sample.x[indx] = (svm_node*) calloc(svm.svm_cont->at(0).size(),sizeof(svm_node));
    sample.W[indx] = svm.weights->at(indx);	    
    for(unsigned indy = 0; indy <svm.svm_cont->at(indx).size() ; indy ++) {
      sample.x[indx][indy] = (svm.svm_cont->at(indx).at(indy));
    }
  }
  if(EVALUATE     == type){
    sample_eval   = sample;
    nsamp_eval    = indx  ;
  } else if(TRAIN == type){
    sample_train  = sample;
    nsamp_train   = indx  ;
  } else if(TEST  == type){
    sample_test   = sample;
    nsamp_test    = indx  ;
  }
}
