/*libsvm interface for the one and two class HEP problems */
#ifndef svm_hint_h
#define svm_hint_h
#include "svm.h"
#include "libsvm_container.h"
#include <cstdlib>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <TH1D.h>
#include <TFile.h>
#include <TGraph.h>
#include <omp.h>
#include <thread>
#include "fom.h"
template<typename T >  using vvector  = std::vector<std::vector<T> >;
template<typename T >  using vvvector = std::vector<std::vector<std::vector<T> > > ;
//svm problem container class
class svm_problem_container {
public:
  void set_training_samp  (const svm_problem & set)   {sample_train = set;};
  void set_test_samp      (const svm_problem & set)   {sample_test  = set;};
  void set_parameter_samp (const svm_parameter & set) {para         = set;};
private:
  svm_model * model;
  svm_problem sample_train;
  svm_problem sample_eval;
  svm_problem sample_test;
  svm_parameter para;
};
//abstract class svm interface
class svm_interface {
protected:
  svm_problem para_scan_problem_test;
  svm_problem para_scan_problem_train;
  std::string svm_type;
  svm_problem_container *svm_pro;
  static const int disc_nbin = 40;
  int max_iter;
  bool precise_scan;
  bool doprobabilitycalc;
  unsigned nsamp_tot, nsamp_train, nsamp_test, nsamp_eval, nbkg_tot, nsig_tot , nbkg_train, nbkg_test, nsig_train, nsig_test, nfeature;
  double nsig_test_w, nbkg_test_w, systematical_unc;
  vvector<double> overall_accuracy;   
  std::vector<double> min_features;
  std::vector<double> max_features;
  svm_interface(std::string type,unsigned n_tot, unsigned bkg_tot, unsigned sig_tot)
    :svm_type(type), nsamp_tot(n_tot), nbkg_tot(bkg_tot), nsig_tot(sig_tot){
    svm_pro             = new svm_problem_container;
    max_iter            = 12;
    precise_scan        = true;
    highest_accur_gamma = 0.1;
    highest_accur_C     = 0.1;
    highest_accur_cut   = 0.5;
    y                   = 0;
    y1                  = 0;
    y2                  = 0;
    systematical_unc    = 0.25;
 };
  void do_precise_scan    (bool precise  = true) {precise_scan = precise;};
  unsigned svm_node_max;
  svm_problem sample_train;
  svm_problem sample_test;
  svm_problem sample_eval;
  svm_parameter para; 
  bool parameters_set;
  double min_signal;
public:  
  //index vectors for the two class problem
  std::vector<double> * y1;
  std::vector<double> * y2;
  std::vector<double> * y ;
  double* y_array;
  double highest_accur_gamma;
  double highest_accur_C;
  double highest_accur_cut;
  enum samp_type {TRAIN,TEST,EVALUATE};
  TGraph* roc;
  TH1D* disc_S;
  TH1D* disc_B;
  TH1D* cuteff;
  virtual void do_probability_calc    (bool probcalc = true) {doprobabilitycalc = probcalc;};
  virtual void set_sample             (const svm_container&, samp_type) = 0;
  virtual bool split_set_sample       (const svm_container &svm) = 0;  
  virtual void set_systematical_unc   (double sys_unc = 0.25) {systematical_unc = sys_unc;}
  virtual void set_parameters         ()=0;  
  virtual void gen_class_index        ()=0;
  virtual void set_indexes            ()=0;
  virtual void obtain_probabilities   (double,double,int) = 0;
  virtual void parameter_scan         (double,double) = 0;
  virtual void classify               (const svm_model *, const svm_problem &, vvector<double> &, vvector<double> &, int) = 0;
  virtual double parameter_comparison (const vvvector<double> &,const vvvector<double> &, const vvvector<double> &,const vvvector<double> &, const int, const double, double) = 0;
  virtual bool parameter_decision     () = 0;
  virtual ~svm_interface              (){delete svm_pro;};
};
class svm_analyze{
private:
  svm_interface * given_svm;  
  bool svm_inter_set;
  int iEval_bkg;
  TString filename;
  void auto_set_core_n(){
    static unsigned int omp_ncore;
    omp_ncore = std::thread::hardware_concurrency();
    if(0 == omp_ncore) {omp_ncore = sysconf( _SC_NPROCESSORS_ONLN );}
    /* oh well, assuming the most CPU s these days have 8 cores */
    if(0 == omp_ncore) {omp_ncore = 8;}
    std::cout << " Number of threads in OMP set to " << omp_ncore - 3 << std::endl; 
    omp_set_num_threads(12);
    //    omp_set_num_threads(1);
  }
  public:
  svm_analyze           ():svm_inter_set(false),filename("svm.root"){ auto_set_core_n();};
  void set_svm_interface(svm_interface *svm_int){
    given_svm = svm_int;
    svm_inter_set=true;
  };
  bool setup_svm        (svm_container& svm){
    if(!svm_inter_set) {
      std::cout<<"ERROR! svm_interface is not set in the analysis class! "<< std::endl;
      exit(2);
    }
    given_svm->set_parameters();
    given_svm->split_set_sample(svm);
    given_svm->set_indexes();
    return true;
  };
  void set_eval(const svm_container &eval, int eval_bkg){
    iEval_bkg = eval_bkg;
    given_svm->set_sample(eval,svm_interface::EVALUATE);
  };
  void set_filename(TString name){filename=name;}
  void Scan_parameters(){
    given_svm->parameter_scan(0.0001, 0.); 
    std:: cout << given_svm->highest_accur_C << " "<<given_svm->highest_accur_gamma << std::endl;
    this->Obtain_probabilities(given_svm->highest_accur_C,given_svm->highest_accur_gamma);
  }
  void Do_probability_calc(){
    given_svm->do_probability_calc();
    given_svm->set_parameters();
  }
  void Obtain_probabilities(double C, double g) const{
    given_svm->obtain_probabilities((double)C,(double)g,iEval_bkg);
    TFile *f = new TFile(filename,"RECREATE");
    given_svm->disc_S->Write();
    given_svm->disc_B->Write();
    given_svm->cuteff->Write();
    given_svm->roc->Write();
    f->Close();
  }
  void Obtain_probabilities(double C, double g, double high_cut){
    given_svm->highest_accur_cut = high_cut; 
    Obtain_probabilities((double)C, (double)g);
  }
};
#endif
