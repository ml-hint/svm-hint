/*libsvm interface for the one and two class HEP problems */
#ifndef csvc_interface_h
#define csvc_interface_h
#include "svm.h"
#include "libsvm_container.h"
#include "svm_hint.h"
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
//c-svc
class csvc_interface:public svm_interface {
private:
  void set_gamma_array               (std::vector<double, std::allocator<double> >&, unsigned int);
  virtual void   classify            (const svm_model* given_model, const svm_problem& given_problem, vvector<double>& vSig, vvector<double>& vBkg, int nbkg);
  virtual double parameter_comparison(const vvvector<double>& vSig,const vvvector<double>& vBkg, const vvvector<double>& VarvSig,const vvvector<double>& VarvBkg, const int iteration, const double pre_accur, double C);  
  virtual bool   parameter_decision  ();
  std::vector<double> gamma; 
  double bkg_weight;
  double sig_weight;
  double C_cut;
  double below_cut_inc;
  double above_cut_inc;
  int gamma_array_size; /* this number should be an odd number */
  std::vector<double>* accuracy;
  std::vector<int   >* max_cut_bin;
  constexpr static double kBkg = 1.0;
  constexpr static double kSig = -1.0;
public:
  virtual void set_sample          (const svm_container&, samp_type);
  virtual void obtain_probabilities(double,double,int);
  csvc_interface                   (unsigned n_tot, unsigned bkg_tot, unsigned sig_tot):svm_interface("c-svc",n_tot,bkg_tot,sig_tot){
    gamma_array_size = 9; /* this number should be an odd number */
    gamma.resize(gamma_array_size);
    below_cut_inc    = sqrt(10.);
    above_cut_inc    = sqrt(10.)/2.;
    C_cut            = 0.5;
    nbkg_train       = 0;
    nbkg_test        = 0;
    bkg_weight       = 1.;
    nsig_train       = 0;
    nsig_test        = 0;
    sig_weight       = 1.;    
    nfeature         = 0;
    parameters_set   = false;
    sample_train.x   = (svm_node **) calloc (nsamp_tot,sizeof(svm_node *));
    sample_test.x    = (svm_node **) calloc (nsamp_tot,sizeof(svm_node *));
    sample_train.W   = (double *   ) calloc (nsamp_tot,sizeof(double));
    sample_test.W    = (double *   ) calloc (nsamp_tot,sizeof(double));
    doprobabilitycalc= false;
  };
  virtual void gen_class_index () {
    if(nbkg_train < 1) {
      std::cout<< "ERROR! while generating the class indexes: # of background events is 0 " << std:: endl; 
      exit(2); 
    }
    if(nsig_train < 1) {
      std::cout<< "WARNING! while generating the class indexes: # of signal events is 0 "   << std:: endl;
    }
    if(y ) {delete(y) ;}
    if(y1) {delete(y1);}
    if(y2) {delete(y2);}
    y1 = new std::vector<double> (nbkg_train,kBkg); 
    y2 = new std::vector<double> (nsig_train,kSig);
    y  = new std::vector<double>;
    //merging the class labels for background (1) and signal (-1) 
    y->insert(y->begin(), y2->begin(), y2->end());
    y->insert(y->begin(), y1->begin(), y1->end());
    //checking the size of the vectors
    if(y->size() != nsamp_train)  {
      std::cout << "ERROR! while generating the class indexes: Check size of y (vector of labels) "  << y->size()   << 
	" total events reserved for the training: "                                                  << nsamp_train << std::endl; 
      exit(2);
    }
  };
  virtual void set_indexes(){
    gen_class_index();
    sample_train.l = (unsigned) nsamp_train;
    sample_test.l  = (unsigned) nsamp_test;
    sample_eval.l  = (unsigned) nsamp_eval;
    y_array        = (double* ) calloc(nsamp_train, sizeof(double));
    for(int indy = 0; indy < nsamp_train; indy ++){ y_array[indy] = y->at(indy);}
    sample_train.y = y_array;
  }
  virtual bool split_set_sample            (const svm_container &svm);  
  virtual void parameter_scan              (double, double);
  void set_para_gamma                      (double g = 0.1)       { 
    if(!parameters_set) {set_parameters();}
    para.gamma          = g;
  }
  void set_para_c                          (double C = 2.)        { 
    if(!parameters_set) {set_parameters();} 
    para.C              = C;
  };
  void set_para_kernel                     (int type = RBF)       { 
    if(!parameters_set) {set_parameters();} 
    para.kernel_type    = type;
  };     
  void set_para_cache                      (double cache = 24000.){ 
    if(!parameters_set) {set_parameters();}
    para.cache_size     = cache;
  };     
  svm_parameter* get_parameters            (double C, double g)   { 
    if(!parameters_set) {set_parameters();}
    para.C              = C;
    para.gamma          = g; 
    return &para;
  };
  void set_training_weights                  (double background_weight, double signal_weight){
    bkg_weight          = background_weight;
    sig_weight          = signal_weight;
  };
  void weight_training_sample              (){
    if(!parameters_set) {set_parameters();}
    para.nr_weight      = 2;
    para.weight_label   = (int*) malloc(2);
    para.weight         = (double*) malloc(2);
    para.weight_label[0]=+1;
    para.weight_label[1]=-1;
    para.weight[0]      = 1.;
    para.weight[1]      = bkg_weight/ sig_weight;
  }
  void deep_copy_svm_pro                   (const svm_problem&,int,int,int,svm_problem&); //add this to svm_interface
  double sum_of_weights                    (const svm_problem& tbs);//move this to svm_interface
  void adjust_weights                      (svm_problem& tbs, double t, double adjust);//move this to svm_interface
  double sum_of_weights                    (const svm_problem& tbs, double t);//move this to svm_interface
  int sum_of_index                         (const svm_problem& tbs, int);//move this to svm_interface
  double sum_of_value                      (const svm_problem& tbs, int);//move this to svm_interface
  double sum_of_y_index                    (const svm_problem& tbs, int);//move this to svm_interface
  void get_extrema_features                (const svm_problem &tbs, std::vector<double> & min,  std::vector<double> & max);
  void scale                               (svm_problem& tbs, const std::vector<double> & min, const  std::vector<double> & max);
  virtual void set_parameters() {
    //setting svm model parameters
    parameters_set      = true;
    set_para_gamma();
    set_para_c();
    set_para_kernel();
    para.svm_type       = C_SVC;
    para.kernel_type    = RBF;
    para.coef0          = 0.;
    para.eps            = 0.001;
    para.nr_weight      = 0;
    para.degree         = 0;
    para.cache_size     = 24000.;
    para.shrinking      = 1;
    para.probability    = doprobabilitycalc;
    svm_pro->set_parameter_samp(para);
  };
  ~csvc_interface(){
    for(unsigned indx = 0; indx < nsamp_train; indx++) { free(sample_train.x[indx]);}
    free(sample_train.x);
    for(unsigned indx = 0; indx < nsamp_test;  indx++) { free(sample_test.x[indx]);}
    free(sample_test.x);
    free(y_array);
    delete y, y1, y2;
  }
  void clean(svm_problem& tbc);
};
#endif
