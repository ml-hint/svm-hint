//simple class to obtain Asimov significance
#ifndef fom_h
#define fom_h
#include <TH2D.h>
#include <TH1D.h>
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "TGraph.h"
static const double epsilon = 5e-1;
class fom{
public:
  enum fom_type{
    asimov,
    stop
  };
  //constructor only with unc initially set
  fom(double uncertain){ unc = uncertain;}
  //destructor
  ~fom(){};
  //set bkg
  void setBackground(double bkg){background = bkg;}
  //set signal
  void setSignal(double sig){signal = sig;}
  TGraph* ROC(TH1D*,TH1D*);
  
  double maxSignificance(TH1D* sig, TH1D* bkg, bool info, double min_signal = 3., TH1D* cuteff = 0);
  float  maxSignificance(TH1F* sig, TH1F* bkg, bool info, float  min_signal = 3., TH1F* cuteff = 0);
  double maxSignificance(const std::vector<std::vector<double> >& sig, const std::vector<std::vector<double> >& bkg,int &bin, double& error,  bool info = true, double min_signal = 3., TH1D* cuteff = 0);
  double getSignificance(const std::vector<std::vector<double> >& sig, const std::vector<std::vector<double> >& bkg,int bin);
  double getSignificance(TH1D* sig, TH1D* bkg, int cut);
  double unbinned_maxSignificance(const std::vector<std::vector<double> >& sig, const std::vector<std::vector<double> >& bkg,int &bin, double& error,  bool info = true, double min_signal = 3., TH1D* cuteff = 0);
  double getSignificance(fom_type f_type) {if (f_type == asimov) asimovZ();
    else if(f_type == stop) Stop();
    else {std::cout << "wrong fom type, please enter one of the following significance types: asimov, stop" << std::endl; std::exit(0);}
    return significance;}
private:
  //asimovZ with uncertainty: http://www.pp.rhul.ac.uk/~cowan/atlas/cowan_statforum_8may12.pdf
  /*  struct sort_with_index {
    size_t i1, i2;
    const  bool operator()(const vvector<double>& sig,const vvector<double>& bkg) {
      return (i1<sig.size()?sig[i1][0]:bkg[i1-sig.size()][0])>(i2<sig.size()?sig[i2][0]:bkg[i2-sig.size()][0]);
    }
    } sort_index;*/

  void asimovZ ();
  void Stop (); //Standard FOM used in sTop analysis
  double errorS, errorB, errorSig;
  double significance;
  double unc;
  double signal;
  double background;
  void check_vars();
  fom(){};
};

#endif
