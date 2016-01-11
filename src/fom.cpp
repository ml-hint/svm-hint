#include "fom.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <iterator>
#include <algorithm>
using namespace std;
void fom::asimovZ(){
  check_vars();
  double varb = background*unc*background*unc; 
  double tot = signal + background;
  double asimovsig = sqrt(2*(tot*log((tot*(varb+background))/((background*background)+tot*varb))-(1/unc/unc)*log(1.+(varb*signal)/(background*(background+varb)))));
  significance = asimovsig;
}
void fom::Stop(){
  check_vars();
  double varb = background*unc*background*unc; 
  double stopsig = signal / sqrt(background+varb);
  significance = stopsig;
}
void fom::check_vars(){
  if(signal     < epsilon) {cout << "WARNING! Number of signal events: " << signal << endl;} 
  if(background < epsilon) {background = epsilon;}
  if(unc        < 0.0001)  {cout << "WARNING! Uncertainty is set to: "   << unc    << endl;}
}
TGraph* fom::ROC(TH1D* sig, TH1D* bkg){
  int max_bin = sig->GetSize() - 2; //substracting over-under flow bins
  if(max_bin != bkg->GetSize() - 2) {cout << " ERROR! Bin numbers are different.[ROC] " << endl; exit(2);} 
  double B[max_bin];
  double S[max_bin];
  double intSig = 0., intBkg = 0.;
  double  integralB = bkg->Integral(1,max_bin);
  double  integralS = sig->Integral(1,max_bin);
  for (int ind = 1; ind < max_bin+1 ; ind++){
      intSig = sig->Integral(ind,  max_bin);
      intBkg = bkg->Integral(ind,  max_bin);
      S[ind] = intSig/integralS; 
      B[ind] = 1. - (intBkg/integralB);  
    }
  TGraph *roc = new TGraph (max_bin, S, B); 
  return roc;
}
double fom::maxSignificance(TH1D* sig, TH1D* bkg, bool info, double min_signal, TH1D* cuteff){
  bool docuteff = false;
  if(cuteff) docuteff = true; 
  double sign = 0.;
  double max_sign = -1.;
  int max_bin = sig->GetSize() - 2; //substracting over-under flow bins
  if(max_bin != bkg->GetSize() - 2) {cout << " ERROR! Bin numbers are different.[maxSignificance] " << endl; exit(2);}
  int min_bin = 0;
  if(docuteff) min_bin = 1;
  else min_bin  =   (int)(bkg->GetMean()*(double)max_bin); // background -> 0 , signal -> 1 
  if(info) cout<<  " mean bin " << min_bin << "   ";
  double intSig = 0., intBkg = 0.;
  int maxBin = 0;
  double given_sig_error = 0.;
  cout<< " max bin " << max_bin << endl;
  for (int ind = min_bin; ind < max_bin+1 ; ind++) {
    intSig = sig->IntegralAndError(ind,  max_bin+1, given_sig_error);
    intBkg = bkg->Integral(ind,  max_bin+1);
    if(intSig - given_sig_error < min_signal) {
      if(info) cout << " integrals sig " << intSig << 
		 " bkg "                      << intBkg << "    " ;
      cout << " Cut for the max significance is on bin " << maxBin   << 
	" with significance "                                 << max_sign <<  endl;
      return max_sign;
    }
    this->setSignal(intSig); this->setBackground(intBkg);
    sign = this->getSignificance(fom::asimov); //AsimovZ is better
    if(docuteff) cuteff->SetBinContent(ind,sign);
    if(max_sign < sign) {max_sign = sign; maxBin = ind;}
    }
  cout << " Cut for the max significance is on bin " << maxBin   << 
    " with significance "                                 << max_sign << endl;
  return max_sign;
}
float fom::maxSignificance(TH1F* sig, TH1F* bkg, bool info, float min_signal, TH1F* cuteff){
  bool docuteff = false;
  if(cuteff) docuteff = true; 
  float sign = 0.;
  float max_sign = -1.;
  int max_bin = sig->GetSize() - 2; //substracting over-under flow bins
  if(max_bin != bkg->GetSize() - 2) {cout << " ERROR! Bin numbers are different.[maxSignificance] " << endl; exit(2);}
  int min_bin = 0;
  if(docuteff) min_bin = 1;
  else min_bin  =   (int)(bkg->GetMean()*(float)max_bin); // background -> 0 , signal -> 1 
  if(info) cout<<  " mean bin " << min_bin << "   ";
  float intSig = 0., intBkg = 0.;
  int maxBin = 0;
  double given_sig_error = 0.;
  cout<< " max bin " << max_bin << endl;
  for (int ind = min_bin; ind < max_bin+1 ; ind++) {
    intSig = sig->IntegralAndError(ind,  max_bin+1, given_sig_error);
    intBkg = bkg->Integral(ind,  max_bin+1);
    if(intSig - given_sig_error < min_signal) {
      if(info) cout << " integrals sig " << intSig << 
		 " bkg "                      << intBkg << "    " ;
      cout << " Cut for the max significance is on bin " << maxBin   << 
	" with significance "                                 << max_sign <<  endl;
      return max_sign;
    }
    this->setSignal(intSig); this->setBackground(intBkg);
    sign = this->getSignificance(fom::asimov); //AsimovZ is better
    if(docuteff) cuteff->SetBinContent(ind,sign);
    if(max_sign < sign) {max_sign = sign; maxBin = ind;}
    }
  cout << " Cut for the max significance is on bin " << maxBin   << 
    " with significance "                                 << max_sign << endl;
  return max_sign;
}
double fom::maxSignificance(const vector<vector<double> >& sig, const vector<vector<double> >& bkg, int &bin, double &error, bool info, double min_signal, TH1D* cuteff){
  const  int max_ind     = 40;
  double cumuS[max_ind]  = {}, cumuB[max_ind] = {};
  double ScumuS          = 0., ScumuB         = 0.;
  double sign            = 0., max_sign       = -1.;
  int maxBin             = 0, bin_position    = 0;
  bool notbump           = false;
  error = 0;
  for(int ind = 0; ind<max_ind; ind++){
    for (auto it = sig.begin(); it != sig.end();it++ ){
      if(it->at(0) >((double)ind+1.0)/(double)max_ind){
	cumuS[ind]+=it->at(1);
	ScumuS+=(it->at(1)*it->at(1));
      }
    }
    for (auto it =bkg.begin();it!=bkg.end();it++){
      if(it->at(0) >((double)ind+1.0)/(double)max_ind){
	cumuB[ind]+=it->at(1);
	ScumuB+=(it->at(1)*it->at(1));
      }
    }    
    if(cumuS[ind]-sqrt(ScumuS) < min_signal) {
      //if(cumuS[ind] < min_signal) {
      bin = maxBin + bin_position/2;
      cout << " cut for the max significance is on bin " << bin+1    <<
	" with significance "                                 << max_sign <<  endl;
      if(info) cout << "Overestimated error for the given significance: " << error << endl;
      return max_sign;
    }
      this->setSignal(cumuS[ind]); this->setBackground(cumuB[ind]);
      sign = this->getSignificance(fom::asimov); //AsimovZ is better
      if(cuteff) { cuteff->SetBinContent(sign,ind+1);}
      if(max_sign < sign) {
	//	error = sqrt(pow(1/(2*sqrt(cumuB[ind]+unc*unc*cumuB[ind]*cumuB[ind])*sqrt(cumuS[ind])),2)*ScumuS+pow(cumuS[ind]*(1+2*cumuB[ind]*unc*unc)/(2*pow(cumuB[ind]+unc*unc*cumuB[ind]*cumuB[ind],1.5)),2)*ScumuB); // this error is obtained assuming Z=sqrt(s/(b+varb^2))
	max_sign = sign; 
	maxBin = ind; 
	bin_position = 0; 
	notbump = true;
      } else if(max_sign == sign && notbump ) {
	bin_position++;
      } else if(max_sign == sign && !notbump) {
	bin_position=0; 
	notbump = true;
      } else if(max_sign  > sign            ) {
	notbump = false;
      }
  }
  bin = maxBin + bin_position/2;
  cout << " cut for the max significance is on bin " << bin+1    << 
    " with significance "                                 << max_sign << endl;
  if(info) cout << "Overestimated error for the given significance: " << error << endl;
  return max_sign;
}
double fom::unbinned_maxSignificance (const std::vector<std::vector<double> >& sig, const std::vector<std::vector<double> >& bkg, int &bin, double &error, bool info, double min_signal, TH1D* cuteff){
  // some little code to sort over 2 vectors
  vector<unsigned> Idx(sig.size()+bkg.size());
  iota( Idx.begin(), Idx.end(), 0 );// 0...nsig-1,nsig,nsig+1,...,nsig+nbkg
  // admittedly unreadable - sort index such that (Idx[i]<sig.size()?sig[Idx[i]][0]:bkg[Idx[i]-sig.size()][0]) is decreasing
  sort(Idx.begin(), Idx.end(), [&sig,&bkg](size_t i1, size_t i2){
      return (i1<sig.size()?sig[i1][0]:bkg[i1-sig.size()][0])>(i2<sig.size()?sig[i2][0]:bkg[i2-sig.size()][0]);
    });
  //  sort(Idx.begin(), Idx.end(), sort_index(sig,bkg));
  double sumS = 0;
  double sumB = 0;
  double sMax = 0;
  int    iMax = 0;
  double sign = 0.;
  for(int i=0;i<Idx.size();i++){
    if(Idx[i]<sig.size()){
      sumS+=sig[Idx[i]][1];
    } else {
      sumB+=bkg[Idx[i]-sig.size()][1];
    }
    if(sumS>min_signal){
      this->setSignal(sumS); this->setBackground(sumB);
      sign = this->getSignificance(fom::asimov);
      if (sign > sMax){
	sMax=sign;
	iMax=i;
      }
    }
  }
  cout << " with significance "                                 << sMax << endl;

  return sMax;
}
