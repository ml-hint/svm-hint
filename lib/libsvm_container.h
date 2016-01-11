#ifndef libsvm_container_h
#define libsvm_container_h
#include"svm.h"
#include <vector>
#include <iostream>
class svm_container{
public:
  svm_container(int nVar,int nEntries):nvar(nVar),nentries(nEntries) 
  {
    svm_cont = new std::vector<svm_node_array>[(nvar+1)*nentries*sizeof(svm_node)];
    feat_counter  = 1;
    weights = new std::vector<double>[nentries];
  };
  void set_feature(double feat_val)
  {
    if(feat_counter == 1) event = new svm_node_array;
    feat.index = feat_counter;
    feat.value = feat_val;
    feat_counter ++;
    event->push_back(feat);
  }
  void set_event (double weight)
  {
    feat.index = -1;
    if (feat_counter != nvar +1){
      std::cout<< "ERROR! Feature number does not match input! " << std::endl;
      exit(2);
    }  
    event->push_back(feat);
    svm_cont->push_back((*event));
    weights->push_back(weight);
    feat_counter = 1;
    delete event;
  }
  typedef std::vector<svm_node> svm_node_array;
  std::vector<svm_node_array>* svm_cont;
  std::vector<double>* weights;
  ~svm_container(){};//{delete svm_cont;};  
private:
  int nvar, nentries; 
  svm_container(){};
  svm_node feat;
  svm_node_array *  event;
  int feat_counter;
};

#endif
