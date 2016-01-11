**** SVM High Energy Physics Interface - SVM-HINT ****

The SVM-HINT software package is a LIBSVM based ROOT interface developed for tackling two-class HEP
classification problems. The software package provides easy to use,
automated methods for SVM hyperparameter optimization by taking Asimov
Significance as the performance measure.
Currently SVM-HINT works with LibSVM v3 [1] and ROOT v5.3-v.6.2 [2] on
scientific linux distributions [3]. 
Therefore, a ROOT installation is required.
An example code is given in the 'svm_hint_tutorial.cpp'. The example
codes gives a basic usage scenario for the SVM-HINT.

Right now SVM-HINT only supports two-class problems with C penalty
parameter. Since it is an abstract factory, other SVM algorithms
included in the LibSVM can be implemented easily. 

Tutorial:
The tutorial provides a simple classification problem to demonstrate
SVM-HINT's capabilities. The ROOT trees (event containers) can be
generated using the command: 
> root maketree.C
This code generates two different ROOT trees: svm_testtrain.root and
svm_eval.root. They contain four features for two different classes. 
SVM-HINT uses a C++ container (svm_container) and requires a size
estimation (number of features/variable and the size of the sample -
if the actual size is more than the estimate the memory is automatically allocated).
> svm_container svm(nvar_est,nevent_est);
Once the variables are read from the trees, they can be assigned to
these containers using 
> svm.set_feature(vars.at(0));
In this example we are using same container for training and test
events (SVM-HINT will automatically seperate them later), and a
seperate container for the evaluation.
After adding all variables of the given event we can push back the
event to the container by using
> svm.set_event(weight);
here the argument of the method ("weight") is the weight of the
individual event (double).
>  svm_interface * csvc = new csvc_interface(nsamp_tot,nbkg_tot,nsig_tot);
is the SVM problem object where we give the total event number, and
training and test event numbers seperately - given that they are in
order (training events + test events).

> svm_analyze stop;
declaration of analyzer object.

> stop.set_filename("tutorial.root");
output file name (note the 'root' extension)
 
> stop.set_svm_interface(csvc);
assigning the problem to the analyzer 

> stop.setup_svm(svm);
assigning the text and training event container
 
> stop.set_eval(svm_eval,nbkg_eval);
assigning the evaluation sample

Now there are two choices if there is an estimate of the hyper
parameters, one can simply:
> stop.Obtain_probabilities(c_value, gamma_value,
estimated_disc_cut_value)

or you can do a grid seach using
> stop.Scan_parameters();


Methods:
  
