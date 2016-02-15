SVM High Energy Physics Interface - SVM-HINT
============================================

The SVM-HINT software package is a LIBSVM based ROOT interface developed for tackling two-class HEP classification problems. The software package provides easy to use, automated methods for SVM hyperparameter optimization by taking Asimov Significance as the performance measure.
Currently SVM-HINT works with LibSVM v3 [1] and ROOT v5.3-v.6.2 [2] on Scientific Linux distributions [3]. 
Therefore, a ROOT installation is required.
An example code is given in the 'svm_hint_tutorial.cpp'. The example code gives a basic usage scenario for the SVM-HINT.

Right now SVM-HINT only supports two-class problems with C penalty parameter. Since it is an abstract factory, other SVM algorithms included in the LibSVM can be implemented easily. 

SVM-HINT's performance and usage scenarios are discussed in [4].

[1] http://www.csie.ntu.edu.tw/~cjlin/libsvm/

[2] https://root.cern.ch/

[3] https://www.scientificlinux.org

[4] http://arxiv.org/pdf/1601.02809v1.pdf

**Installation:**
Make sure that the ROOT path is set correctly:

[sahin@naf-uhhcms06:~] echo $ROOTSYS

/afs/desy.de/products/root/amd64_rhel60/5.34.00

inside of the svm_interface directory: 

[sahin@naf-uhhcms06:~/svm_interface] make

g++ -O3 -fPIC -fopenmp -std=c++0x  -c   ./libsvm-weights-3.20//svm.cpp

g++ -O3 -fPIC -fopenmp -std=c++0x  svm_hint_tutor.cpp svm.o fom.cpp csvc_interface.cpp -I/afs/desy.de/products/root/amd64_rhel60/5.34.00/include -L/afs/desy.de/products/root/amd64_rhel60/5.34.00/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic  -o tutor_svm -lm 

to run tutorial:

[sahin@naf-uhhcms06:~/svm_interface] root maketree.C

[sahin@naf-uhhcms06:~/svm_interface] ./tutor_svm


**Tutorial**

The tutorial provides a simple classification problem to demonstrate
SVM-HINT's capabilities. The ROOT trees (event containers) can be
generated using the command: 

[sahin@naf-uhhcms06:~/svm_interface] root maketree.C

This code generates two different ROOT trees: svm_testtrain.root and
svm_eval.root. They contain four features for two different classes. 
SVM-HINT uses a C++ container (svm_container) and requires a size
estimation (number of features/variable and the size of the sample -
if the actual size is more than the estimate the memory is automatically allocated).

$ svm_container svm(nvar_est,nevent_est);

Once the variables are read from the trees, they can be assigned to
these containers using 

$ svm.set_feature(vars.at(0));

In this example we are using same container for training and test events (SVM-HINT will automatically separate them later), and a
separate container for the evaluation.

After adding all variables of the given event we can push back the
event to the container by using

$ svm.set_event(weight);

here the argument of the method ("weight") is the weight of the
individual event (double).

$  svm_interface * csvc = new csvc_interface(nsamp_tot,nbkg_tot,nsig_tot);

nsamp_tot: Total number of events to be trained and tested

nbkg_tot: Total number of background events

nsig_tot: Total number of signal events

is the SVM problem object where we give the total event number, and
training and test event numbers separately - given that they are in
order (training events + test events).

$ svm_analyze stop;

declaration of analyzer object.

$ stop.set_filename("tutorial.root");

output file name (note the 'root' extension)
 
$ stop.set_svm_interface(csvc);

assigning the problem to the analyzer 

$ stop.setup_svm(svm);

assigning the text and training event container
 
$ stop.set_eval(svm_eval,nbkg_eval);

assigning the evaluation sample

Now there are two choices if there is an estimate of the hyper
parameters, one can simply:

$ stop.Obtain_probabilities(c_value, gamma_value, estimated_disc_cut_value)

or you can do a grid search using

$ stop.Scan_parameters();


 **Short class definitions**

svm_interface: 
Abstract base class deriving: csvc_interface 

csvc_interface : 
The LIBSVM-csvc interface class. It holds the model containers of the given SVM instance. 

svm_analyze:
SVM problem analyzer

fom:
Calculates various figure of merits (i.e. Z_Asimov, Z_{S/sqrt{S+B}}).
