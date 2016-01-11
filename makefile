CXX ?= g++
CFLAGS = -O3 -fPIC -fopenmp -std=c++0x 
CC = icpc
SHVER = 2
OS = $(shell uname)

LIBSVM_DIR = ./libsvm-weights-3.20/

SVMHINT_SRC = ./src/
SVMHINT_LIB = ./lib/
SVMHINT_TEST = ./test/

all: tutorial
lib: svm.o
	if [ "$(OS)" = "Darwin" ]; then \
		SHARED_LIB_FLAG="-dynamiclib -Wl,-install_name,libsvm.so.$(SHVER)"; \
	else \
		SHARED_LIB_FLAG="-shared -Wl,-soname,libsvm.so.$(SHVER)"; \
	fi; \
	$(CXX) $${SHARED_LIB_FLAG} svm.o -o libsvm.so.$(SHVER)

svm.o:$(LIBSVM_DIR)/svm.cpp $(LIBSVM_DIR)/svm.h
	$(CXX) $(CFLAGS) -c   $(LIBSVM_DIR)/svm.cpp
tutorial:$(SVMHINT_TEST)/svm_hint_tutor.cpp $(SVMHINT_SRC)/fom.cpp $(SVMHINT_LIB)/fom.h $(SVMHINT_LIB)/libsvm_container.h  $(SVMHINT_LIB)/svm_hint.h $(SVMHINT_LIB)/csvc_interface.h svm.o
	$(CXX) $(CFLAGS) $(SVMHINT_TEST)/svm_hint_tutor.cpp svm.o $(SVMHINT_SRC)/fom.cpp $(SVMHINT_SRC)/csvc_interface.cpp -I$(LIBSVM_DIR) -I$(SVMHINT_LIB) -I$(ROOTSYS)/include -L$(ROOTSYS)/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic  -o tutor_svm -lm 
clean:
	rm -f *~ svm.o tutor_svm libsvm.so.$(SHVER)