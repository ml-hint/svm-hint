CXX ?= g++
CC = icpc
SHVER = 2
OS = $(shell uname)

export OMP_ENABLE:=0

LIBSVM_DIR = ./libsvm-weights-3.20/

SVMHINT_SRC = ./src/
SVMHINT_LIB = ./lib/

ifndef ROOTSYS
$(error ROOTSYS is not defined!)
endif

ROOT_LIBS = `root-config --libs` -lGenVector -lMathMore -lMinuit2
ifeq ($(OMP_ENABLE),1)
	CXX_FLAGS =  -O3 -fPIC -fopenmp `root-config --cflags`
else
	CXX_FLAGS =  -O3 -fPIC `root-config --cflags`
endif


ifeq ("$(OS)","Darwin")
	CXX_FLAGS += -std=c++14
	INCS = -I/usr/include -I${ROOTSYS}/include -I$(LIBSVM_DIR) -I$(SVMHINT_LIB) -lz -lCore
else
	CXX_FLAGS += -std=c++0x
	INCS = -I/usr/include -I${ROOTSYS}/include -I$(LIBSVM_DIR) -I$(SVMHINT_LIB) -L/usr/lib64 -lz -lCore		
endif


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
	@echo "Compiling LIBSVM..." 	
	$(CXX) $(CXX_FLAGS) -c  $(LIBSVM_DIR)/svm.cpp -D"OMP_ENABLE=${OMP_ENABLE}"
	@echo "Done." 	
tutorial:$(SVMHINT_TEST)/svm_hint_tutor.cpp  $(SVMHINT_SRC)/fom.cpp $(SVMHINT_LIB)/fom.h $(SVMHINT_LIB)/libsvm_container.h  $(SVMHINT_LIB)/svm_hint.h $(SVMHINT_LIB)/csvc_interface.h svm.o 
	@echo "Compiling $@..." 	
	$(CXX) $(CXX_FLAGS) $(SVMHINT_TEST)/svm_hint_tutor.cpp svm.o $(SVMHINT_SRC)/fom.cpp  -D"OMP_ENABLE=${OMP_ENABLE}" $(SVMHINT_SRC)/user_csvc_interface.cpp  $(SVMHINT_SRC)/csvc_interface.cpp $(INCS) $(ROOT_LIBS) -o tutor_svm -lm  
	@echo "Done." 		
TTH:$(SVMHINT_TEST)/svm_TTH.cpp  $(SVMHINT_SRC)/fom.cpp $(SVMHINT_LIB)/fom.h $(SVMHINT_LIB)/libsvm_container.h  $(SVMHINT_LIB)/svm_hint.h $(SVMHINT_LIB)/csvc_interface.h svm.o 
	@echo "Compiling $@..." 	
	$(CXX) $(CXX_FLAGS) $(SVMHINT_TEST)/svm_TTH.cpp svm.o $(SVMHINT_SRC)/fom.cpp  -D"OMP_ENABLE=${OMP_ENABLE}" $(SVMHINT_SRC)/user_csvc_interface.cpp  $(SVMHINT_SRC)/csvc_interface.cpp $(INCS) $(ROOT_LIBS) -o TTH_svm -lm  
	@echo "Done." 		
asimov: $(SVMHINT_TEST)asimov.cpp $(SVMHINT_SRC)\fom.cpp $(SVMHINT_LIB)\fom.h 
	@echo "Compiling $@..." 	
	$(CXX) $(CXX_FLAGS) $(SVMHINT_TEST)\asimov.cpp $(SVMHINT_SRC)\fom.cpp -I$(LIBSVM_DIR) -I$(SVMHINT_LIB) $(INCS) $(ROOT_LIBS) -o asimov -lm	
	@echo "Done." 	
clean:
	rm -f *~ svm.o tutor_svm asimov libsvm.so.$(SHVER)