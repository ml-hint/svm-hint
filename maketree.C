void maketree()
{
  //create a Tree file tree1.root

  //create the file, the Tree and a few branches
  TFile f("svm_testtrain.root","recreate");

 
  TTree t1("tree","example input tree for the SVM-interface tutorial Train and Test");
  Double_t var1 = 0., var2 = 0., var3 = 0., var4 = 0., weight = 0.; int sep = 0;
  Double_t random;
  TRandom3 rand;
  t1.Branch("var1",&var1,"var1/D");
  t1.Branch("var2",&var2,"var2/D");
  t1.Branch("var3",&var3,"var3/D");
  t1.Branch("var4",&var4,"var4/D");
  t1.Branch("weight",&weight,"weight/D");
  t1.Branch("sep",&sep,"sep/I");

  //fill the tree
  for (Int_t i=0;i<20000;i++) {
    if(i < 17500){
      var1 = (float)sin(rand.Gaus(60.,2.)); //first var
      var2 = (float)rand.Exp(10. ); // second var 
      var3 = (float)rand.Gaus(11.,8.); // third var
      var4 = sqrt((float)rand.Exp(1. )); 
      sep = -1;
      weight = rand.Gaus(0.5,0.000001);
    } else {
      var1 = (float)sin(rand.Gaus(90.,2.)); //first var
      var2 = (float)rand.Exp(9.)*rand.Gaus(2.,0.2); // second var 
      var3 = (float)rand.Gaus(11.,4.); // third var 
      var4 = sqrt((float)rand.Exp(5. )); 
      sep = 1;
      weight = rand.Gaus(0.01,0.000001);
    }
    
    t1.Fill();
  }
  t1.Write();
  TFile f2("svm_eval.root","recreate");
  TTree t2("tree","example input tree for the SVM-interface tutorial Evaluation");

  t2.Branch("var1",&var1,"var1/D");
  t2.Branch("var2",&var2,"var2/D");
  t2.Branch("var3",&var3,"var3/D");
  t2.Branch("var4",&var4,"var4/D");
  t2.Branch("weight",&weight,"weight/D");
  t2.Branch("sep",&sep,"sep/I");

  for (Int_t i=0;i<20000;i++) {
    if(i < 17500){
      var1 = (float)sin(rand.Gaus(60.,2.)); //first var
      var2 = (float)rand.Exp(10. ); // second var 
      var3 = (float)rand.Gaus(11.,8.); // third var 
      var4 = sqrt((float)rand.Exp(1. )); 
      sep = -1;
      weight = rand.Gaus(0.5,0.000001);
    } else {
      var1 = (float)sin(rand.Gaus(90.,2.)); //first var
      var2 = (float)rand.Exp(9.)*rand.Gaus(2.,0.2); // second var 
      var3 = (float)rand.Gaus(11.,4.); // third var 
      var4 = sqrt((float)rand.Exp(5. )); 
      sep = 1;
      weight = rand.Gaus(0.01,0.000001);
    }
    
    t2.Fill();
  }
  t2.Write();
  //  TFile f("svm_testtrain.root","recreate");
}
