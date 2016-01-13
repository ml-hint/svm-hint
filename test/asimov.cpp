#include "fom.h"
#include <iostream>
using namespace std;
int main(int argc, char* argv[])
{

  if ( argc != 4) {
    cout << " Usage ./asimov signal background uncertainty" << endl; 
    return 0;
  }
  fom FOM(atof(argv[3]));
  FOM.setSignal(atof(argv[1]));
  FOM.setBackground(atof(argv[2]));
  cout << FOM.getSignificance(fom::asimov) << endl;

  return 0;
}
