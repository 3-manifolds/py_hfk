// Sample application using hfk_class.cpp.
//
// This code is derived, with the permission of the authors,
// from code written by John A. Baldwin and William D. Gillam,
// who retain all rights to the code.

#include "hfk_class.cpp"

int Progress(const char *message, int percent){
  switch(percent) {
  case 0:
    cout << message << "       ";
    break;
  case -1:
    cout << "\n" << message << "\n";
    break;
  default:
    cout << "\b\b\b\b\b\b\b" << setw(4) << percent << "%  " << flush;
  }
  return 0;
}


int main(int argc, char *argv[]){
  int m, a;
  //12n_679
  //const int gridsize = 13;
  //int black[13] = {5,3,1,2,0,12,11,8,10,9,7,6,4};
  //int white[13] = {9,6,5,4,3,1,0,2,12,11,10,8,7};

  //Conway Mutant C_{2,1}
  const int gridsize = 11;
  int black[11]={10,9,3,4,5,8,6,7,1,2,0};
  int white[11]={6,1,7,0,3,10,9,2,4,8,5};

  //const int gridsize=12;
  //int black[12] = {4,1,3,7,8,9,11,2,10,6,5,0};
  //int white[12] = {11,5,0,10,1,7,8,4,3,2,9,6};

 // Check that the grid is valid.
  if(!ValidGrid(gridsize, black, white)) {
    cout << "Invalid grid!!\n";
    return 0;
  }

  Link link(gridsize, black, white, &Progress);
  if ( link.aborted ) {
    cout << "\n*****   Computation aborted!  *****\n";
    return -1;
  }

  cout << "\nHFK^ ranks:\n";
  for(a = link.HFK_maxA; a >= -link.HFK_maxA; a--) {
    for(m = link.HFK_minM; m <= link.HFK_maxM; m++)
      cout << setw(3) << link.HFK_Rank(m,a);
    cout << "\n";
  }
  cout << "\n";
  return 0;
}
