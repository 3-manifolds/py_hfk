// $Author$ $Date$ $Revision$
// C++ class to compute Heegaard Floer knot homology.
//
// This code is derived, with the permission of the authors,
// from code written by John A. Baldwin and William D. Gillam
// who retain all rights to the code.

#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <list>
using namespace std;

#if defined(__GNUC__)
typedef long long int Int64;
#elif defined(_MSC_VER)
typedef __int64 Int64;
#endif

// Class Declarations

class Vertex {
public:
  Int64 perm;
  list<int> out;
  list<int> in;
  bool alive;
  Vertex(Int64 p);
  ~Vertex();
};

class Generator {
public:
  int m;
  int a;
  Generator(int m, int a);
  ~Generator();
};

class Homology {
  friend class Link;
public:
  Homology();
  ~Homology();
private:
  list < Generator > gens;
  void add(Generator g);
  int maxM;
  int maxA;
};

class BiArray {
public:
  int M;
  int A;
  int& operator()(int m, int a);
  int operator()(int m, int a) const;
  BiArray(int M, int A);
  ~BiArray();
private:
  int *data;
};

class Link {
public:
  int gridsize;
  int *black;
  int *white;
  int Ashift;
  int HFK_maxA;
  int HFK_minM;
  int HFK_maxM;
  int aborted;
  int WindingNumber(int x, int y);
  int MaslovGrading(int g[]);
  int MOS_Rank(int m, int a);
  int HFK_Rank(int m, int a);
  Link(int gridsize, 
       int *black, 
       int *white, 
       int (*Progress)(const char *, int),
       bool quiet=1);
  ~Link();
private:
  int graphsize;
  vector < vector <Int64> > Generators;
  vector < vector <Vertex> > Graph;
  Homology MOS;
  unsigned char Rectangles[16*16*16*16];
  BiArray *MOS_Array;
  BiArray *HFK_Array;
  int (*Progress)(const char *message, int percent);
  bool quiet;
  short *counter;
  int *g;
  int *gij;
  int  AlexanderShift();
  int  RectDotFree(int xll, int yll, int xur, int yur, int which);
  int  BuildVertices();
  int  BuildEdges(int MM, const char *msg, int *count);
  int  Reduce(int MM, const char *msg, int *count);
  int  FindHomology();
  void BuildRectangles();
  void ComputeMOSRanks();
  void ComputeHFKRanks();
};

// Function Prototypes
void      Int2Perm(Int64 k, int h [], int size);
Int64     Perm2Int(int y [], int size);
void      NextPerm(short *counter, int *h, int size);
bool      ValidGrid(int gridsize, int black[], int white[]);
int       Find(vector <Vertex> &V, Int64 x);

// Global Data

#if defined(__GNUC__)
Int64 Factorial[16] = {
  1LL,1LL,2LL,6LL,24LL,120LL,720LL,5040LL,40320LL,362880LL,3628800LL,
  39916800LL,479001600LL,6227020800LL,87178291200LL,1307674368000LL};
#elif defined(_MSC_VER)
Int64 Factorial[16] = {
  1,1,2,6,24,120,720,5040,40320,362880,3628800,
  39916800,479001600,6227020800,87178291200,1307674368000};
#endif

#define Binomial(n,k) (Factorial[n]/(Factorial[k]*Factorial[n-k]))

// Class Methods
Vertex::Vertex(Int64 p): perm(p){
  alive = 1;
}

Vertex::~Vertex(){
}

Generator::Generator(int m, int a): m(m), a(a){
}

Generator::~Generator(){
}

Homology::Homology(){
  maxM = 0;
  maxA = 0;
}

Homology::~Homology(){};

void Homology::add(Generator g){
    if (abs(g.m) > maxM) maxM = abs(g.m);
    if (abs(g.a) > maxA) maxA = abs(g.a);
    gens.push_back(g);
}

BiArray::BiArray(int M, int A): M(M), A(A) {
    int i;
    data = new int[(2*M+1)*(2*A+1)];
    for (i=0; i<(2*M+1)*(2*A+1); i++) data[i]=0;
}

BiArray::~BiArray(){
  delete[] data;
}

inline
int& BiArray::operator()(int m, int a) {
  if ( (m < - M) || (m > M) || (a < -A) || (a > A) )
    cout <<"Out of bounds ("<< m <<","<< a <<") in ("<< M <<","<< A <<")\n";
    return data[(2*A+1)*(M+m) + (A+a)];
}

inline
int BiArray::operator()(int m, int a) const {
  if ( (m < - M) || (m > M) || (a < -A) || (a > A) )
    cout <<"Out of bounds ("<< m <<","<< a <<") in ("<< M <<","<< A <<")\n";
  return data[(2*A+1)*(M+m) + (A+a)];
}

Link::Link(int gridsize, 
	   int *black, 
	   int *white, 
	   int (*Progress)(const char *, int),
	   bool quiet) :
  gridsize(gridsize),
  black(black),
  white(white),
  Progress(Progress),
  quiet(quiet){
    MOS_Array = NULL;
    HFK_Array = NULL;
    Graph.resize(4*gridsize);
    Generators.resize(4*gridsize);
    Ashift = AlexanderShift();
    BuildRectangles();
    counter = new short [gridsize-1];
    g = new int [gridsize];
    gij = new int [gridsize];
    aborted = ( gridsize > 16 || BuildVertices() || FindHomology() );
    if (aborted )return;
    MOS_Array = new BiArray(MOS.maxM, MOS.maxA);
    ComputeMOSRanks();
    HFK_Array = new BiArray(MOS.maxM+gridsize, MOS.maxA+gridsize);
    ComputeHFKRanks();
  };

Link::~Link(){
  if (MOS_Array) delete MOS_Array;
  if (HFK_Array) delete HFK_Array;
  delete [] counter;
  delete [] g;
  delete [] gij;
};

// Compute winding number around (x,y)
int Link::WindingNumber(int x, int y)
{
  int i, ret=0;
  for(i=0; i<x; i++) {
    if ((black[i] >= y) && (white[i] < y)) ret++;
    if ((white[i] >= y) && (black[i] < y)) ret--;
  }
  return ret;
}

// Compute the Alexander grading shift.
int Link::AlexanderShift(){
  int i, shift, WN=0;
  for(i=0; i<gridsize; i++) {
    WN += WindingNumber(i, black[i]);
    WN += WindingNumber(i, (black[i]+1)%gridsize);
    WN += WindingNumber((i+1)%gridsize, black[i]);
    WN += WindingNumber((i+1)%gridsize, (black[i]+1)%gridsize);
  } 
  for(i=0; i<gridsize; i++) {
    WN += WindingNumber(i , white[i]);
    WN += WindingNumber(i ,(white[i]+1)%gridsize);
    WN += WindingNumber((i+1)%gridsize, white[i]);
    WN += WindingNumber((i+1)%gridsize, (white[i]+1)%gridsize);
  }
  /*
  if ( WN > 0 ) {
    int *temp;
    temp = black;
    black = white;
    white = temp;
    WN = -WN;
  }
  */
  shift = (WN - 4 * gridsize + 4)/8;
  if (!quiet)
    cout << "Alexander Grading Shift: " << shift << "\n";
  return shift;
}

// Returns 1 if the given rectangle has no dot; 0 if it has a dot;
// or -1 on error.
int Link::RectDotFree(int xll, int yll, int xur, int yur, int which){
  int dotfree = 1;
  switch (which) {
  case 0: 
    for(int x=xll; x<xur && dotfree; x++) {
      if (white[x] >= yll && white[x] < yur) dotfree = 0;
      if (black[x] >= yll && black[x] < yur) dotfree = 0;
    }
    return dotfree;
  case 1:
    for(int x=0; x<xll && dotfree; x++) {
      if (white[x] < yll || white[x] >= yur) dotfree = 0;
      if (black[x] < yll || black[x] >= yur) dotfree = 0;
    }
    for(int x=xur; x<gridsize && dotfree; x++) {
      if (white[x] < yll || white[x] >= yur) dotfree = 0;
      if (black[x] < yll || black[x] >= yur) dotfree = 0;
    }
    return dotfree;
  case 2:
    for(int x=xll; x<xur && dotfree; x++) {
      if (white[x] < yll || white[x] >= yur) dotfree = 0;
      if (black[x] < yll || black[x] >= yur) dotfree = 0;
    }
    return dotfree;
  case 3:
    for(int x=0; x<xll && dotfree; x++) {
      if (white[x] >= yll && white[x] < yur) dotfree = 0;
      if (black[x] >= yll && black[x] < yur) dotfree = 0;
    }
    for(int x=xur; x<gridsize && dotfree; x++) {
      if (white[x] >= yll && white[x] < yur) dotfree = 0;
      if (black[x] >= yll && black[x] < yur) dotfree = 0;
    }
    return dotfree;
  }
  return -1; //Error!
}

 // Compute Maslov grading via the formula:
 // 4M(y)
 // = 4M(white)+4P_y(R_{y, white})+4P_{white}(R_{y, white})-8W(R_{y, white})
 // = 4-4*gridsize+4P_y(R_{y, white})+4P_{x_0}(R_{y, white})-8W[j](R_{y, white})

int Link::MaslovGrading(int g[]) {
  int i, j, k, wi, Wi, Wj, gi, Gi, Gj;
  int P = 4 - 4*gridsize; // Four times the Maslov grading of x_0
  for(i=0; i<gridsize; i++) {
    Wi = white[i];
    wi = Wi - 1;
    if (wi < 0) wi += gridsize;
    Gi = g[i];
    gi = Gi - 1;
    if (gi < 0) gi += gridsize;

    // Calculate incidence number R_{y x_0}.S for each of the four
    // squares S having (i,white[i]) as a corner and each of the 
    // four squares having (i,g[i]) as a corner and shift P appropriately


    for(j=0; j<=i; j++) {
      Wj = white[j];
      Gj = g[j];
      // Squares whose BL corners are (i,white[i]) and (i,g[i])
      // multiply by 7 because of the -8W(R_{y, white}) contribution
      P -= 7*((Wj > Wi) & (Gj <= Wi));
      P += 7*((Gj > Wi) & (Wj <= Wi)); 
      P +=   ((Wj > Gi) & (Gj <= Gi));
      P -=   ((Gj > Gi) & (Wj <= Gi)); 
      // Squares whose TL corners are (i,white[i]) and (i,g[i]) (mod gridsize)
      P += ((Wj > wi) & (Gj <= wi));
      P -= ((Gj > wi) & (Wj <= wi)); 
      P += ((Wj > gi) & (Gj <= gi));
      P -= ((Gj > gi) & (Wj <= gi)); 
    }
    k = i-1;
    if (k < 0) k += gridsize;
    // Squares whose BR corners are ...
    for(j=0; j<= k; j++) {
      Wj = white[j];
      Gj = g[j];
      P += ((Wj > Wi) & (Gj <= Wi));
      P -= ((Gj > Wi) & (Wj <= Wi)); 
      P += ((Wj > Gi) & (Gj <= Gi));
      P -= ((Gj > Gi) & (Wj <= Gi)); 
    // Squares whose TR corners are ...
      P += ((Wj > wi) & (Gj <= wi));
      P -= ((Gj > wi) & (Wj <= wi)); 
      P += ((Wj > gi) & (Gj <= gi));
      P -= ((Gj > gi) & (Wj <= gi)); 
    }
  }
  return (P/4);
}

// Add a vertex for each permutation with non-negative Alexander
// grading.
int Link::BuildVertices() {
  int i, AGrading, percent;
  Int64 count = 0, end = 0, step;
  stringstream msg;
  int winding_numbers[16][16];

  msg << "Constructing generators ... ";

  // Build a table of winding numbers.
  for(int x=0; x < gridsize; x++) {
    for(int y=0; y < gridsize; y++) {
      winding_numbers[x][y] = WindingNumber(x, y);
    }
  }
  if (!quiet) {
    cout << "Matrix of winding numbers:\n";
    for(int y=gridsize - 1; y>=0; y--) {
      for(int x=0; x<gridsize; x++) {
	cout << setw(2) << winding_numbers[x][y];
      }
      cout << "\n";
    }
  }

  for(i=0; i < gridsize-1; i++)
    counter[i] = 0;
  
  step = 100000 < Factorial[gridsize] ? 100000 : Factorial[gridsize];
  Progress(msg.str().c_str(), 0);
  while (end < Factorial[gridsize]) {
    end = count + step;
    if (end > Factorial[gridsize])
      end = Factorial[gridsize];
    percent = (int)(100.0*(count/(float)Factorial[gridsize]));
    if ( Progress(msg.str().c_str(), max(percent,1)) )
      return -1;
     //This loop accounts for most of the computation.
    for(; count < end; count++) {
      NextPerm(counter, g, gridsize);
      AGrading = Ashift;
      for(i=0; i<gridsize; i++)
	AGrading -= winding_numbers[i][g[i]];
      if (AGrading >= 0) {
	Generators[gridsize + gridsize + MaslovGrading(g)].push_back(count);
      }
    }
  }
  Progress(msg.str().c_str(), 100);
  msg.str("");
  graphsize = 0;
  for(i=0; i < 4*gridsize; i++)
    graphsize += Generators[i].size();
  msg << "Number of generators: " << graphsize << "\n";
  Progress(msg.str().c_str(), -1);
  return 0;
}

void Link::BuildRectangles(){
  // Build a table showing whether a rectangle on the torus contains a dot.
  for(int xll=0; xll < gridsize; xll++) {
    for(int xur=xll+1; xur < gridsize; xur++) {
      for(int yll=0; yll < gridsize; yll++) {
	for(int yur=yll+1; yur < gridsize; yur++) {
	  Rectangles[xll<<12 | yll<<8 | xur<<4 | yur] = 
	    RectDotFree(xll,yll,xur,yur,0)    |
	    RectDotFree(xll,yll,xur,yur,1)<<1 |
	    RectDotFree(xll,yll,xur,yur,2)<<2 |
	    RectDotFree(xll,yll,xur,yur,3)<<3;
	}
      }
    }
  }
}

// Add edges describing the boundary map.
int Link::BuildEdges(int MM, const char *msg, int *count) {
  int i, j, k, Gi, Gj, indexgij, rectinfo; 
  int index = 0, edges = 0, end = 0, step, percent;
  bool firstrect = 0;
  bool secondrect = 0;

  step = min(20000, (int)Graph[MM].size());

  // Add the edges.
  while (end < (int)Graph[MM].size()) {
    end = min(index + step, (int)Graph[MM].size());
<<<<<<< local
    percent = (*count + index)/(1 + graphsize/50);
=======
    percent = int(50.0*(*count + index)/(float)graphsize);
>>>>>>> other
    if ( Progress(msg, max(1, percent)) )
      return -1;
    for(; index < end; index++) {
      Int2Perm(Graph[MM][index].perm, g, gridsize);
      for(i=0; i<gridsize; i++) {
	Gi = g[i];
	for(j=i+1; j<gridsize; j++) {
	  Gj = g[j];
	  if(Gi < Gj) {
	    rectinfo = Rectangles[i<<12 | Gi<<8 | j<<4 | Gj];
	    firstrect = (bool)(rectinfo & 1);
	    for(k=i+1; k<j && firstrect; k++) {
	      firstrect = !(bool)((Gi < g[k]) & (g[k] < Gj));
	    }
	    secondrect = (bool)(rectinfo & 2);
	    for(k=0; k<i && secondrect; k++) {
	      secondrect = !(bool)((g[k] < Gi) | (g[k] > Gj));
	    }
	    for(k=j+1; k<gridsize && secondrect; k++) {
	      secondrect = !(bool)((g[k] < Gi) | (g[k] > Gj));
	    }
	  }
	  if(Gj < Gi) {
	    rectinfo = Rectangles[i<<12 | Gj<<8 | j<<4 | Gi];
	    firstrect = rectinfo & 4;
	    for(k=i+1; k<j && firstrect; k++) {
	      firstrect = !((g[k] < Gj) | (g[k] > Gi));
	    }
	    secondrect = rectinfo & 8;
	    for(k=0; k<i && secondrect; k++) {
	      secondrect = !((g[k] > Gj) & (g[k] < Gi));
	    }
	    for(k=j+1; k<gridsize && secondrect; k++) {
	      secondrect = !((g[k] > Gj) & (g[k] < Gi));
	    }
	  }
	  if(firstrect != secondrect) { // Exactly one rectangle is a boundary
	    for(k=0; k<gridsize; k++) {
	      gij[k] = g[k];
	    }
	    gij[i] = Gj;
	    gij[j] = Gi;
	    indexgij = Find(Graph[MM-1], Perm2Int(gij, gridsize));
	    if(indexgij==-1) {
	      Progress("Error with Alexander grading!!", -1);
	      return -1;
	    }
	    Graph[MM][index].out.push_back( indexgij );
	    Graph[MM-1][indexgij].in.push_back( index );     
	    edges++;
	  }
	}
      }
    }
  }
  *count += end;
  return 0;
}
  
// Eliminate edges to find a homology basis.
int Link::Reduce(int MM, const char *msg, int *count){
  int source = 0, target;
  int end = 0, step, percent;

  step = min(20000, (int)Graph[MM].size());

  while (end < (int)Graph[MM].size()) {
    end = min(source + step, (int) Graph[MM].size());
<<<<<<< local
    percent = (*count + source)/(1 + graphsize/50);
=======
    percent = (int)50.0*(*count + source)/((float)graphsize);
>>>>>>> other
    if ( Progress(msg, max(1,percent)) )
      return -1;
    for(; source<end; source++) {
      if( (!Graph[MM][source].alive) || Graph[MM][source].out.size()==0)
	continue;

      target = Graph[MM][source].out.front();
      Graph[MM][source].alive = 0;
      
      // For every j->target, j != source
      for(list<int>::iterator j = Graph[MM-1][target].in.begin();
	  j != Graph[MM-1][target].in.end(); j++) {
	    if( !Graph[MM][*j].alive ) continue;
	    // For every source->k
	    for(list<int>::iterator k = Graph[MM][source].out.begin();
		k != Graph[MM][source].out.end(); k++) {
	      // Search for j->k. If found, remove it; if not, add it.
	      list<int>::iterator search;
	      for(search = Graph[MM][*j].out.begin();
		  search != Graph[MM][*j].out.end(); search++ )
		if (*search == *k) {
		  Graph[MM][*j].out.erase( search );
		  if( *k != target) Graph[MM-1][*k].in.remove(*j);
		  break;
		}
	      if (search == Graph[MM][*j].out.end()){
	        Graph[MM][*j].out.push_back(*k);
	        Graph[MM-1][*k].in.push_back(*j);
	      } 
	    }
	  }
      
      // For each source->j, remove source from j.in
      for(list<int>::iterator j = Graph[MM][source].out.begin();
	  j != Graph[MM][source].out.end(); j++)
	Graph[MM-1][*j].in.remove(source);
      
      //Kill source and target
      Graph[MM-1][target].alive = 0;
      Graph[MM-1][target].out.clear();
      Graph[MM-1][target].in.clear();
      Graph[MM][source].out.clear();
      Graph[MM][source].in.clear();
    }
  }
  *count += end;
  return 0;
}

int Link::FindHomology() {
  int i, j, count=0, MM;
  const char *msg="Computing homology ...";

  Progress(msg, 0);

  MM = 4*gridsize;
  for(i=0; i < (int)Generators[MM-1].size(); i++)
    Graph[MM-1].push_back(Vertex(Generators[MM-1][i])); 


  for (MM = 4*gridsize -1; MM > 0; MM--) {
    for(i=0; i < (int)Generators[MM-1].size(); i++) {
      Graph[MM-1].push_back(Vertex( Generators[MM-1][i])); 
    }
    Generators[MM-1].clear();

    aborted = BuildEdges(MM, msg, &count) || Reduce(MM, msg, &count);
    if ( aborted ) return aborted;

    for(i=0; i < (int)Graph[MM].size(); i++) {
      if(Graph[MM][i].alive) {
	Int2Perm(Graph[MM][i].perm, g, gridsize);
	int AGrading = Ashift;
	for(j=0; j<gridsize; j++)
	  AGrading -= WindingNumber(j, g[j]);
	MOS.add(Generator(MaslovGrading(g),AGrading));
      }
    }

    Graph[MM].clear();
  }

  Progress(msg, 100);
  Progress("", -1);
  return 0;
}

// Compute MOS ranks.
void Link::ComputeMOSRanks(){
  for (list<Generator>::iterator 
	 G=MOS.gens.begin(); G != MOS.gens.end(); G++)
    (*MOS_Array)((*G).m, (*G).a)++;
  if (!quiet) {
    cout << "MOS ranks for non-negative Alexander grading:\n";
    for(int a = MOS.maxA; a >= 0; a--) {
      for(int m = -MOS.maxM; m <= MOS.maxM; m++) {
	cout << setw(5) << (*MOS_Array)(m,a);
      }
      cout << "\n";
    }
  }
}

// Compute HFK^ ranks from MOS ranks.
void Link::ComputeHFKRanks(){
  int a, m, i;

  for(a = 0; a <= MOS.maxA; a++)
    for(m = -MOS.maxM; m <= MOS.maxM; m++)
      (*HFK_Array)(m,a) = MOS_Rank(m,a);

  // MOS =  \hat HFK \otimes K^{gridsize-1}
  for(a = MOS.maxA + gridsize ; a >= -MOS.maxA - gridsize; a--) {
    for(m = MOS.maxM + gridsize; m >= -MOS.maxM - gridsize; m--) {
      if( (*HFK_Array)(m,a) > 0) {
	for(i=1; i <= gridsize-1; i++)
	  (*HFK_Array)(m-i,a-i) -= (*HFK_Array)(m,a)*Binomial(gridsize-1,i);
      }
    }
  }

  // Symmetrize, since MOS was only computed for a >= 0. 
  for(a = 1; a <= MOS.maxA; a++)
    for(m = -MOS.maxM - gridsize; m <= MOS.maxM + gridsize - 2*MOS.maxA; m++)
      (*HFK_Array)(m,-a) = (*HFK_Array)(m + 2*a, a);

    HFK_maxA = MOS.maxA;
    HFK_minM = 0;
    HFK_maxM = 0;
    for (m = -MOS.maxM - gridsize; m <=0; m++) {
      for (a = -MOS.maxA; a <= MOS.maxA; a++)
	if ( (*HFK_Array)(m,a) ) {
	  HFK_minM = m;
	  break;
	}
      if (HFK_minM) break;
    }
    for (m = MOS.maxM + gridsize; m >=0; m--) {
      for (a = MOS.maxA; a >= MOS.maxA; a--)
	if ( (*HFK_Array)(m,a) ) {
	  HFK_maxM = m;
	  break;
	}
      if (HFK_maxM ) break;
    }
}

int Link::MOS_Rank(int m, int a) {
  return (*MOS_Array)(m,a);
}

int Link::HFK_Rank(int m, int a) {
  return (*HFK_Array)(m,a);
}

// Returns i if V[i]=x or -1 if x isn't a member of V.
int Find(vector <Vertex> &V, Int64 x) {
  int above=V.size()-1;
  int below=0;
  while(above - below > 1) {
    if (x >= V[below+(above-below)/2].perm) below += (above-below)/2;
    else above = below+(above-below)/2;
  }
  if (V[below].perm == x) return below;
  if (V[above].perm == x) return above;
  return -1;
}

// Check that we have exactly 2 dots of each color in each row and column.
bool ValidGrid(int gridsize, int black[], int white[]) {
 int numwhite=0;
 int numblack=0;
 for(int i=0; i<gridsize; i++) {
  for(int j=0; j<gridsize; j++) {
   if (white[j]==i) numwhite++;
   if (black[j]==i) numblack++;
  }
  if (numwhite != 1 || numblack != 1)
   return 0;
  numwhite=0;
  numblack=0;
 }
 return 1;
}

// Maps a permutation of size n to an integer < n!
// See: Knuth, Volume 2, Section 3.3.2, Algorithm P
Int64 Perm2Int(int Perm[], int size) {
 Int64 Int=0;
 int temp, m, r = size;
 while (r > 0){
   for (m=0; m < r; m++){
     if (r - Perm[m] == 1)
       break;
   }
   Int = Int*r + m;
   r -= 1;
   temp = Perm[r];
   Perm[r] = Perm[m];
   Perm[m] = temp;
 }
 return Int;
}

// Inverse mapping, from integers < n! to permutations of size n
// Writes the permutation corresponding to Int into the array Perm.
void Int2Perm(Int64 Int, int Perm[], int size) {
  int r, m, temp;
  for(int i=0; i<size; i++) Perm[i]=i;
  r = 1;
  while (r < size)  {
    m = Int%(r+1);
    Int = Int/(r+1);
    temp = Perm[r];
    Perm[r] = Perm[m];
    Perm[m] = temp;
    r += 1;
  }
  return;
}

// Generator for permutations.  Inputs a factorial based counter and
// an array.  Writes the permutation indexed by the counter into the
// array, and then increments the counter.
void NextPerm(short counter[], int P[], int size) {
  int i, m, temp;
  for(i=0; i<size; i++) P[i]=i;
  for (i=1; i < size; i++) {
    m = counter[i-1];
    temp = P[i];
    P[i] = P[m];
    P[m] = temp;
  }
  for (i=0; i<size-1; i++) {
    counter[i] += 1;
    if (counter[i] == i+2)
      counter[i] = 0;
    else
      break;
  }
  return;
}
