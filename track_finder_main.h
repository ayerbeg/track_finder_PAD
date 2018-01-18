#ifndef track_finder_main_h
#define track_finder_main_h 1

#define ndim 5000//so big? 
#define MAX_NUM_CHAINS 48000 //in the original is 100, should studied
#define MAX_HITS_ON_CHAIN 5500 //in the original is 100, should studied
#define MAX_LINK_SEP 11 //This is the maximum separation to be included in the chain

#include "TVector3.h"
#include <vector>
#include <map>


using namespace std;
typedef struct {
  Int_t Status;
  Int_t PAD;
  Double_t X;
  Double_t Y;
  Double_t Z;
  Int_t TID;
  Int_t pad;
 unsigned int tdc;
  Double_t adc;
  Int_t timeslot;
  Int_t Index;
} HitStruct;

HitStruct hitevent[ndim];//The structure has dimensions of how many hits in the event
HitStruct readevent[ndim];//it is declared here, but should be local defined because it is used in the readout


//************************
//INPUT ROOTFILE VARIABLES
TTree   *RTPCTree;
Int_t    event;
  vector<int>* Vpad;
  vector<int>* Vtdc;
  vector<double>* Vadc;
Int_t    nHits;
//************************

//************************
// HIT VECTOR
vector< vector<int> > vHit;//vector with the size of time slots with vectors containing the pads firing on each time slot
vector<int> vPAD;
//************************

vector< int > vChain;

vector<int> pad;


//************************
// DAVID'S VARIABLES
int maxconcpads = 0;
int concpads = 0;
int maxconctime = 0;

int Pad;

int PadPhi;
int PadZ;

unsigned int padindexmax = 0;

int TrigWindSize = 120;

int StepSize = 1;

int TID;
bool breakTIDloop = false;

double checkpadphi = 0;
double checkpadphiprev = 0;
double checkpadz = 0;
double checkpadzprev = 0;


Double_t Num_of_Rows = 180.;
int adjthresh = 2;

vector<int> TIDVec;
map<int, map<int,vector<int>> > TIDMap;
//************************


int tid_1;
int tid_2;
int tid_3;


Int_t    NoH;//Number of Hits

Int_t t = 0;

int next_hit_index_loop = 1;

std::vector<int> time_scale_array ;


vector<int> fTimePad;

Int_t    max_entries;//temporary variable to control number of entries to read

Int_t    adj_chain_hits[MAX_NUM_CHAINS][MAX_HITS_ON_CHAIN];

Int_t Entries;

//Int_t max_entries;



//************************************************
//PAD SYSTEM
Double_t ncol = 100.;
Double_t nrow = 180.;

int fRow(Int_t );
int fCol(Int_t );

Int_t rowT;
Int_t colT;

Int_t row[ndim];
Int_t col[ndim];

//************************************************


Int_t time_slot;


TString inputfile;
TString outputfile;


Int_t mcts;

TFile *rootoutfile;

TTree *chaintree ;

void chain_finder();

void store_data();

void variable(TString);

//void readout(Int_t, Int_t);
void search(Int_t, Int_t);
void readout(HitStruct[], Int_t);

void bubbleSort(HitStruct[] , int);
void time_label(HitStruct[] , int);
void swap(int , int );
  
#endif 
