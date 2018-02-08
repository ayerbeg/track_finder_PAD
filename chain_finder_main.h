#ifndef chain_finder_main_h
#define chain_finder_main_h 1

#define ndim 500//so big? 
#define MAX_NUM_CHAINS 100 //in the original is 100, should studied
#define MAX_HITS_ON_CHAIN 500 //in the original is 100, should studied
#define MAX_LINK_SEP 11 //This is the maximum separation to be included in the chain

#include "TVector3.h"


typedef struct {
  Double_t X;
  Double_t Y;
  Double_t Z;
  Int_t Status;
} HitStruct;


typedef struct {
  Int_t ID; //This is the chain ID
  Double_t X_rec[MAX_HITS_ON_CHAIN];
  Double_t Y_rec[MAX_HITS_ON_CHAIN];
  Double_t Z_rec[MAX_HITS_ON_CHAIN];
  Int_t Hit;
} ChainStruct;

HitStruct hitevent[ndim];//The structure has dimensions of how many hits in the event

ChainEvStruct ChainArray[MAX_NUM_CHAINS];

Int_t    EventID;

Int_t    NoH;//Number of Hits

Int_t    anchor_hit, seed_hit, next_hit, seed_index;
Int_t    num_chains;

Int_t    num_hits_this_chain[MAX_NUM_CHAINS];
Int_t    chain_hits[MAX_NUM_CHAINS][MAX_HITS_ON_CHAIN];

Double_t separation, acceptance;

TVector3 gra_vec[MAX_HITS_ON_CHAIN];
TVector3 pnext_pre;
TVector3 pnext;

Int_t    index_hits;
Int_t    maxin; //max number of hits in an event. Calculated during the readout of the event


Int_t Entries;

double space;
double max_ang;
double min_ang;
double ang_sep;

Int_t temp_eve;
Int_t real_eve_counter;
Int_t double_counter;
Int_t mcts;

Double_t Max_Link_Sep;
Double_t Max_Ang;
Double_t Min_Ang;
Double_t Ang_Sep;

Double_t rad2deg;

void accept_hit(Int_t);
void search(HitStruct[], Int_t, Int_t);
void store_data(Int_t);

#endif 
