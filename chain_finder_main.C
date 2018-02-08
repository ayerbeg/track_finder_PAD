#include <iostream>
#include "string.h"
#include "string"
#include "TFile.h"
#include "TTree.h"
#include "math.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include <TBenchmark.h>
#include <fstream>
#include <sstream>

#include "chain_finder_main.h"

#include "TVector3.h"

/* The status bits for hits in the hitlists */
// I keep these definitions from the original code from Howard
// they were much more, but it was defined in the rtpc.h code
// 

#define HISUSED 1  /* - (was 1) used */

#define HITUNAV 1 /* reasons to NOT use a hit in the chain-linker */

using namespace std;


void search(HitStruct hitevent[], Int_t max_hits, Int_t event_ii )
{
  rad2deg = 180/(4*atan(1));

  // *************************************
  // THESE VARIABLES WERE READ FROM A FILE
  
  Max_Link_Sep = space;
  Max_Ang      = max_ang;
  Min_Ang      = min_ang;
  Ang_Sep      = ang_sep;
  
  // *************************************
  
  TVector3 nil_vec(0,0,0);//null vector

  for (anchor_hit = 0; anchor_hit < max_hits; anchor_hit++) //index of anchor hit
    {
      if (event_ii==0)	  cout<<hitevent[anchor_hit].Status<<" "<<anchor_hit<<endl;
	  
      if ( hitevent[anchor_hit].Status & HITUNAV )
	//If Status=1 (used), do not do anything. 
	{
	  //
	}
	 
      else
	{
	  acceptance = 0;
	  num_hits_this_chain[num_chains] = 1;//THERE IS ALWAYS ONE HIT IN A CHAIN
	      
	  chain_hits[num_chains][0] = anchor_hit;// **(A)**
	      
	  hitevent[anchor_hit].Status |= HISUSED;// |= --> or eq;  a = a | b. It just assign the USED status to the hit

	  if (event_ii==0) cout<< hitevent[anchor_hit].Status<<endl;

	  //SEARCH ALGORITHM----->
	      
	  for (seed_hit = 0; seed_hit < num_hits_this_chain[num_chains]; seed_hit++)
	    {
	      if (event_ii==0) cout<<" HERE 1: "<<hitevent[seed_hit].X<<" "<< hitevent[seed_hit].Y<<" "<< hitevent[seed_hit].Z<<" seed_hit: "<<seed_hit<<endl;
		  
	      seed_index = chain_hits[num_chains][seed_hit];// **(A)** (it is related to the anchor hit)
	
	      TVector3 pseed(hitevent[seed_index].X, hitevent[seed_index].Y, hitevent[seed_index].Z);

	      for(next_hit = 0; next_hit<max_hits; next_hit++)
		//max_hits: number of hits inside the event
		//NOTE: we have to take care about the variables, where are defined, how are handled

		{
		  if( !(hitevent[next_hit].Status & HITUNAV) )
		    {
		      TVector3 pnext(hitevent[next_hit].X, hitevent[next_hit].Y, hitevent[next_hit].Z);
			  
		      if (next_hit>0 && (pnext!=nil_vec))
			{
			  TVector3 pnext_pre(hitevent[next_hit-1].X, hitevent[next_hit-1].Y, hitevent[next_hit-1].Z);
			  if (!(pnext-pnext_pre).Mag() == 0) // definitely removes the same hits (perhaps redundant)
			    {
			      TVector3 dif_vec(pnext - pseed);
				  
			      separation = (dif_vec).Mag();//REMEMBER!! They can be treated as vectors for the calculus we use
				  
			      //  we differenciate the separation in the first to the second
			      if(
				 ( separation <= Max_Link_Sep && separation > 0.0 )
				 )
				//by definition, separation is always>0 but with this condition,
				//we remove 0's (i.e. repeated hits) IT MUST BE SOLVED WITH THE PREVIOUS CONDITION
				{
				  index_hits = num_hits_this_chain[num_chains];//Just rename the variable to handle better
				      
				  //gra_vec[0] is always 0,0,0
				  gra_vec[index_hits].SetX(dif_vec.X());
				  gra_vec[index_hits].SetY(dif_vec.Y());
				  gra_vec[index_hits].SetZ(dif_vec.Z());
				      
				  if(index_hits>1)
				    {
				      acceptance = dif_vec.Angle(gra_vec[index_hits-1]) *rad2deg;
					  
				    }//if(index_hits>1...	  
				      
				  if(acceptance>90.) acceptance = 180. - acceptance;
				      
				  if (separation <= Ang_Sep && acceptance < Max_Ang-15)
				    {
				      accept_hit(next_hit);
				    }
				      
				  if (separation > Ang_Sep && acceptance < Min_Ang-10)
				    {					
				      accept_hit(next_hit);
				    }
				      
				  //					}
				      
				}// if(separation<=
				  
			      else
				{}
				  
			    }//if(pnext
			}//if(next_hit>0			  
		    }// if( !(hitevent
		      
		}// for(next_hit
		  
	    }//  for (seed_hit

	  //**************************************************************
	  // The conditions breaks when there are no more hits in the pool
	  // susceptible to be part of the chain
	  //**************************************************************
	  
	  if( num_hits_this_chain[num_chains] >= 6) 
	    {

	      store_data(event_ii);// sending the index of event (ii) is just to control how many times the event was studied, i.e. the split tracks. it can be rid out. 
	
	    }
	  else
	    {
	      hitevent[next_hit].Status &= ~HISUSED;
	    }
	      
	}//else (from if( hitevent[anchor_hit].Status & HITUNAV)
	  
    }//for(anchor_hit...


  
  //***********************************************
  // THESE LINES ARE USING FOR EFFICIENCY/DEBUGGING
  // COULD BE DELETED

  cout<<"Total Events: "<<real_eve_counter<<endl;//The number of first time entering into the search
  cout<<"double counts (split tracks): " <<double_counter<<endl;//Although is called double, could be more times analyzing the same event, i.e. split tracks
  cout<<"more than double counts (split tracks): " <<mcts<<endl;//how many tracks are splited in more than two
  //***********************************************

  
  
}
  

void accept_hit(Int_t next_hit)
{
  	  
  if (num_hits_this_chain[num_chains] >= MAX_HITS_ON_CHAIN)
    {        
      printf("Too many hits for the chain list. Aborting.\n"); 
      return;//It was -1 but is a void function
    }
    
  //THE HIT INDEX WHICH ACOMPLISH THE CONDITION, IS STORED-->
  chain_hits[num_chains][num_hits_this_chain[num_chains]] = next_hit;
  
  /* mark it as used */
  hitevent[next_hit].Status |= HISUSED; //this kind of assignment is for 'historical' reasons
  num_hits_this_chain[num_chains]++;//ADD HIT TO THE CHAIN, INCREASE INDEX

  //Its final number (when finished) is the total number of hits in the chain
  
}





void store_data(Int_t ii)
{
  //I just move the algorithm to a function. SHOULD BE CLEANED. I think I don't need to send any argument as
  // is written now. 
  
	      
  printf("....KEEPING THIS CHAIN #%d. %d hits \n", num_chains,num_hits_this_chain[num_chains]);
  
  //***********************************************
  // THESE LINES ARE USING FOR EFFICIENCY/DEBUGGING
  // COULD BE DELETED
  
  if(ii!=temp_eve)
    {
      real_eve_counter++;
      cout<<"EVENT: "<<ii<<endl;
      mcts = 0;
    }
  else
    {
      cout<<"DOUBLE!!! <<<<<<<<------------------------------------------"<< ii<<endl;
      double_counter++;
      
      if(mcts==1)mcts++;
      if(mcts==0)mcts = 1;
    }
  //***********************************************

  
  temp_eve = ii;
  
  
  for (Int_t jj=0; jj<num_hits_this_chain[num_chains]; jj++)
    {
      printf(" %d", chain_hits[num_chains][jj]);
      
      //Fill the variables value of the found chain
      ChainArray[num_chains].X_rec[jj] = hitevent[chain_hits[num_chains][jj]].X;
      ChainArray[num_chains].Y_rec[jj] = hitevent[chain_hits[num_chains][jj]].Y;
      ChainArray[num_chains].Z_rec[jj] = hitevent[chain_hits[num_chains][jj]].Z;
      
    }

//NOTE: THERE WAS A CLEAN-UP LOOP, BUT WAS REMOVED DUE TO THE ARRAY
  
  ChainArray[num_chains].Hit = num_hits_this_chain[num_chains];//number of hits in the chain
  ChainArray[num_chains].ID = num_chains;//Event Index. A Chain per index
  
  printf("\n\n");

  num_chains++; //NEW CHAIN, INCREASE THE INDEX

  // THE FINAL RESULT HERE IS A STRUCTURE WITH:
  // X,Y,Z=space coordinates, ID=chain ID, HIT=number of hits in the chain
 
  // BUT THIS IS THE STORE OF THE STRUCTURE, THE CODE IS STILL RUNNING
  // UNTIL THE LOOP (OVER HITS AND EVENTS) FINISH
  
}




