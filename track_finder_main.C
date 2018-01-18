#include <iostream>
#include "string.h"
#include "string"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "math.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include <TBenchmark.h>

#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>

#include "track_finder_main.h"

#include "TVector3.h"

/* The status bits for hits in the hitlists */
// I keep these definitions from the original code from Howard
// they were much more, but it was defined in the rtpc.h code
// 

#define HISUSED 1  /* - (was 1) used */

#define HITUNAV 1 /* reasons to NOT use a hit in the chain-linker */

using namespace std;

/*THIS IS THE EXECUTABLE VERSION OF CHAIN FINDER*/

void chain_finder()
{
  //- ********************************
  // DEFINE VARIABLES FROM ROOT FILE

  
 
  TString inputf = inputfile;
  cout<<"INPUT FILE: " <<inputf<<endl;
  
  TFile *infile = new TFile(inputf);

  infile->GetObject("tvec", RTPCTree);

  
  Entries = RTPCTree->GetEntries();
  cout<<"Entries: "<<Entries<<endl;

  RTPCTree ->SetBranchAddress("event", &event);
  RTPCTree ->SetBranchAddress("nHits", &nHits);
  RTPCTree ->SetBranchAddress("pad", &Vpad);
  RTPCTree ->SetBranchAddress("tdc", &Vtdc);
  RTPCTree ->SetBranchAddress("adc", &Vadc);  
  
  
  //***********************************
  // DEFINE VARIABLES TO SAVE ROOT FILE
  
  TString outputfile  = "ChainEvents_1.root";
  rootoutfile  = new TFile(outputfile , "recreate");
  
  cout<<"OUTPUT TMP FILE: " <<outputfile <<endl;
  
   
  //***********************************

  
 
  // Create a TTree
  chaintree =  new TTree("tvec","Tree with vectors");
  chaintree -> Branch("event",&event);
  chaintree -> Branch("nHits",&nHits);
  chaintree -> Branch("pad",&pad);
      
   
  //LOOP TO READ THE VARIABLES FROM THE ROOT FILE
  
  //here we are reading the whole file, which has many superevents
  //perhaps we want to read one event at the time, for test purposes.
  
  //In general, this is a track/chain finder, so this code only FILTER
  // and its output should be individuals chains to be analyzed later
  // with the Kalman Filter or something else.
  


  // max_entries = 1;
  max_entries = Entries;
  //  Int_t ii = 1 ;

    
  for(Int_t ii = 0; ii < max_entries; ii++)  //uncomment for the whole file.
    {
      readout(hitevent, ii);//Read the event and store in the struct variables

      bubbleSort(hitevent, NoH);//reorder the data according time data

      time_label(hitevent, NoH);//assign time slot to each hit (USELESS if using Jixie data)

      
      cout<<"NoH: "<<NoH<<endl;
      
      search(NoH, ii);
	     
	        
    }

  store_data();
  
  rootoutfile ->Close();

  cout<<"RETURN"<<endl;
    
}





// WORKING ON THE PARAMETERS
void search(Int_t max_hits, Int_t event_ii )
{

  TIDVec.clear();
  TIDVec.push_back(1);


 
  for(int t = 0; t <time_slot; t += StepSize) //change it!!
    {
      concpads = 0;
      //    cout<<"ENTER LOOP "<<t<<endl;
      //loop over all pads for each time slice

      for(unsigned int p = 0; p <  vHit[t].size(); p++)
	{
	  Pad = vHit[t][p];//PadNum.get(p);
		

	  //I HAVE THE HITS ALREADY DISCRIMINATED

	  // ADC = ADCMap.get(Pad)[t];
	  // //only pads which have an ADC value above threshold will be assigned a TID
	  // if(ADC > thresh)


	  //returns the row and column of the current Pad
	  PadPhi = fRow(Pad);
	  PadZ   = fCol(Pad);
	
	  //loop through all TID's in a vector which will grow to include all future TID's
	  for(unsigned int i = 0; i < TIDVec.size(); i++)
	    {
	      breakTIDloop = false;
	      TID = TIDVec.at(i);
		
	     
				
	      //if TID is already in the map
	      if(TIDMap.count(TID))
		{
		  //  cout<<"TIDMap.count(TID) ENTER "<< TID<<endl;

		  if(t>0)
		    {
		      padindexmax = max(TIDMap[TID][t].size(), TIDMap[TID][t-StepSize].size());
		    }
		  else
		    {
		      padindexmax = TIDMap[TID][t].size();
		    }

		  //loop through all pads in TIDMap and compare there row and column to current Pad

		  
		  for(unsigned int padindex = 0; padindex < padindexmax; padindex++)
		    {

		      //Check previous time slice for adjacency
		      if(t>0)
			{
			  if( padindex < TIDMap[TID][t-StepSize].size() )
			    {
			      checkpadphiprev = fRow( TIDMap[TID][t-StepSize][padindex] );
			      checkpadzprev   = fCol( TIDMap[TID][t-StepSize][padindex] );
										
			      if(abs( checkpadphiprev - PadPhi ) >= ( Num_of_Rows - adjthresh ) )
				{
				  if(checkpadphiprev > PadPhi)
				    {
				      checkpadphiprev -= Num_of_Rows;
				    }
				  else
				    {
				      PadPhi -= Num_of_Rows;
				    }
				}
							
			      if( ( abs( checkpadphiprev - PadPhi ) <= adjthresh ) && 
				  ( abs( checkpadzprev   - PadZ )   <= adjthresh ) )
				{
				  TIDMap[TID][t].push_back(Pad);
				  breakTIDloop = true;
				  break;
				}
			    }
			} //if(t>0)

		      if(padindex < TIDMap[TID][t].size())
			{
			  checkpadphi = fRow(TIDMap[TID][t][padindex]);								
			  checkpadz   = fCol(TIDMap[TID][t][padindex]);
									
			  if(abs(checkpadphi - PadPhi) >= (Num_of_Rows - adjthresh))
			    {
			      if(checkpadphi > PadPhi)
				{
				  checkpadphi -= Num_of_Rows;
				}
			      else
				{
				  PadPhi -= Num_of_Rows;
				}
			    }
							
			  //Check current time slice for adjacency
			  if( ( abs( checkpadphi - PadPhi) <= adjthresh ) && 
			      ( abs( checkpadz   - PadZ)   <= adjthresh ) )
			    {
			      TIDMap[TID][t].push_back(Pad);
			      breakTIDloop = true;
			      break;
			    }
			}
		    }//  for(int padindex = 0						
		} // if(TIDMap.containsKey(TID))


		  //TID not already in map

	      else
	
		{
		  // TIDMap.insert(TID, new HashMap<Integer, Vector<Integer>>());
		
		  // for(int time = 0; time < TrigWindSize; time += StepSize)
		  //   {
		  //     TIDMap[TID].put(time, new Vector<>());//add TID to map
		  //   }
		  cout<<"TID: "<<TID<<" t: "<< t <<endl;
		  TIDMap[TID][t].clear();
		  TIDMap[TID][t].push_back(Pad);
		  TIDVec.push_back(TID+1);
		  break;

		}// if(TIDMap.containsKey(TID))

	      if(breakTIDloop){break;}

	    }//for(int i = 0; i < TIDVec.size(); i++)

	  concpads++;
	  if(concpads > maxconcpads) 
	    { 
	      maxconcpads = concpads;
	      maxconctime = t;					
	    }
	    
				
	}// 	for(unsigned int p
		
    }// for(int t = 0; 


  

  // I NEED TO UNDERSTAND WHY THE LOOP DO NOT BEHAVE AS I EXPECT  

  cout<<"size TIDMAP: "<<TIDMap.size()<<endl;
  cout<<"size TIDMAP[1]: "<<TIDMap[1].size()<<endl;
  cout<<"size TIDMAP[2]: "<<TIDMap[2].size()<<endl;


}



void store_data()
{

  tid_1 = TIDMap.size();
  for(int sz = 1; sz <= tid_1; sz++)
    {
      tid_2 = TIDMap[sz].size();
      for(int sz2 = 0; sz2 < tid_2; sz2++)
	{
	  tid_3 = TIDMap[sz][sz2].size();
	  for(int sz3 = 0; sz3 < tid_3; sz3++)
	    {
	      if(TIDMap[sz][sz2][sz3]>=0)
		pad.push_back(TIDMap[sz][sz2][sz3]);
	      cout<<"pad("<<sz<<","<<sz2<<","<<sz3<<"): "<< TIDMap[sz][sz2][sz3]<<endl;
	    }
	}
      event = sz-1;
      chaintree -> Fill();
    }
  
  rootoutfile -> Write();
  
}


  



void readout(HitStruct readevent[], Int_t ii)
{
  
  //This function is designed to readout two events and store
  //its contents in a expanded variable array with hit_1+hit_2 elements.
  //The second event is marked with a flag 'Ev_pos'
  //when its switch to the first position, first is checked
  //if the hit was used before. The first hit in the new iteration,
  //corresponds to the hit_1 position (from 0) in the previous iteration
  //and this is the hit we want to consider if was used or not.
  //Then it is filled as normal.
   

  // CLEAN THE VECTORS BEFORE READ THE ENTRY
    
  Vpad->clear();
  Vtdc->clear();
  Vadc->clear();

  RTPCTree->GetEntry(ii);
  
  NoH = nHits;//Number of Hits (NoH).From the event

  cout<< Vtdc->size()<<endl;
  cout<< Vadc->size()<<endl;
  cout<< Vpad->size()<<endl;


  for(Int_t p = 0; p<NoH; p++)
    {
      readevent[p].PAD  = Vpad->at(p);
      readevent[p].tdc  = Vtdc->at(p);
      readevent[p].adc  = Vadc->at(p);
      readevent[p].Status = 0;


      if(p>0)
	{
	  if( readevent[p-1].tdc != readevent[p].tdc)
	    {
	      time_scale_array. push_back(readevent[p-1].tdc);
	    }
	}
    }
    


  // These lines sort the time value and later remove the repeated values, keeping only simply values
  std::sort(time_scale_array.begin(), time_scale_array.end() );
  
  time_scale_array.erase( unique(time_scale_array.begin(), time_scale_array.end() ), time_scale_array.end() );
  
  
}




void time_label(HitStruct readevent[], int n)
{
  time_slot = 0;

  readevent[0].timeslot = time_slot;
  
  for(Int_t p = 1; p<n; p++)
    {
      if (readevent[p].tdc != readevent[p-1].tdc)
	{
	  readevent[p].timeslot = ++time_slot;
	}
      else
	{
	  readevent[p].timeslot = time_slot;
	}
  

    }
  cout<<"max Time Slot: "<<time_slot<<endl;
}



//*********************************************************************
//Algorithm from http://www.algolist.net/Algorithms/Sorting/Bubble_sort

void bubbleSort(HitStruct readevent[], int n)
{

  cout<<"ENTERING BUBBLESORT"<<endl;

  bool swapped = true;
  int j = 0;
  int tmp_tdc, tmp_PAD;
  double  tmp_adc;

  while (swapped)
    {
      swapped = false;
      j++;
      for (int i = 0; i < n - j; i++)
	{

	  //we sort only in time (tdc)
	  if(readevent[i].tdc > readevent[i+1].tdc)
	    {
	      
	      tmp_tdc = readevent[i].tdc;
	      readevent[i].tdc = readevent[i+1].tdc;
	      readevent[i+1].tdc = tmp_tdc;

	      tmp_adc = readevent[i].adc;
	      readevent[i].adc = readevent[i+1].adc;
	      readevent[i+1].adc = tmp_adc;

	      tmp_PAD = readevent[i].PAD;
	      readevent[i].PAD = readevent[i+1].PAD;
	      readevent[i+1].PAD = tmp_PAD;
	      
	      swapped = true;
	    }
	}
    }

  Int_t r = 0;

  vHit.clear();


    
  for (unsigned int t = 0; t <= time_scale_array.size(); t++)
    {
      while((readevent[r].tdc-840)/120 == t)
	{
	  //	  cout<<t<<" "<<readevent[r].PAD<<endl;
	  vPAD.push_back(readevent[r].PAD);
	  r++;
	}

      //      cout<<"vPAD size: "<<vPAD.size()<<endl;
      vHit.push_back(vPAD);
      vPAD.clear();  
    }

  Int_t coi = 0;
  
  for(unsigned int qr = 0; qr < vHit.size(); qr++)
    {

      coi = coi + vHit[qr].size();
	
    }
  cout<< "total number: " << coi << endl;


  cout<<"LEAVING BUBBLESORT"<<endl;


}

//*********************************************************************






//*********************************************************
// deconvolute the PAD information in (col,row) information
// it is supposed that the index grows in the Z-direction

int fRow(Int_t pad_info)
{
  rowT = pad_info/ncol;

  if (pad_info<0) rowT = -99;

  return rowT;
}

int fCol(Int_t pad_info)
{	
  colT = fmod(pad_info,ncol);

  if (pad_info<0) colT = -99;
  
  return colT;
}
//*********************************************************





void variable(TString FileName)
{

    
  ifstream infile;
  string line;
  infile.open(FileName);
  if (infile)
    {
      while (!infile.eof())
	{
	  if (getline(infile,line))
	    {
	      string buf;
	      string wert;
	      stringstream ss(line);
	      //		cout << ">" << line << endl;
	      ss >> buf;
		
	      if (buf == "inputfile")	 
		{			
		  ss >> inputfile;
		}

	      if (buf == "outputfile")	 
		{			
		  ss >> outputfile;
		  if (outputfile == "") outputfile = "ChainCandidates.root"; //DEFAULT VALUE
		}


		
	    }
	}
    }
}



int main(int argc, char** argv)
{
  
  TString dataname;
  dataname = "chain.ini";
  cout << "Compiled on: " __DATE__ " " __TIME__ "." << endl;

  cout<<"\n\nINI file: "<<dataname<<endl;
  variable(dataname);

  
  TStopwatch timer;
  timer.Start();
  
  if(argc < 1) return 1;

  //    tree(argv[1]);
  //  for(Int_t k=0;k<1000;k++)
  {
    chain_finder();
    //    cout<<"cycle: "<<k<<endl;
  }
    
  timer.Stop();
  
  Double_t rtime = timer.RealTime();
  Double_t ctime = timer.CpuTime();
  printf("RealTime=%f seconds, CpuTime=%f seconds\n",rtime,ctime);
  
  return 0;
}
