// EVENT DISPLAY FOR BONuS RTPC
// ORIGINAL BY HUAN YAO 
// MODIFIED BY SCOTT BARCUS, EVAN MCCLELLAN AND CARLOS AYERBE 2016-17
// FOR THE GRINCH PROTOTYPE DETECTOR


#include "TTree.h"
#include "TGLHistPainter.h"
#include "TGLBoxPainter.h"



//#define file "../files/nt_P-10T-10_200ns.root"

//#define file "../files/1000SuperEvents_PAD.root"

//#define file "../files/nt_P-10T-10_200ns_0.6mm.root"


#define test 0; //0 for simul data



int Event;

const Int_t maxHitsPerEvent = 250;
int* Hit_Events=NULL;
TCanvas* c;

int N_evts;

int N_ROW = 180;
int N_COL = 100;

// Number of PADs on RTPC
const int N_PMT= N_ROW*N_COL;

Int_t i,j;

Int_t next_evt=0;
Int_t usr_evt;

TBox**PAD = new TBox*[N_PMT];

Int_t    eventEntries;
TTree   *eventtree;
TBranch *branch;
Int_t    nHits;
Int_t    channelID;
int      tpad;

Int_t event_cnt; 
std::vector<int>* pad=0;

/*
  Color_t color_wh[45]={kRed,
  kViolet,
  kGreen-2,
  kOrange+7,
  kBlue,
  kOrange+7,
  kCyan+2,
  kSpring-5,
  kMagenta,
  kRed+3,
  kViolet-7,
  kGreen-6,
  kOrange+3,
  kBlue-2,
  kPink-9,
  kCyan+2,
  kYellow-6,
  kMagenta-9,
  kBlack,
  kOrange,
  kRed,
  kViolet,
  kGreen-2,
  kOrange+7,
  kBlue,
  kOrange+7,
  kCyan+2,
  kSpring-5,
  kMagenta,
  kRed+3,
  kViolet-7,
  kGreen-6,
  kOrange+3,
  kBlue-2,
  kPink-9,
  kCyan+2,
  kYellow-6,
  kMagenta-9,
  kBlack,
  kOrange,
  kRed,
  kViolet,
  kGreen-2,
  kOrange+7,
  kBlue};
*/

Color_t color_wh[9]={kRed,
		     kOrange+7,
		     kYellow,
		     kGreen,
		     kCyan,
		     kBlue,
		     kViolet,
		     kMagenta+3,
		     kBlack};


Color_t color_wh_2[5]={kRed+3,
		       kYellow-2,
		       kGreen-2,
		       kAzure-7,
		       kPink-6};

void event_display()
{
  //*************************************
  //DOUBLE LOOP FOR EVENT DISPLAY CREATION
 
  double x,y;//Lower corner pad

  Int_t k=0;

  double pad_z = 4;
  double pad_phi = 2.8;

  Int_t color = 0;

  //  cout<<"Creating Event Display"<<endl;

  //IN ORDER TO CHANGE THE NUMBER ORDER OF THE PADS
  //SWAP THE 'for' LOOPS LINES

  for (j=0; j<N_ROW; j++ ) //<----
    {
      for (i=0; i<N_COL; i++ ) //<----
	{
	  double shift = j%4;
	  x = i*pad_z+shift;
	  y = j*pad_phi;
	  
	  PAD[k] = new TBox (x, y,x+pad_z,y+pad_phi);
	  PAD[k]->SetFillStyle(0);

	  if (j%4==0 ) color++; 
	  
	  //	  PAD[k]->SetLineColor(color_wh[color%9]);
	  PAD[k]->SetLineColor(kCyan+4);
	  //	  cout<<"k: "<<k<<" X: "<<x<<" Y: "<<y<<endl;
	  if ( i==0 && j==0 )
	    PAD[k]->Draw();
	  else
	    PAD[k]->Draw("same");
	  k++;
	}
    }  //End building event display geometry.

  TText *x_tit = new TText(200,-20,"Z-direction");
  x_tit->SetTextAlign(22);
  x_tit->SetTextColor(kRed+2);
  x_tit->SetTextFont(43);
  x_tit->SetTextSize(20);
  x_tit->Draw("same");
  
  TText *y_tit = new TText(-15,280,"Phi-direction");
  y_tit->SetTextAlign(22);
  y_tit->SetTextColor(kRed+2);
  y_tit->SetTextFont(43);
  y_tit->SetTextSize(20);
  y_tit->SetTextAngle(90);
  y_tit->Draw("same");

  //  cout<<"End Event Display Creation"<<endl;
  //*************************************
}


void show(int option=0)
{
  /*
  // ************************************************************************************************************
  // command line input  (USED IN GRINCH)
  // ************************************************************************************************************
  //  string pre = "../decodefiles/minimal_"; //this is the path line to the datafiles and the prefix of the file
  
  string pre = "./minimal_";
  string post = ".root";
  
  pre.append(filenumber);
  string file_str = pre.append(post);
  const char *file = file_str.c_str();
  
  cout<< "File to analyze: " << file_str << endl;
  
  TFile *infile= new TFile(file);
  // ************************************************************************************************************
  // end command line input
  // ************************************************************************************************************
  */


  
  // hard code data file 

  if(test==0)
    {
      //   TFile *infile=new TFile("./ODU_files/hvector.root");
      //      TFile *infile=new TFile("../1000SuperEvents.root"); 
          TFile *infile=new TFile("ChainEvents_1.root");
       infile->GetObject("tvec", eventtree);
    }
  else
    {
      TFile *infile=new TFile("ChainEvents_1.root");
      eventtree=(TTree*)infile->Get("chaintree"); //change tree according your root file
    }

  //    eventtree=(TTree*)infile->Get("ep"); //change tree according your root file
  

  cout<<" Option: "<<option<<endl;   
  eventEntries = eventtree->GetEntries();
 
 
  //Get entries from eventtree
  
  cout<<"eventEntries: "<<eventEntries<<endl;
  
  N_evts = eventEntries-1;//Max index of events is number of entries-1  
  
  cout<<"PREVIOUS EVENT: #"<<event_cnt<<endl;


  //SELECTION

  if ( option==-2 )   //First Event
    { 
      fill_event(0);
      event_cnt = 0;
    }
  else if ( option==-3 ) //Last Event 
    { 
      fill_event(N_evts);
      event_cnt = N_evts;
    }
  else if ( option==0 ) //Ask Event
    { 
      Printf("Please Enter Event Number:");
      cin>>usr_evt;
      fill_event(usr_evt);
      event_cnt = usr_evt;
    }

 
  else if ( option == 1)
    {
      event_cnt++;
      if ( event_cnt<0 ) 
	{
	  Printf("You reach the first event");
	  fill_event(N_evts)
	    ;	  event_cnt = N_evts;
	}
      
      else if ( event_cnt>N_evts ) 
	{
	  Printf("You reach the last event");
	  fill_event(0);
	  event_cnt = 0;
	}
      else 
	{ 
	  fill_event(event_cnt);
	}
    }

  else if ( option == -1)
    {
      event_cnt--;
      if ( event_cnt<0 ) 
	{
	  Printf("You reach the first event");
	  fill_event(N_evts);
	  event_cnt = N_evts;
	}
      
      else if ( event_cnt>N_evts ) 
	{
	  Printf("You reach the last event");
	  fill_event(0);
	  event_cnt = 0;
	}
      else 
	{ 
	  fill_event(event_cnt);
	}
    }

   
}



void fill_event(Int_t Event)
{
  cout<<"Event # = "<<Event<<endl;  

  Int_t pad_test[5]={0, 179, 360, 17999, 17820};

  // Int_t pad_test[5]={157, 57, 14957, 14958, 14858};
  
  nHits = 0;
  
  if ( c )
    {
      c -> Clear();
    }
  else 
    {
      c = new TCanvas("c");
    }
  
  c->Range(-30,-40, N_COL*4 +10, N_ROW*2.8+10);

  event_display();

  if(test==0)
    {
      eventtree ->SetBranchAddress("pad", &pad); 
    }
  
  eventtree->GetEntry(Event); 
  
  if(test==0)
    {
      //      Int_t numHits = eventtree->GetLeaf("nHits")->GetValue();
      Int_t numHits = pad->size();
    }
  else
    {
      Int_t numHits = eventtree->GetLeaf("Hit")->GetValue();
    }
  
  cout << "numHits = " << numHits << endl;

  
  for(int j = 0; j < numHits; j++)
    {  
      
      
      //  
      if(test==0)
	{
	  
	  tpad  = pad->at(j);
	}
      else
	{
	  tpad = eventtree->GetLeaf("StepID") ->GetValue(j); //CHANGE NAME OF VARIABLE ACCORDING YOUR ROOT DATA
	}    
      
            
      if(tpad > 0 && tpad < 18000)//THIS IS NECESSARY SINCE THERE ARE PADs NUMBERS OUT OF RANGE
	{
	  Int_t d = tpad;	
	  
	  PAD[d]->SetFillStyle(1001);	
	  //	  PAD[d]->SetFillColor(color_wh_2[j%5]);
	  PAD[d]->SetFillColor(kRed);
	  PAD[d]->Draw("same");
	}
    }
  
  c->SetTitle(Form("Event %d",Event));


}

//Global variable (datafile number ONLY IF COMMAND INPUT IS USED)
string filenumber;

void BONuS_Event_Display() 
{

  // command line input.
  // cout<<"File number to analyze: ";
  // cin>>filenumber;
  
  
  //Create GUI
  bar = new TControlBar("vertical", "Event Displayer");
  bar->AddButton("Event Displayer", "");
  bar->AddButton("", "");
  bar->AddButton("First", "show(-2)");
  bar->AddButton("Next", "show(1)");
  bar->AddButton("Prev", "show(-1)");
  bar->AddButton("Last", "show(-3)");
  bar->AddButton("", "");
  bar->AddButton("Event #", "show(0)");
  bar->AddButton("", "");
  bar->AddButton("Exit",".q");
  bar->Show();
 
}
