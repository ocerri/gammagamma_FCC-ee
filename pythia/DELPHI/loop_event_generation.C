#include "pythia6_gammagamma_hadrons.h"
#include "pythia6_gammagamma_hadrons.C"


void loop_event_generation(){
  for(int i=0;i<10;i++){
    pythia6_gammagamma_hadrons(10000,195.5,30);
    //gSystem->Sleep(1000);
    pythia6_gammagamma_hadrons(10000,195.5,10);
    cout << endl<< endl<< endl<< endl<<"Run " << i << " finished."<< endl<< endl<< endl<< endl;
  }
}
    
    
