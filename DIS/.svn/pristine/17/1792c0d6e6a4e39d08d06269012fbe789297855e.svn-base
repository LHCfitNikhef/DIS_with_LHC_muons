// main01.cc is a part of the PYTHIA event generator.
// Copyright (C) 2012 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program. It fits on one slide in a talk. 
// It studies the charged multiplicity distribution at the LHC.

#include "Pythia8/Pythia.h"
#include "DISPowhegHooks.h"
#include "DISPowhegHooksVincia.h"
// #include "Pythia8Plugins/PowhegHooks.h"
// #include "Pythia8Plugins/PowhegHooksVincia.h"
#include "Pythia8Plugins/LHAFortran.h"
#include <sstream>

#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "Pythia8Plugins/HepMC2.h"
#include <cassert>

const int dglap=0, dipole=1, vincia=2, vinciaGlobal=3;

using namespace Pythia8;


extern "C" {
  extern struct {
    int tune;
    int had;
    int hepmc;
    int MPI;
    int shower;   // Pythia shower used: 0=default shower, 1=dipole recoil, 2=vincia
    int QED;
    int LO;

  } cpy8pars_;
}

class myLHAupFortran : public LHAupFortran {
protected:
  // User-written routine that does the intialization and fills heprup.
  virtual bool fillHepRup() {return true;}

  // User-written routine that does the event generation and fills hepeup.
  virtual bool fillHepEup() {return true;}

};


class MyUserHooks : public UserHooks {

public:

  MyUserHooks() {cout << "Setting up Hook";}

  // Destructor deletes anti-kT jet finder.
  ~MyUserHooks() {;}

  // Allow process cross section to be modified..
private:

};

Pythia pythia;

// Add in user hooks for shower vetoing
shared_ptr<UserHooks> powhegHooks;
bool loadHooks = false;

shared_ptr<HepMC::IO_GenEvent> ptrHepMC = NULL;
shared_ptr<HepMC::Pythia8ToHepMC> toHepMC  = NULL;


std::shared_ptr<myLHAupFortran> LHAinstance{new myLHAupFortran()};


extern "C" {
  // F77 interface to pythia8
  void pythia_option0_(char *string) {
    pythia.readString(string);
  }

  void pythia_init_() {

    // Vincia
    if(cpy8pars_.shower >= vincia){
      std::cout<<"Using Vincia shower";
      //use a massive charm
      pythia.readString("Vincia:nFlavZeroMass = 3");

      if(cpy8pars_.shower == vinciaGlobal){
	// transverse momentum imbalance due to ISR shared globally
	pythia.readString("Vincia:kineMapIF = 2");  //Gluon emissions use a probabilistic selection between the global and local maps.
	                                            //Antennae that only contain initial-state singularities always use the global one.
	                                            //Antennae that only contain final-state singularities always use the local one.
	std::cout<<" with global recoil for IF \n";
      }else{
	std::cout<<"\n";
      }
      pythia.readString("PartonShowers:model = 2"); // Use Vincia’s shower model.
      if(cpy8pars_.QED ==0) {
	std::cout<<"Switching OFF QED shower in Vincia \n"; 
	pythia.readString("Vincia:ewMode = 0");  // Switch off QED/EW showers.
      }else if(cpy8pars_.QED == 2) {
	std::cout<<"Using the most ``coherent'' QED shower in Vincia \n"; 
	pythia.readString("Vincia:ewMode = 2");  // Use most coherent QED
      }else{
	std::cout<<"Using the fastest QED shower in Vincia \n"; 
	pythia.readString("Vincia:ewMode = 1");  // Default
      }
  
    }
    //Dipole shower
    else{
      std::cout<<"Using standard Dipole shower";
      if(cpy8pars_.shower == dipole){
	std::cout<<" with fully local recoil \n";
	pythia.readString("SpaceShower:dipoleRecoil = on");
      }else{
	std::cout<<" with global recoil for IF\n";
      }	
      if(cpy8pars_.QED ==0) {
	std::cout<<"Switching OFF QED shower in Standard dipole shower \n"; 
	pythia.readString("SpaceShower:QEDshowerByQ = off"); // From quarks.        
	pythia.readString("SpaceShower:QEDshowerByL = off"); // From Leptons.      
	pythia.readString("TimeShower:QEDshowerByQ = off"); // From quarks.         
	pythia.readString("TimeShower:QEDshowerByL = off"); // From Leptons.
	pythia.readString("TimeShower:QEDshowerByGamma = off"); // Prevent gamma -> f fbar
	pythia.readString("TimeShower:QEDshowerByOther  = off");
      }else{
	std::cout<<"Using QED shower in Standard dipole shower \n"; 
      }
    }

    // tune
    if(cpy8pars_.tune > 0) {
      std::stringstream ss;
      ss << cpy8pars_.tune;
      pythia.readString("Tune:pp = "+ss.str());
      cout << "pythia8F377: setting pythia tune "<<cpy8pars_.tune<<"\n";
    }

    // Change C hadronization parameters
    /*pythia.readString("StringZ:rFactC  = 1.0");
    pythia.readString("StringZ:useNonstandardC = on");
    pythia.readString("StringZ:aNonstandardC = 1.827");
    pythia.readString("StringZ:bNonstandardC = 0.837");*/

    

   if(cpy8pars_.had == 0) {
      pythia.readString("HadronLevel:All = off");
      std::cout<<"Switching off hadronization \n";
   }else{
     if(cpy8pars_.had == 1){
       std::cout<<"Including hadronization \n";
     }else{
       std::cout<<"Including hadronization, but not the unstable hadrons decay \n";
       pythia.readString("HadronLevel:Decay = off");  //switch off unstable hadron decay
     }
    }
    // uncomment to set pi0 stable (this avoid a plethora
    // of photons in the final state...)
    //pythia.readString("111:mayDecay = off");
    // same for the higgs
    //pythia.readString("25:mayDecay = off");
    // //D,Ds hadrons
    // pythia.readString("411:mayDecay = off");
    // pythia.readString("421:mayDecay = off"); 
    // pythia.readString("413:mayDecay = off"); 
    // pythia.readString("423:mayDecay = off"); 
    // pythia.readString("431:mayDecay = off"); 
    // pythia.readString("433:mayDecay = off"); 
    // pythia.readString("-411:mayDecay = off");
    // pythia.readString("-421:mayDecay = off"); 
    // pythia.readString("-413:mayDecay = off"); 
    // pythia.readString("-423:mayDecay = off"); 
    // pythia.readString("-431:mayDecay = off"); 
    // pythia.readString("-433:mayDecay = off"); 
    
    if(cpy8pars_.MPI == 1){
      std::cout<<"Include the underlying event \n";
    }else{
      std::cout<<"Switch off the underlying event \n";
      pythia.readString("PartonLevel:MPI = off");
    }



    // Load configuration file
    if(cpy8pars_.LO ==1){

      cout<<"LO. powheg hooks not loaded. Use pTmaxMatch = 2 [max starting scale]"<<endl;
      if (cpy8pars_.shower >= vincia){
	pythia.readString("Vincia:pTmaxMatch = 2");
      }else{
	pythia.readString("SpaceShower:pTmaxMatch = 2");
	pythia.readString("TimeShower:pTmaxMatch = 2");
      }
      pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
      
    }else{

          //POWHEG TAILORED SETTINGS
      // -----> from main31.cc
      // Add further settings that can be set in the configuration file
      pythia.settings.addMode("POWHEG:nFinal",    2, true, false, 1, 10);  // 2 born level partons (undoing W decay)
      pythia.settings.addMode("POWHEG:veto",      0, true, true,  0, 2);
      pythia.settings.addMode("POWHEG:vetoCount", 10, true, false, 0, 0);  //check some emissions
      pythia.settings.addMode("POWHEG:pThard",    0, true, true,  0, 2);
      pythia.settings.addMode("POWHEG:pTemt",     0, true, true,  0, 2);
      pythia.settings.addMode("POWHEG:emitted",   0, true, true,  0, 3);
      pythia.settings.addMode("POWHEG:pTdef",     0, true, true,  0, 2);
      pythia.settings.addMode("POWHEG:MPIveto",   0, true, true,  0, 1);
      pythia.settings.addMode("POWHEG:QEDveto",   0, true, true,  0, 2);   // 0= do not veto QED radiation
      pythia.settings.addMode("POWHEG:dis_map",   1, true, true,  0, 2);
      
      pythia.readFile("main31.cmnd");



      int dis_map = 0;
      if (powheginput("#flg_dis") == 0){
	dis_map = 0;
      }else{
	if(powheginput("#iupperisr") == 2.0){
	  dis_map = 2;
	}else{
	  dis_map = 1;
	}
      }
      std::stringstream ds;
      ds << dis_map;
      pythia.readString("POWHEG:dis_map = "+ds.str());
      

    // Read in main settings
      int nEvent      = pythia.settings.mode("Main:numberOfEvents");
      int nError      = pythia.settings.mode("Main:timesAllowErrors");
      // Read in POWHEG settings
      int vetoMode    = pythia.settings.mode("POWHEG:veto");
      int MPIvetoMode = pythia.settings.mode("POWHEG:MPIveto");
      loadHooks  = (vetoMode > 0 || MPIvetoMode > 0);
      
    
      if (loadHooks) {
	cout<<"powheg hooks loaded. Use pTmaxMatch = 2 [scalap is a veto scale]"<<endl;
	
	// Set ISR and FSR to start at the kinematical limit
	if (vetoMode > 0) {
	  if(cpy8pars_.shower >= vincia){
	    pythia.readString("Vincia:pTmaxMatch = 2");
	  }else{
	    pythia.readString("SpaceShower:pTmaxMatch = 2");
	    pythia.readString("TimeShower:pTmaxMatch = 2");
	  }
	}
	// Set MPI to start at the kinematical limit
	if (MPIvetoMode > 0) {
	  pythia.readString("MultipartonInteractions:pTmaxMatch = 2");
	}
	
	if (cpy8pars_.shower >= vincia){
	  powhegHooks = make_shared<PowhegHooksVincia>();
	}else{
	  powhegHooks = make_shared<PowhegHooks>();
	}
	
	pythia.setUserHooksPtr((UserHooksPtr)powhegHooks);
      } else {
      // do this (wimpy shower) if hooks are not loaded.
      // (which can be obtained setting POWHEG:veto=0)
	cout<<"powheg hooks not loaded. Use pTmaxMatch = 1 [scalup = starting scale]"<<endl;
	if (cpy8pars_.shower >= vincia){
	  pythia.readString("Vincia:pTmaxMatch = 1");
	}else{
	  pythia.readString("SpaceShower:pTmaxMatch = 1");
	  pythia.readString("TimeShower:pTmaxMatch = 1");
	}
      }
    }

    if(cpy8pars_.hepmc > 0) {
      toHepMC  = make_shared<HepMC::Pythia8ToHepMC>();
      ptrHepMC = make_shared<HepMC::IO_GenEvent>();
    }


    pythia.readString("Beams:frameType = 5");
    pythia.readString("Check:beams False");
    // Set up incoming beams, for frame with unequal beam energies.
    //pythia.readString("Beams:frameType = 2");
    //// BeamA = proton.
    //pythia.readString("Beams:idB = 2212");
    //pythia.settings.parm("Beams:eB", 820.0);
    //// BeamB = electron.
    //pythia.readString("Beams:idA = 11");
    //pythia.settings.parm("Beams:eA", 27.5);
    ////    pythia.setUserHooksPtr(MyHook);


    pythia.setLHAupPtr(LHAinstance); 
    LHAinstance->setInit();  

 

    pythia.init();

    
  }

  void pythia_next_(int & iret){
  // Begin event loop. Generate event. Skip if error. List first one.
    iret = pythia.next();
    if(iret == 0) {
      // bad event
      if (pythia.info.atEndOfFile()) {
	cout << "Info: end of Les Houches file reached" << endl;
	iret=-1;
	return ;
      }
    }else {

      if(toHepMC) {
	// Construct new empty HepMC event and fill it.
	// Units will be as chosen for HepMC build, but can be changed
	// by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
	HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
	toHepMC->fill_next_event( pythia , hepmcevt );
	// // FIX BELOW to allow multiple weights  ---> NOT TESTED
	// HepMC::WeightContainer& weights = hepmcevt->weights();
	// weights.clear();
	// weights["000"] = _lastWeight.central();
	// Write out the HepMC event.
	*ptrHepMC << hepmcevt;
	delete hepmcevt;
      }
    
  }

    // RIVET
    //rivet(); 
    
    // // uncomment the following is want to use pipe to analyze events with Rivet.
    // // Construct new empty HepMC event and fill it.
    // // Units will be as chosen for HepMC build, but can be changed by arguments, e.g. GenEvt( HepMC::Units::GEV, HepMC::Units::MM)
    // HepMC::GenEvent* hepmcevt = new HepMC::GenEvent();
    // ToHepMC.fill_next_event( pythia , hepmcevt );
    // // Write the HepMC event to file. Done with it.
    // ascii_io << hepmcevt;
    // delete hepmcevt;             
  }

  void pythia_to_hepevt_(const int &nmxhep, int & nhep, int * isthep,
			 int * idhep,
			 int  (*jmohep)[2], int (*jdahep)[2],
			 double (*phep)[5], double (*vhep)[4]) {
    nhep = pythia.event.size();
    if(nhep>nmxhep) {cout << "too many particles!" ; exit(-1); }
    for (int i = 0; i < pythia.event.size(); ++i) {
      *(isthep+i) = pythia.event[i].statusHepMC(); //!ER: !!!!!!
      //*(isthep+i) = pythia.event[i].status();
      
      *(idhep+i) = pythia.event[i].id();
      // All pointers should be increased by 1, since here we follow
      // the c/c++ convention of indeces from 0 to n-1, but i fortran
      // they are from 1 to n.
      (*(jmohep+i))[0] = pythia.event[i].mother1() + 1;
      (*(jmohep+i))[1] = pythia.event[i].mother2() + 1;
      (*(jdahep+i))[0] = pythia.event[i].daughter1() + 1;
      (*(jdahep+i))[1] = pythia.event[i].daughter2() + 1;
      (*(phep+i))[0] = pythia.event[i].px();
      (*(phep+i))[1] = pythia.event[i].py();
      (*(phep+i))[2] = pythia.event[i].pz();
      (*(phep+i))[3] = pythia.event[i].e();
      (*(phep+i))[4] = pythia.event[i].m();
    }
    // override mother of very first event, set to 0
    *(jmohep)[0] = 0 ;
    *(jmohep)[1] = 0 ;

    // print every event 
    //pythia.event.list(); 
  }

  void pythia_stat_() {
    pythia.stat();

  }

}

