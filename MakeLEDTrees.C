#define MakeLEDTrees_cxx
#include "MakeLEDTrees.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <algorithm>

void MakeLEDTrees::Loop()
{

   const double pedestal = -0.005;
   const double triggerThreshhold = -0.01;
   
   TFile *ntupleFile = new TFile("LEDTriggerPed.root", "RECREATE");
   
   TTree *tree0 = new TTree("wavedata","Wavedata");
   
   ntupleFile->cd();


   
   //const int pulseStartTick = 980;
   
   
   
   //declare all branch addresses
   int currentFile = -1;
   int waveNumber  = 0;
   long runDate     = 0;
   int startTick   = 0; 
   int endTick     = 0;
   double area      = 0.0;
   double amplitude = 0.0;
   double Pedestal[N_SAMPLES];
   int pulseNumber = 0;
   
   for(int tick = 1; tick < N_SAMPLES; ++tick){
	    Pedestal[tick] = 0.0;
   }
   
   tree0->Branch("WaveNo", &waveNumber, "WaveNo/I");
   tree0->Branch("PulseNo", &pulseNumber, "PulseNo/I");
   tree0->Branch("StartTick", &startTick, "StartTick/I");
   tree0->Branch("EndTick", &endTick, "EndTick/I");
   tree0->Branch("Date", &runDate, "Date/L");
   tree0->Branch("Amplitude", &amplitude, "Amplitude/D");
   tree0->Branch("Area", &area, "Area/D");
   tree0->Branch("Pedestal", Pedestal, Form("Pedestal[%d]/D", N_SAMPLES) );      
   //tree0->Branch("Waveform", Waveform, Form("Waveform[%d]/D", N_SAMPLES));
   
   
   std::string currentFileName;
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;

   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //copy waveforms into vectors (make my life easier)
      ++waveNumber;
      std::vector<double> waveForm;
      waveForm.assign(ch2wfms, ch2wfms+N_SAMPLES);
      //getAmplitude(waveForm);
      //std::cout << "Waveform amplitude: " << getAmplitude(waveForm) << "V" << std::endl; 
      pulseNumber = 0;
      
      if(fCurrent != currentFile){
        currentFile = fCurrent;
	currentFileName = getFileName();
	std::size_t pos = currentFileName.find("2016");
	std::string timeStamp = currentFileName.substr(pos, 16);
	//strip the timeStamp down to time
	std::string::iterator end_pos = std::remove(timeStamp.begin(), timeStamp.end(), '-');
	timeStamp.erase(end_pos, timeStamp.end());
	
	end_pos = std::remove(timeStamp.begin(), timeStamp.end(), '_');
	timeStamp.erase(end_pos, timeStamp.end());
	
	//finally, convert the string time to int time
	runDate = std::stol(timeStamp);
	
	
	//std::cout << runDate << std::endl;
      }	
      
      startTick = FindPulseStartTime(waveForm, 0, triggerThreshhold);
      
      endTick   = FindPulseEndTime(waveForm, startTick, pedestal);
      
      area    = getArea(waveForm, startTick, endTick);
      amplitude = getAmplitude(waveForm, startTick, endTick);
      //std::cout << startTick << " " << endTick << " " << amplitude << " " << area << std::endl;
      if(endTick == -1){
         for(int tick = 1; tick < N_SAMPLES; ++tick){
            Pedestal[tick] = ch2wfms[tick];	    
         }
	 
	 tree0->Fill();
      } 
      
      else{
         for(int tick = 1; tick < N_SAMPLES; ++tick){
            Pedestal[tick] = -10.0;
      
         }
      }

	  
      while(endTick < N_SAMPLES && endTick != -1 && startTick != -1){

	  tree0->Fill();
	  startTick = FindPulseStartTime(waveForm, endTick+1, triggerThreshhold);
          endTick   = FindPulseEndTime(waveForm, startTick, pedestal);
          area      = getArea(waveForm, startTick, endTick);
          amplitude = getAmplitude(waveForm, startTick, endTick);
	  
	  ++pulseNumber;
	  	
	  	              
      }//end of while
     
   }//end of big loop
   tree0->Write();
   ntupleFile->Close();
   
   
}

void MakeLEDTrees::ViewWave(){
   
   TH1F *waveForm = new TH1F("waveFrom", "", N_SAMPLES, 0, N_SAMPLES);
   TH1F *h_startTick = new TH1F("startTick", "", N_SAMPLES, 0, N_SAMPLES);
   TH1F *h_endTick = new TH1F("endTick", "", N_SAMPLES, 0, N_SAMPLES);
   gStyle->SetOptStat(0);
   waveForm->SetLineColor(kBlack);
   h_startTick->SetLineColor(kBlue);
   h_endTick->SetLineColor(kRed);
   
   h_startTick->SetLineWidth(2);
   h_endTick->SetLineWidth(2);
   
   TLegend *legend = new TLegend(0.7, 0.85, 0.9, 0.75);
   legend->AddEntry(h_startTick, "Pulse Start Time", "l");
   legend->AddEntry(h_endTick, "Pulse End Time", "l"); 
   
   const double pedestal = -0.005;
   const double triggerThreshhold = -0.01;
   
   int startTick   = 0; 
   int endTick     = 0;
   double area      = 0.0;
   double amplitude = 0.0;
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   std::cout << nentries << std::endl;

   Long64_t nbytes = 0, nb = 0;
   
   
   TCanvas *can = new TCanvas("can", "", 1500, 800);
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   Long64_t jentry=58917;
   int waveNo = 0;
   while(jentry < nentries && jentry >= 0){
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      std::string cmd;
      ++waveNo;
      for(int tick = 1; tick < N_SAMPLES; ++tick){
         waveForm->Fill(tick, ch2wfms[tick]);
      
      }
      
      std::vector<double> waveData;
      waveData.assign(ch2wfms, ch2wfms+N_SAMPLES);
      
      startTick = FindPulseStartTime(waveData, 0, triggerThreshhold);
      endTick   = FindPulseEndTime(waveData, startTick, pedestal);
      area    = getArea(waveData, startTick, endTick);
      amplitude = getAmplitude(waveData, startTick, endTick);
      
      h_startTick->Fill(startTick-1, -100);
      h_endTick->Fill(endTick, -100);
      
      while(endTick < N_SAMPLES && endTick != -1 && startTick != -1){
	  startTick = FindPulseStartTime(waveData, endTick+1, triggerThreshhold);
          endTick   = FindPulseEndTime(waveData, startTick, pedestal);
          area      = getArea(waveData, startTick, endTick);
          amplitude = getAmplitude(waveData, startTick, endTick);
	  
	  h_startTick->Fill(startTick-1, -100);
          h_endTick->Fill(endTick+1, -100);
      
	  	
	  	              
      }//end of while
      
      //waveForm->GetXaxis()->SetRangeUser(1900, 5000);
      waveForm->Draw("hist");
      h_endTick->Draw("hist same");
      h_startTick->Draw("hist same");
      legend->Draw("same");
      waveForm->SetTitle(Form("Wave Number: %d", waveNo) );
      can->Update();
      std::cout << "Move to next waveform: [n]" <<  std::endl;
      std::cout << "Move to previous waveform: [p]" << std::endl; 
      std::cout << "Quit like a coward: [q]" << std::endl;
      std::cin >> cmd;
      
      if(cmd == "n"){
        std::cout << "yay!" << std::endl;
	++jentry;
      }
      
      else if(cmd == "p"){
        std::cout << "boo!" << std::endl;
	--jentry;

      }
      
      else if(cmd == "q"){
        std::cout << "buh bye!" << std::endl;
	break;
      }
      
      can->Clear();
      can->Update();
      waveForm->Reset();
      h_startTick->Reset();
      h_endTick->Reset();
      
      //delete can;
                     
   }

  
}//end of ViewWave


double MakeLEDTrees::getAmplitude(std::vector<double> wave, int startTick, int endTick){
    if(startTick == -1 || endTick == -1)
      return 0.0;
    double max = 0.0;
    for(int point = startTick; point < endTick; ++point){
       if(wave[point] < max){
         max  = wave[point];
       } 
    }
    //std::cout << "Max tick: " << tick << std::endl; 
    return max;

} // end of getAmplitude

double MakeLEDTrees::getArea(std::vector<double> wave, int startTick, int endTick){
    if(startTick == -1 || endTick == -1)
      return 0.0;
    
    double sum = 0.0;
    for(int point = startTick; point < endTick; ++point){
       sum += fabs(wave[point]);
    }
    
    //std::cout << "Max tick: " << tick << std::endl; 
    return sum;

} // end of getAmplitude

//Try to find when the pulse began. Define this as when the pulse falls above the trigger level
int MakeLEDTrees::FindPulseStartTime(std::vector<double> ADCs, int MPulseStartTime, double Threshold){


   
   for(int i = MPulseStartTime; i < N_SAMPLES; ++i){
     if( ADCs[i] <= Threshold)
       return i;
      
   }//end of for loop
   
   return -1; //indicates a failure

} //End of FindEndTime

//Try to find when the pulse ended. Define this as when the pulse drops to pedestal levels
int MakeLEDTrees::FindPulseEndTime(std::vector<double> ADCs, int MPulseStartTime, double Pedestal){

if(MPulseStartTime == -1)
  return -1;

int returnTick = 0;

const int max_ticks = 12;
const int max_width = 5000; //about 40 microseconds 

int quietTicks = 0;
int currentwidth = 0;



for(int i = MPulseStartTime; i < N_SAMPLES; ++i){
   ++currentwidth;
   if( ADCs[i] >= Pedestal){
        if(quietTicks == 0) returnTick = i;
	++quietTicks;	
      }
      	
   else
      quietTicks = 0;	
   
      
   if(quietTicks >= max_ticks)
      return returnTick;
      
   else if(currentwidth > max_width)
      return i;
		  
   }//end of for loop
   return -1; //indicates a failure

} //End of FindEndTime
