#define MakeTrees_cxx
#include "MakeTrees.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TGaxis.h>
#include <vector>
#include <iostream>
#include <algorithm>

void MakeTrees::GetTriggerSlopes(){
   const int slopeStartTick = 1900;
   const int slopeEndTick   = 2900;
   const int maxWidth     = slopeEndTick - slopeStartTick;
   const int minWidth     = 125;
   const double minSlope = 0.005;
   const double minAmp = -0.4;
   
   TFile *file = new TFile("SlopeHistogramsSmall.root", "RECREATE");
   
   const double pedestal = -0.005;
   //const double pedestal = -0.008; //use for 10 mv / Div
   const double triggerThreshhold = -0.01;
   double amplitude = 0.0;
   int startTick   = 0; 
   int endTick     = 0;
   
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   TH1F *waveSlope       = new TH1F("waveSlope", "", N_SAMPLES, 0, N_SAMPLES);
   TH1F *triggerSlopes   = new TH1F("triggerSlopes", "", 200, 0, 200);
   TH1F *triggerSlopesOneUs   = new TH1F("triggerSlopesOneUs", "", 200, 0, 200);
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      //copy waveforms into vectors (make my life easier)
      
      std::vector<double> waveForm;
      //Restricted waveform of the trigger pulse      
      waveForm.assign(ch1wfms, ch1wfms+N_SAMPLES);
      startTick = FindPulseStartTime(waveForm, 0, triggerThreshhold);      
      endTick   = startTick + maxWidth;
      double amplitude = getAmplitude(waveForm, startTick, endTick);
      
      
      if( amplitude <= minAmp){
         int quietTicks    = 0;
	 int maxQuietTicks = 0;
	 for(int tick = startTick; tick < endTick + 1; ++tick){
           double slope = getSlope(waveForm, tick);
	   
	   if(fabs(slope) < minSlope){	   	   
	     ++quietTicks;
	     if(quietTicks > maxQuietTicks)
	       maxQuietTicks = quietTicks;	       	     
           }
	   
	   else
	     quietTicks = 0;	   	   
         }
	 triggerSlopes->Fill(maxQuietTicks);
	 
	 maxQuietTicks = 0;
	 quietTicks = 0;
	 
	 for(int tick = startTick; tick < startTick + minWidth + 1; ++tick){
           double slope = getSlope(waveForm, tick);
	   
	   if(fabs(slope) < minSlope){	   	   
	     ++quietTicks;
	     if(quietTicks > maxQuietTicks)
	       maxQuietTicks = quietTicks;	       	     
           }
	   
	   else
	     quietTicks = 0;	   	   
         }
	 triggerSlopesOneUs->Fill(maxQuietTicks);
	 
      //std::cout << "Most Consecutive quiet ticks = " << maxQuietTicks << std::endl;
      
      
     }//end of if	 
   }
   file->Write();   

}//end of GetTriggerSlopes

void MakeTrees::Loop()
{

   const double pedestal = -0.005;
   const double triggerThreshhold = -0.01;
   
   TFile *ntupleFile = new TFile("SlowFill50mvDivNewPulse.root", "RECREATE");
   
   TTree *tree0 = new TTree("wavedata","Wavedata");
   
   ntupleFile->cd();

   int currentFile = -1;
   int waveNumber  = 0;
   long runDate     = 0;
   int startTick   = 0; 
   int endTick     = 0;
   double area      = 0.0;
   double amplitude = 0.0;
   int pulseNumber = 0;
   
   tree0->Branch("WaveNo", &waveNumber, "WaveNo/I");
   tree0->Branch("PulseNo", &pulseNumber, "PulseNo/I");
   tree0->Branch("StartTick", &startTick, "StartTick/I");
   tree0->Branch("EndTick", &endTick, "EndTick/I");
   tree0->Branch("Date", &runDate, "Date/L");
   tree0->Branch("Amplitude", &amplitude, "Amplitude/D");
   tree0->Branch("Area", &area, "Area/D");
   //tree0->Branch("TriggerSlope", TriggerSlope, Form("TriggerSlope[%d]/D", totalWidth);      
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
      waveForm.assign(ch1wfms, ch1wfms+N_SAMPLES);
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
	
	
	std::cout << runDate << std::endl;
      }	
      
      startTick = FindPulseStartTime(waveForm, 0, triggerThreshhold);
      endTick   = FindPulseEndTime(waveForm, startTick, pedestal);
      area    = getArea(waveForm, startTick, endTick);
      amplitude = getAmplitude(waveForm, startTick, endTick);
      if(endTick == -1){
         tree0->Fill();
      } 
      

	  
      while(endTick < N_SAMPLES && endTick != -1 && startTick != -1){
      
          /*std::cout << "Start tick: " << startTick << std::endl;
          std::cout << "End tick: "   << endTick << std::endl;
          std::cout << "Width: "      << (endTick - startTick) << std::endl;
          std::cout << "Area: "       << area << std::endl;
          std::cout << "Amplitude: "       << amplitude << std::endl;
          std::cout << "Area / Amp: "       << area/amplitude << std::endl;
          std::cout << "Wavenumber: " << waveNumber << std::endl;*/

	  
	  startTick = FindPulseStartTime(waveForm, endTick+1, triggerThreshhold);
          endTick   = FindPulseEndTime(waveForm, startTick, pedestal);
          area      = getArea(waveForm, startTick, endTick);
          amplitude = getAmplitude(waveForm, startTick, endTick);
	  
	  ++pulseNumber;
	  
	  tree0->Fill();
	  	
	  	              
      }//end of while
     
   }//end of big loop
   tree0->Write();
   ntupleFile->Close();
   
   
}

void MakeTrees::ViewWave(){
   
   TH1F *waveForm    = new TH1F("waveFrom", "", N_SAMPLES, 0, N_SAMPLES);
   TH1F *waveSlope   = new TH1F("waveSlope", "", N_SAMPLES, 0, N_SAMPLES);
   TH1F *h_startTick = new TH1F("startTick", "", N_SAMPLES, 0, N_SAMPLES);
   TH1F *h_endTick   = new TH1F("endTick", "", N_SAMPLES, 0, N_SAMPLES);
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

   Long64_t nbytes = 0, nb = 0;
   
   
   TCanvas *can = new TCanvas("can", "", 1500, 800);
   //for (Long64_t jentry=0; jentry<nentries;jentry++) {
   Long64_t jentry=0;
   int waveNo = 0;
   
   
   while(jentry < nentries && jentry >= 0){
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      
      std::string cmd;
      ++waveNo;
      std::vector<double> waveData;
      waveData.assign(ch1wfms, ch1wfms+N_SAMPLES);
      
      for(int tick = 1; tick < N_SAMPLES; ++tick){
         waveForm->Fill(tick, ch1wfms[tick]);
	 waveSlope->Fill(tick, getSlope(waveData, tick) );
	 //cout << getSlope(waveData, tick) << endl;
      
      }
      
      
      
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
      can->Divide(1,2);
      can->cd(1);
      waveForm->GetXaxis()->SetRangeUser(1900, 2900);
      waveForm->SetMaximum(0.01);
      //waveForm->SetMinimum(-0.01);
      waveSlope->SetMarkerColor(kGreen+2);
      waveSlope->SetMarkerStyle(kFullDotMedium);
      
      //waveForm->SetMinimum(-0.1);
      waveForm->Draw("hist");
      //legend->Draw("same");
      
      h_endTick->Draw("hist same");
      h_startTick->Draw("hist same");
                
      waveForm->SetTitle(Form("Wave Number: %d", waveNo) );
      //gPad->SetLogy(1);
      can->cd(2);
      //waveSlope->GetXaxis()->SetRangeUser(1900, 2900);
      waveSlope->GetYaxis()->SetTitle("Slope");
      waveSlope->SetMinimum(-0.05);
      waveSlope->SetMaximum(0.05);
      waveSlope->Draw("P same");
      
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
      waveSlope->Reset();
      
      //delete can;
                     
   }

  
}//end of ViewWave


double MakeTrees::getAmplitude(std::vector<double> wave, int startTick, int endTick){
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

double MakeTrees::getArea(std::vector<double> wave, int startTick, int endTick){
    if(startTick == -1 || endTick == -1)
      return 0.0;
    
    double sum = 0.0;
    for(int point = startTick; point < endTick; ++point){
       sum += fabs(wave[point]);
    }
    
    //std::cout << "Max tick: " << tick << std::endl; 
    return sum;

} // end of getAmplitude

double MakeTrees::getSlope(std::vector<double> wave, int tick){
    if(tick != N_SAMPLES-1)
      return (wave[tick + 1] - wave[tick]);
    else
      return 0.0;  

} // end of getAmplitude

//Try to find when the pulse began. Define this as when the pulse falls above the trigger level
int MakeTrees::FindPulseStartTime(std::vector<double> ADCs, int MPulseStartTime, double Threshold){


   
   for(int i = MPulseStartTime; i < N_SAMPLES; ++i){
     if( ADCs[i] <= Threshold)
       return i;
      
   }//end of for loop
   
   return -1; //indicates a failure

} //End of FindEndTime

//Try to find when the pulse ended. Define this as when the pulse drops to pedestal levels
int MakeTrees::FindPulseEndTime(std::vector<double> ADCs, int MPulseStartTime, double Pedestal){

if(MPulseStartTime == -1)
  return -1;

int returnTick = 0;

const int max_ticks = 18;
const int max_width = 999999999; //do not set a maximum

int quietTicks = 0;
int currentwidth = 0;



for(int i = MPulseStartTime; i < N_SAMPLES; ++i){
   ++currentwidth;
   if( ADCs[i] >= Pedestal && ADCs[i] < 0.002){
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
