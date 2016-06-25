#define RadonFinder_cxx
#include "RadonFinder.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <math.h>
#include <TF1.h>
#include <TPaletteAxis.h>
#include <TAxis.h>
#include <TLatex.h>
#include <vector>
#include <TLegend.h>
#include <TColor.h>
#include <TProfile.h>
#include <TMatrixD.h>

/*Random Fit code
TF1 f1("f1","[0]*([1]*exp(-[1])/([2]*sqrt(2*TMath::Pi()))*exp(-pow(x-[3],2)/(2*pow([2],2)))+pow([1],2)*exp(-[1])/2./([2]*sqrt(4*TMath::Pi()))*exp(-pow(x-2*[3],2)/(4*pow([2],2))))\
",0.6,2.0);
  f1.SetParameters(4000.,0.2,0.3,0.8);
*/  


void RadonFinder::MakeHists(std::string fillRate)
{
 
   const double tick_to_ns = 8.0;
   const double tick_to_us = 8.0*1e-3;
   const double tick_to_s = 8.0*1e-9;
   const int maxDepth = 20;
   const double minAmp = 0.01;
   const double maxAmp = 0.1;
   const int maxWidth = 5;
   const int maxDays = 15;
   const double startDay = 0.0;
   const double minSPEArea = 0.0;
   const double maxSPEArea = 0.92;
   const bool  rejectTriggerPulse = false;
   //const double gainSF = 0.35;
   const double gainSF = 1.0;
   const int holdoff_time = 1250;
   //const int holdoff_time = 0;
   
   std::string histoName = fillRate + "Histos.root";
   
   TFile *histoFile = new TFile(histoName.c_str(), "RECREATE");
   histoFile->cd();
   
   TH1F *h_widths   = new TH1F("width", "Pulse Widths", 100, 0 , 100);
   //TH1F *h_areas    = new TH1F("area", "Pulse Area", 2000, 0 , 20.0);
   TH1F *h_areas    = new TH1F("area", "Pulse Area", 1000, 0 , 80.0);
   TH1F *h_amps     = new TH1F("amp", "Pulse Amplitude", 250, 0 ,0.5);
   TH1F *h_ampOArea = new TH1F("ampoverarea", "Pulse Area / Amplitude", 100, 0 , 0.2);
   TH1F *h_timeDiff = new TH1F("timediff", "Time Between Neighboring Pulses", 250*tick_to_ns, 0, 50000*tick_to_us);
   TH1F *h_triggerDiff = new TH1F("triggerdiff", "Time Between Trigger and Subsequent Pulses", 250*tick_to_ns, 0, 50000*tick_to_us);
   TH1F *h_triggerStartTick = new TH1F("triggerstarttick", "", 3500, 2000, 5500);
   TH1F *h_triggerEndTick = new TH1F("triggerendtick", "", 3500, 2000, 5500);
   TH1F *h_triggerWidth   = new TH1F("triggerwidth", "", 1500, 0, 1500);
   TH1F *h_activeTicks    = new TH1F("activeTicks", "", 10000, 0, 50000);
   
   TH1F *h_minute = new TH1F("minute", "", 4, 0, 60);
   
   TH2F *h_areaVWidth = new TH2F("areavwidth", "Pulse Area vs Pulse Width", 30, 0, 30, 3000, 0, 30);
   TH2F *h_ampVArea   = new TH2F("ampvarea", "Pulse Amplitude vs Pulse Area", 250, 0 , 20.0, 50, 0, 0.5);
   TH2F *h_wavesVTime = new TH2F("wavesvtime", "Number of Pulses vs Day", 30, 0 , 30.0, 50, 0, 50.0);
   TH2F *h_rateVTime  = new TH2F("ratevtime", "1 PE Pulse Rate vs Day", 30, 0 , 30, 250, 0, 500.0);
   TH2F *h_volumeVTime = new TH2F("volumevtime", "Volume (cm^{3}) vs. Time", 133, 0, 1995, 54, 200, 470);
   TH2F *h_wavesVVolume = new TH2F("wavesvvolume", "Number of Pulses vs Volume", 270, 200, 470, 50, 0 , 50.0);
   
   TH1F *h_areaPerWidth[maxDepth];
   TH1F *h_ampOverAreaPerWidth[maxDepth];
   
   TH1F *h_speNumberPerDay[maxDays];
   TH1F *h_speRatePerDay[maxDays];
   
   for(int i=0; i < maxDepth; ++i){
      h_areaPerWidth[i] = new TH1F(Form("areaPerWidth%d", i),  Form("Pulse Area %d Ticks Wide", i), 250, 0 , 20.0);
      h_ampOverAreaPerWidth[i] = new TH1F(Form("ampOverAreaPerWidth%d", i),  Form("Pulse Area / Amplitude %d Ticks Wide", i), 200, 0 , 0.2);
   
   }
   
   for(int day = 0; day < maxDays; ++day){
     const int nBins = 50;
     const int nBinsRate = 250;
     h_speNumberPerDay[day] = new TH1F(Form("speNDay%d", day), "Number of SPE per Day", nBins, 0, 50.0);
     h_speRatePerDay[day]   = new TH1F(Form("speRDay%d", day), "Rate of SPE per Day", nBinsRate, 0, 500.0);
   }
   
   
   int currentWave = -1;
   double currentTime = -1.0;
   double triggerTime = 0.0;
   int triggerStartTick = 0;
   int triggerEndTick = 0; 
   double currentAmp = 0.0;
   int pulsesPerForm = 0;
   double daysSince  = -1.0;
   double currentDay = 0.0;
   
   double cumMinutes = 0.0;
   double currentMinute = -1.0;
   
   unsigned long long int totalExposureTicks = 0;
 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   bool isBigPulse = false;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      const int pulseWidth = (EndTick - StartTick);
      const double pulseArea = Area*tick_to_ns;
      const double ampOverArea = fabs(Amplitude) / pulseArea;
      double day = floor((Date % 1000000) / 10000);
      
      
      
      runTime runDate = runTime(Date);
      
      if(runDate.minute != currentMinute){
         cumMinutes += runDate.minute;
	 currentMinute = runDate.minute;
      }
      
      h_minute->Fill(runDate.minute);
      
      //lArVolume in Liters
      double lArVolume = 0.001*MinToVolume(cumMinutes);
      h_volumeVTime->Fill(cumMinutes, lArVolume);            

      const bool onePEPulse = (gainSF*pulseArea > minSPEArea && gainSF*pulseArea <= maxSPEArea);
      //const bool onePEPulse = true;
      
      //Igonre the trigger pulse cuz it's stoopid
      if(PulseNo == 1){
         triggerStartTick = StartTick;
	 triggerEndTick = EndTick;
	 h_triggerStartTick->Fill(StartTick);
	 h_triggerEndTick->Fill(EndTick);
	 h_triggerWidth->Fill(EndTick - StartTick);
	 
	 
	 
	 if(rejectTriggerPulse)
	    continue;

      }
      
      //Only look at pulses significantly after the trigger pulse
      if( (StartTick - triggerEndTick) < holdoff_time)
        continue;
      
      
      
      if(onePEPulse)
        ++pulsesPerForm;
      
      
      //indicates we've moved to a new waveform
      if(currentWave != WaveNo){
        
	if(day != currentDay){
	  daysSince = daysSince + 1.0;
	  currentDay = day;
	  cout << daysSince << " " << day << " " << runDate.day << endl;
	  //if(daysSince == 6)
	  //  cout << day << endl;
	}  
	 
        currentWave = WaveNo;
	triggerTime = StartTick*tick_to_us;
	
	currentAmp  = Amplitude;
	if(currentWave != -1) {
	  //For internal trigger; we subtracted of the holdoff time and trigger width from the readout window
	  int activeTicks = N_SAMPLES - ( holdoff_time+(triggerEndTick - triggerStartTick)+triggerStartTick );
	  h_activeTicks->Fill(activeTicks);
	  
	  //For random trigger: the entire readout window is open
	  //int activeTicks = N_SAMPLES;
	  
	  totalExposureTicks += activeTicks;
	  
	  
	  double activeTime = TICK_TO_NS*activeTicks*1e-6; 
	  double rate = (double)pulsesPerForm/activeTime; //rate in khz
	  //std::cout << rate  << std::endl;
	  h_wavesVTime->Fill(daysSince, pulsesPerForm);
	  h_rateVTime->Fill(daysSince, rate);
	  h_wavesVVolume->Fill(lArVolume, pulsesPerForm);
	  
	}  
	pulsesPerForm = 0;
	//cout << day <<endl;
      }
      
      else{
         if(onePEPulse){
	   if(currentTime > 0.0)
	      h_timeDiff->Fill((StartTick*tick_to_us - currentTime));
	   h_triggerDiff->Fill((StartTick*tick_to_us - triggerTime));
	 }
      
      }
      
      
      //cout << "All: " << onePEPulse << " " << isBigPulse << endl;
      //if(1){
      if(onePEPulse && !isBigPulse){
        	
        //cout << "Passed: " << onePEPulse << " " << isBigPulse << endl;
	h_widths->Fill(pulseWidth);	  
        h_areaVWidth->Fill( pulseWidth, fabs(gainSF*pulseArea) );
        h_amps->Fill(fabs(gainSF*Amplitude));
	
	if(fabs(Amplitude) >= minAmp && pulseWidth < maxWidth){
	  h_areaPerWidth[pulseWidth]->Fill(fabs(gainSF*pulseArea) );
	  h_ampOverAreaPerWidth[pulseWidth]->Fill( ampOverArea );
	
	}
	    
	if(fabs(gainSF*Amplitude) >= minAmp && pulseWidth > 1){
	   h_areas->Fill(fabs(gainSF*pulseArea));
	   h_ampVArea->Fill(fabs(gainSF*pulseArea), fabs(gainSF*Amplitude));
	   //cout << ampOverArea << endl;	
	   h_ampOArea->Fill( ampOverArea );
	}  
	
      }
      currentTime = StartTick*tick_to_us;
      currentAmp  = Amplitude;
   }//end of big loop
   
   for(int day = 0; day < maxDays; ++day){
     const int nBins = 50;
     const int nBinsRate = 250;
     for(int bin = 1; bin < nBins+1; ++bin){
        double binVal = h_wavesVTime->GetBinContent(day+1, bin);
	h_speNumberPerDay[day]->SetBinContent(bin,  binVal);
	h_speNumberPerDay[day]->SetBinError(bin, sqrt(binVal) ); 
     
     }
     
     for(int bin = 1; bin < nBins+1; ++bin){
        double binVal = h_rateVTime->GetBinContent(day+1, bin);
	h_speRatePerDay[day]->SetBinContent(bin,  binVal);
	h_speRatePerDay[day]->SetBinError(bin, sqrt(binVal) ); 
     
     }
            
  }
   cout << "Activity Summary: " << endl;
   cout << "Total Exposure ticks: "<< totalExposureTicks << endl;
   cout << "Total Exposure time: " << totalExposureTicks * TICK_TO_NS * 1e-9 << " seconds" << endl;
   
   
   histoFile->Write();
}//end of makeHists

void RadonFinder::MakePlots(std::string fillRate){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.1);
  gStyle->SetLabelSize(0.03, "x");
  gStyle->SetLabelSize(0.03, "y");
  gStyle->SetLineWidth(3);
  gStyle->SetLineColor(1);
  
  const double loge = 0.434;
  const double lnTwo = 0.693;
  const int maxDepth = 20;
  const int maxDays = 15;
  
  std::vector<int> good_colors;
  good_colors.clear();
  good_colors.push_back( kBlack );
  good_colors.push_back( kRed+2 );
  good_colors.push_back( kBlue+2 );
  good_colors.push_back( kYellow+2 );
  good_colors.push_back( kMagenta+2 );
  good_colors.push_back( kGreen+3 );
  good_colors.push_back( kOrange-3 );
  good_colors.push_back( kPink+2 );
  good_colors.push_back( kCyan+2 );
  
  std::string histoName = fillRate + "Histos.root";
  const string topDir = "/uboone/data/users/joelam/" + fillRate +"Plots/";
  
  //Get the Histogram File
  TFile *histoFile = new TFile(histoName.c_str(), "READ");
  TH1F *h_timeDiff  = (TH1F*)histoFile->Get("timediff");
  TH1F *h_ampOArea  = (TH1F*)histoFile->Get("ampoverarea");
  TH1F *h_amp       = (TH1F*)histoFile->Get("amp");
  TH1F *h_area      = (TH1F*)histoFile->Get("area");
  TH1F *h_width     = (TH1F*)histoFile->Get("width");
  
  TH2F *h_ampVArea  = (TH2F*)histoFile->Get("ampvarea");
  TH2F *h_pVTime    = (TH2F*)histoFile->Get("wavesvtime");
  TH2F *h_rVTime    = (TH2F*)histoFile->Get("ratevtime");
  TH2F *h_pVVolume  = (TH2F*)histoFile->Get("wavesvvolume");
  
  TH1F *h_ampOverAreaPerWidth[maxDepth];
  TH1F *h_speNumberPerDay[maxDays];
  TH1F *h_speRatePerDay[maxDays];
  
  for(int i = 0; i < maxDepth; ++i){
     h_ampOverAreaPerWidth[i] = (TH1F*)histoFile->Get(Form("ampOverAreaPerWidth%d", i));
  
  }
  
  for(int day = 0; day < maxDays; ++day){
     h_speNumberPerDay[day] = (TH1F*)histoFile->Get(Form("speNDay%d", day) );
     h_speRatePerDay[day]   = (TH1F*)histoFile->Get(Form("speNDay%d", day) );
                 
  }
  
  if(h_timeDiff == 0){
    std::cout << "Can't get Histogram" << std::endl;
    return;    
  }
  
  TH1F *h_timeDiffLog = (TH1F*)h_timeDiff->Clone("timeDiffLog");
  h_timeDiffLog->Reset();
  for(int i = 1; i < h_timeDiffLog->GetNbinsX()+1; ++i){
    
    double error;
    if(h_timeDiff->GetBinContent(i) > 0 ){
      h_timeDiffLog->SetBinContent(i, log10(h_timeDiff->GetBinContent(i)) );
      error = CalcLogError(h_timeDiff->GetBinContent(i), h_timeDiff->GetBinError(i));
    
    }
    
    else{
      h_timeDiffLog->SetBinContent(i, 0.0 );
      error = 0.0;  
    } 
    
    
    h_timeDiffLog->SetBinError(i, error); 
  
  }
  
  TCanvas *can = new TCanvas("can", "", 1500, 1500);
  h_timeDiffLog->SetLineColor(kBlack);
  h_timeDiffLog->SetMarkerStyle(kFullDotMedium);
  h_timeDiffLog->SetMarkerColor(kBlack);
  h_timeDiffLog->GetXaxis()->SetTitle("Time difference (#mu s)");
  h_timeDiffLog->Draw("E0");
  double m;
  double b;
  double m_err;
  double b_err;
  TF1 *fit = new TF1("fit", "pol1", 50.0 , 400.0 );
  //fit->SetParameters(b, m);
  //ratio->GetXaxis()->SetRangeUser(0, 0.1);
  h_timeDiffLog->Fit(fit, "QRN");
  double pars[2];
  fit->GetParameters(pars);
  double *errs = fit->GetParErrors();
  b = pars[0];
  b_err = errs[0];
  
  m = pars[1];
  m_err = (1/pow(m, 2))*pow( (errs[1]/m), 2);	
           
  std::string fitParam = Form( "y = (%.2e #pm %.2e)x + %.2f #pm %.2f", m, errs[1], b, b_err);
  std::string chi2 = Form("#chi^{2} / ndf: %.2f / %d = %.3f",fit->GetChisquare(), fit->GetNDF(), (double)fit->GetChisquare()/fit->GetNDF() );
  double tau = lnTwo/(m*loge);
  double tau_err = (m_err*lnTwo) / (loge);
  cout << "Tau error: " << tau_err << endl;
  std::string lifeTime = Form("#tau_{1/2} = %.2f #pm %.3f #mus", tau, tau_err);
  fit->SetLineColor(kRed);
  fit->Draw("same");
  TLatex *fitLabel = new TLatex(0.4, 0.8, fitParam.c_str());
  fitLabel->SetNDC();
  fitLabel->SetTextSize(0.03);
  fitLabel->SetTextColor(kRed);
  fitLabel->Draw("same");
  TLatex *chiLabel = new TLatex(0.4, 0.75, chi2.c_str());
  chiLabel->SetNDC();
  chiLabel->SetTextSize(0.03);
  chiLabel->SetTextColor(kRed);
  chiLabel->Draw("same");
  TLatex *tauLabel = new TLatex(0.4, 0.7, lifeTime.c_str());
  tauLabel->SetNDC();
  tauLabel->SetTextSize(0.03);
  tauLabel->SetTextColor(kRed);
  tauLabel->Draw("same");
  
  string outFile = topDir + "TimeDiff.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  
  h_ampOArea->SetLineColor(kBlack);
  h_ampOArea->SetTitle("Pulse Amplitude / Area");
  h_ampOArea->GetXaxis()->SetTitle("ns^{-1}");
  h_ampOArea->GetXaxis()->SetRangeUser(0, 0.14);
  h_ampOArea->SetLineWidth(3);
  h_ampOArea->Draw("hist");
  outFile = topDir + "AreaOverAmp.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  
  h_amp->SetLineColor(kBlack);
  h_amp->SetTitle("Pulse Amplitude");
  h_amp->GetXaxis()->SetTitle("V");
  h_amp->SetLineWidth(3);
  h_amp->Draw("hist");
  outFile = topDir + "Amp.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  
  h_width->SetLineColor(kBlack);
  h_width->SetTitle("Pulse Width");
  h_width->GetXaxis()->SetTitle("Number of Ticks (1 tick = 8.0 ns)");
  h_width->GetXaxis()->SetRangeUser(0, 50.0);
  h_width->SetLineWidth(3);
  h_width->Draw("hist");
  outFile = topDir + "Width.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  
  h_area->SetLineColor(kBlack);
  h_area->SetTitle("Pulse Area ");
  h_area->GetXaxis()->SetRangeUser(0, 10.0);
  h_area->GetXaxis()->SetTitle("V*ns");
  h_area->SetLineWidth(3);
  h_area->Draw("hist");
  outFile = topDir + "Area.pdf";
  
  can->Print(outFile.c_str(), "pdf");  
  can->Clear();
  can->Update();
  
  TH1D *h_areaUnitNorm = (TH1D*)h_amp->Clone("areaUnitNorm");
  //double norm = (double) 1 / h_areaUnitNorm->Integral();
  //h_areaUnitNorm->Scale(norm);
  h_areaUnitNorm->Sumw2();
  h_areaUnitNorm->SetLineColor(kBlack);
  h_areaUnitNorm->SetTitle("Pulse Area ");
  h_areaUnitNorm->GetXaxis()->SetRangeUser(0, 0.1);
  h_areaUnitNorm->GetXaxis()->SetTitle("V*ns");
  h_areaUnitNorm->SetLineWidth(3);
  //h_areaUnitNorm->SetMaximum(1.5e5);
  h_areaUnitNorm->Draw("hist");
  outFile = topDir + "AreaUnitNorm.pdf";
  TF1 f1("f1","[0]*([1]*exp(-[1])/([2]*sqrt(2*TMath::Pi()))*exp(-pow(x-[3],2)/(2*pow([2],2)))+pow([1],2)*exp(-[1])/(2*[2]*sqrt(4*TMath::Pi()))*exp(-pow(x-2*[3],2)/(4*pow([2],2))))",0.01,0.05);
  f1.SetParameters(800000.,0.02,0.03,0.02);
  //TF1 f1("f1", "pow([0], x*10)*exp(-[0]) / TMath::Gamma(x*10+1)", 0.0, 10.0);
  //double parsOne[1] = {6.0};
  //f1.SetParameters(parsOne);
  h_areaUnitNorm->Fit(&f1, "QRN");
  f1.SetLineColor(kRed);
  //f1.Draw("same");
  chi2 = Form("#chi^{2} / ndf: %.2f / %d = %.3f",f1.GetChisquare(), f1.GetNDF(), (double)f1.GetChisquare()/f1.GetNDF() );
  TLatex *chiLabel1 = new TLatex(0.4, 0.75, chi2.c_str());
  chiLabel1->SetNDC();
  chiLabel1->SetTextSize(0.03);
  chiLabel1->SetTextColor(kRed);
  //chiLabel1->Draw("same");  
  can->Print(outFile.c_str(), "pdf");  
  can->Clear();
  can->Update();
  
  
  h_area->GetXaxis()->SetRangeUser(0, 80.0);
  h_area->Draw("hist");
  gPad->SetLogy(1);
  outFile = topDir + "AreaLog.pdf";
  can->Print(outFile.c_str(), "pdf");
  gPad->SetLogy(0);
  
  can->Clear();
  can->Update();
  
  TGaxis::SetMaxDigits(3);
  //SetRedHeatPalette();
  
  h_ampVArea->GetYaxis()->SetTitle("Pulse Amplitude");
  h_ampVArea->GetYaxis()->SetTitleOffset(1.35);
  h_ampVArea->GetXaxis()->SetTitle("Pulse Area");
  h_ampVArea->Draw("colrz");
  gPad->SetLogz();
  gPad->Update();
  TPaletteAxis *palette = (TPaletteAxis*)h_ampVArea->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(0.85);
  palette->SetX2NDC(0.9);
  
  outFile = topDir + "AmpVArea.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  
  h_ampVArea->GetYaxis()->SetTitle("Pulse Amplitude");
  h_ampVArea->GetXaxis()->SetTitle("Pulse Area");
  h_ampVArea->GetXaxis()->SetRangeUser(0, 1.0);
  h_ampVArea->Draw("colrz");
  gPad->Update();
  gPad->SetLogz(1);
  TPaletteAxis *palette0 = (TPaletteAxis*)h_ampVArea->GetListOfFunctions()->FindObject("palette");
  palette0->SetX1NDC(0.85);
  palette0->SetX2NDC(0.9);
  
  outFile = topDir + "AmpVAreaZoom.pdf";
  can->Print(outFile.c_str(), "pdf");  
  
  can->Clear();
  can->Update();
  
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.4);
  legend->SetHeader("Pulse Width:");
  
  for(int i = 2; i < 7; ++i){
    h_ampOverAreaPerWidth[i]->SetLineWidth(3);
    h_ampOverAreaPerWidth[i]->SetTitle("Pulse Amplitude / Area");
    h_ampOverAreaPerWidth[i]->SetLineColor(good_colors[i-2]);
    legend->AddEntry(h_ampOverAreaPerWidth[i], Form("%d ticks wide", i), "l");
    if(i == 2)
      h_ampOverAreaPerWidth[i]->Draw("hist");
    else  
      h_ampOverAreaPerWidth[i]->Draw("hist same");
  
  }
  legend->Draw("same");
  outFile = topDir + "AmpOverAreaBins.pdf";
  can->Print(outFile.c_str(), "pdf"); 

  can->Clear();
  can->Update();
  
  h_pVTime->GetYaxis()->SetTitle("Number of 1 PE Pulses");
  TProfile *prof_pVTime = h_pVTime->ProfileX("prof_tmp"); // gotta do this before scaling
  h_pVTime->GetYaxis()->SetRangeUser(0, 50.0);
  h_pVTime->GetXaxis()->SetTitle("Days Since Start of Data Taking");
  h_pVTime->SetMaximum(11e3);
  //normalize to unit area
  double scale = (double)(1 / h_pVTime->Integral() );
  
  
  //h_pVTime->SetMinimum(5e3);
  //h_pVTime->GetXaxis()->SetRangeUser(0, 1.0);
  h_pVTime->Scale(scale);
  h_pVTime->Draw("colrz");
  gPad->SetLogz(1);
  gPad->Update();
  palette0 = (TPaletteAxis*)h_pVTime->GetListOfFunctions()->FindObject("palette");
  
  if(palette0 != 0){
    palette0->SetX1NDC(0.85);
    palette0->SetX2NDC(0.9);
  }

  outFile = topDir + "1PEvTime.pdf";
  can->Print(outFile.c_str(), "pdf");   
  
  can->Clear();
  can->Update();  
  
  prof_pVTime->GetYaxis()->SetTitle("Mean Number of 1 PE Pulses");

  prof_pVTime->GetYaxis()->SetRangeUser(0.0, 20.0);
  prof_pVTime->GetXaxis()->SetRangeUser(0.0, 12.0);
  prof_pVTime->GetXaxis()->SetTitle("Days Since Start of Data Taking");
  prof_pVTime->SetMarkerStyle(kFullDotLarge);
  prof_pVTime->SetMarkerSize(1.75);
  prof_pVTime->SetMarkerColor(kBlack);
  prof_pVTime->SetLineColor(kBlack);
  prof_pVTime->Draw();
  outFile = topDir + "1PEProf.pdf";
  can->Print(outFile.c_str(), "pdf"); 

  h_rVTime->GetYaxis()->SetTitle("Rate of 1 PE Pulses");
  TProfile *prof_rVTime = h_rVTime->ProfileX("prof_tmp"); //gotta do this before scaling
  h_rVTime->GetYaxis()->SetRangeUser(0, 120.0);
  h_rVTime->GetXaxis()->SetTitle("Days Since Start of Data Taking");
  //normalize to unit area
  scale = (double)(1 / h_rVTime->Integral() );
  
  //h_pVTime->SetMinimum(5e3);
  //h_pVTime->GetXaxis()->SetRangeUser(0, 1.0);
  h_rVTime->Scale(scale);
  h_rVTime->Draw("colrz");
  gPad->SetLogz(1);
  gPad->Update();
  palette0 = (TPaletteAxis*)h_rVTime->GetListOfFunctions()->FindObject("palette");
  
  if(palette0 != 0){
    palette0->SetX1NDC(0.85);
    palette0->SetX2NDC(0.9);
  }

  outFile = topDir + "1PERatevTime.pdf";
  can->Print(outFile.c_str(), "pdf");   
  
  can->Clear();
  can->Update();  
  
  
  prof_rVTime->GetYaxis()->SetTitle("Mean Rate of 1 PE Pulses (kHz)");
  
  prof_rVTime->GetYaxis()->SetRangeUser(0.0, 40.0);
  prof_rVTime->GetXaxis()->SetRangeUser(0.0, 12.0);
  prof_rVTime->GetXaxis()->SetTitle("Days Since Start of Data Taking");
  prof_rVTime->SetMarkerStyle(kFullDotLarge);
  prof_rVTime->SetMarkerSize(1.75);
  prof_rVTime->SetMarkerColor(kBlack);
  prof_rVTime->SetLineColor(kBlack);
  prof_rVTime->Draw();
  outFile = topDir + "1PERateProf.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.4);
  legend->SetHeader("Pulse Width:");
  
  for(int i = 0; i < 8; ++i){
    
    h_speNumberPerDay[i]->SetLineWidth(3);
    double scaleSF = 1 / h_speNumberPerDay[i]->Integral();
    h_speNumberPerDay[i]->Scale(scaleSF);
    h_speNumberPerDay[i]->SetMaximum(0.23);
    h_speNumberPerDay[i]->GetXaxis()->SetRangeUser(0, 15.0);
    h_speNumberPerDay[i]->SetTitle("Pulse Amplitude / Area");
    h_speNumberPerDay[i]->SetLineColor(good_colors[i]);
    h_speNumberPerDay[i]->SetMarkerColor(good_colors[i]);
    h_speNumberPerDay[i]->SetMarkerStyle(20);
    legend1->AddEntry(h_speNumberPerDay[i], Form("Day %d", i), "l");
    if(i == 0)
      h_speNumberPerDay[i]->Draw("hist");
    else  
      h_speNumberPerDay[i]->Draw("hist same");
  
  }
  legend1->Draw("same");
  outFile = topDir + "NPEBins.pdf";
  can->Print(outFile.c_str(), "pdf"); 

  can->Clear();
  can->Update();
  
  h_pVVolume->GetYaxis()->SetTitle("Number of 1 PE Pulses");
  TProfile *prof_pVVolume = h_pVVolume->ProfileX("prof_tmp0"); //gotta do this before scaling
 // h_rVTime->GetYaxis()->SetRangeUser(0, 600.0);
  h_pVVolume->GetXaxis()->SetTitle("LAr in TallBo (L)");
  
  h_pVVolume->Draw("colrz");
  gPad->SetLogz(1);
  gPad->Update();
  palette0 = (TPaletteAxis*)h_pVVolume->GetListOfFunctions()->FindObject("palette");
  
  if(palette0 != 0){
    palette0->SetX1NDC(0.85);
    palette0->SetX2NDC(0.9);
  }

  outFile = topDir + "1PENumbervVolume.pdf";
  can->Print(outFile.c_str(), "pdf");   
  
  can->Clear();
  can->Update();
  
  prof_pVVolume->GetYaxis()->SetTitle("Mean Number of 1 PE Pulses");
  
  prof_pVVolume->GetYaxis()->SetRangeUser(0.0, 12.0);
  prof_pVVolume->GetXaxis()->SetRangeUser(200.0, 512.0);
  prof_pVVolume->GetXaxis()->SetTitle("LAr in TallBo (L)");
  prof_pVVolume->SetMarkerStyle(kFullDotLarge);
  prof_pVVolume->SetMarkerSize(1.75);
  prof_pVVolume->SetMarkerColor(kBlack);
  prof_pVVolume->SetLineColor(kBlack);
  prof_pVVolume->Draw();
  outFile = topDir + "1PENumbervVolumeProf.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  
     

}//end of MakePlots

void RadonFinder::MakeComparisonPlots(bool areaNorm){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.1);
  gStyle->SetLabelSize(0.03, "x");
  gStyle->SetLabelSize(0.03, "y");
  gStyle->SetLineWidth(3);
  gStyle->SetLineColor(1);
  
  const double loge = 0.434;
  const double lnTwo = 0.693;
  const int maxDepth = 20;
  const int maxDays = 15;
  const double exposure_factor = 0.877; 
  
  std::vector<int> good_colors;
  good_colors.clear();
  good_colors.push_back( kBlack );
  good_colors.push_back( kRed+2 );
  good_colors.push_back( kBlue+2 );
  good_colors.push_back( kYellow+2 );
  good_colors.push_back( kMagenta+2 );
  good_colors.push_back( kGreen+3 );
  good_colors.push_back( kOrange-3 );
  good_colors.push_back( kPink+2 );
  good_colors.push_back( kCyan+2 );
  
  //Get the Histogram Files
  TFile *histoFileSlow = new TFile("SlowPRandomHistos.root", "READ");
  TFile *histoFileFast = new TFile("FastPRandomHistos.root", "READ");   
  
  TH1F *h_timeDiff_slow    = (TH1F*)histoFileSlow->Get("timediff");
  TH1F *h_ampOArea_slow    = (TH1F*)histoFileSlow->Get("ampoverarea");
  TH1F *h_amp_slow         = (TH1F*)histoFileSlow->Get("amp");
  TH1F *h_area_slow        = (TH1F*)histoFileSlow->Get("area");
  TH1F *h_width_slow       = (TH1F*)histoFileSlow->Get("width");
  
  TH2F *h_ampVArea_slow    = (TH2F*)histoFileSlow->Get("ampvarea");
  TH2F *h_timeVPulse_slow  = (TH2F*)histoFileSlow->Get("wavesvtime");
  TH2F *h_rateVPulse_slow  = (TH2F*)histoFileSlow->Get("ratevtime");
  
  TH1F *h_timeDiff_fast    = (TH1F*)histoFileFast->Get("timediff");
  TH1F *h_ampOArea_fast    = (TH1F*)histoFileFast->Get("ampoverarea");
  TH1F *h_amp_fast         = (TH1F*)histoFileFast->Get("amp");
  TH1F *h_area_fast        = (TH1F*)histoFileFast->Get("area");
  TH1F *h_width_fast       = (TH1F*)histoFileFast->Get("width");
  
  TH2F *h_ampVArea_fast    = (TH2F*)histoFileFast->Get("ampvarea");
  TH2F *h_timeVPulse_fast  = (TH2F*)histoFileFast->Get("wavesvtime");
  TH2F *h_rateVPulse_fast  = (TH2F*)histoFileFast->Get("ratevtime");
  
  TH1F *h_speNumberPerDay_fast[maxDays];
  TH1F *h_speNumberPerDay_slow[maxDays];
  
  for(int day = 0; day < maxDays; ++day){
     h_speNumberPerDay_fast[day]   = (TH1F*)histoFileFast->Get(Form("speNDay%d", day) );
     h_speNumberPerDay_slow[day]   = (TH1F*)histoFileSlow->Get(Form("speNDay%d", day) );
                 
  }
  
  //area normalize the two runs
  int startBin = h_amp_slow->GetXaxis()->FindBin(0.01);
  double scale = 1.0;
  string normString;
  if(areaNorm){
    scale = ( h_amp_fast->Integral(startBin, h_amp_fast->GetNbinsX()-1) / h_amp_slow->Integral(startBin, h_amp_fast->GetNbinsX()-1) );
    normString = "Area Normalized";
    
  }
    
  else{
    scale = exposure_factor;  
    normString = "Exposure Normalized";
  }  
  
  cout << scale << endl;
  
  const string topDir = "/uboone/data/users/joelam/ComparisonPlots/";
  
  TLatex *normLabel = new TLatex(0.6, 0.7, normString.c_str());
  normLabel->SetNDC();
  normLabel->SetTextSize(0.03);
  normLabel->SetTextColor(kBlue-2);
  
  TCanvas *can = new TCanvas("can", "", 1500, 1500);
  h_amp_slow->SetLineColor(kBlack);
  h_amp_fast->SetLineColor(kRed+2);
  h_amp_slow->GetXaxis()->SetTitle("V");
  h_amp_slow->SetLineWidth(3);
  h_amp_fast->SetLineWidth(3);
  h_amp_slow->Scale(scale);
  h_amp_slow->Draw("hist");
  h_amp_fast->Draw("hist same");
  TLegend *legend = new TLegend(0.7, 0.85, 0.9, 0.75);
  legend->AddEntry(h_amp_fast, "Fast Fill", "l");
  legend->AddEntry(h_amp_slow, "Slow Fill", "l");
  legend->Draw("same");

  normLabel->Draw("same");  
  string outFile = topDir + "AmpComparison.pdf";  
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  
  
  
  h_area_slow->SetLineColor(kBlack);
  h_area_fast->SetLineColor(kRed+2);
  h_area_slow->GetXaxis()->SetTitle("V*ns");
  h_area_slow->SetTitle("Pulse Area");
  h_area_fast->GetXaxis()->SetTitle("V*ns");
  h_area_fast->SetTitle("Pulse Area");
  h_area_slow->SetLineWidth(3);
  h_area_fast->SetLineWidth(3);
  
  h_area_slow->GetXaxis()->SetRangeUser(0, 5.0);
  h_area_fast->GetXaxis()->SetRangeUser(0, 5.0);  
  h_area_slow->Scale(scale);
  h_area_slow->SetMaximum(500e3);
  h_area_slow->Draw("hist");
  h_area_fast->Draw("hist same");
  legend->Draw("same");
  outFile = topDir + "AreaComparison.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  
  h_area_slow->GetXaxis()->SetRangeUser(0, 80.0);
  h_area_fast->GetXaxis()->SetRangeUser(0, 80.0);
  h_area_slow->Scale(scale);
  h_area_slow->Draw("hist");
  h_area_fast->Draw("hist same");
  gPad->SetLogy(1);
  legend->Draw("same");
    normLabel->Draw("same");
  outFile = topDir + "AreaComparisonLog.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  gPad->SetLogy(0);
  
  h_ampOArea_slow->SetLineColor(kBlack);
  h_ampOArea_fast->SetLineColor(kRed+2);
  h_ampOArea_slow->GetXaxis()->SetTitle("ns^{-1}");
  h_ampOArea_slow->SetTitle("Pulse Amplitude / Area");
  h_ampOArea_fast->GetXaxis()->SetTitle("ns");
  h_ampOArea_fast->SetTitle("Pulse Amplitude / Area");
  h_ampOArea_slow->SetLineWidth(3);
  h_ampOArea_fast->SetLineWidth(3);

  h_ampOArea_slow->GetXaxis()->SetRangeUser(0, 0.15);
  h_ampOArea_fast->GetXaxis()->SetRangeUser(0, 0.15);
  h_ampOArea_slow->Scale(scale);
  h_ampOArea_slow->Draw("hist");
  h_ampOArea_fast->Draw("hist same");
  legend->Draw("same");
    normLabel->Draw("same");
  outFile = topDir + "AmpOverAreaCompare.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  
  h_width_slow->SetLineColor(kBlack);
  h_width_fast->SetLineColor(kRed+2);
  h_width_slow->GetXaxis()->SetTitle("ns");
  h_width_slow->SetTitle("Pulse Width");
  h_width_fast->GetXaxis()->SetTitle("ns");
  h_width_fast->SetTitle("Pulse Width");
  h_width_slow->SetLineWidth(3);
  h_width_fast->SetLineWidth(3);  
  h_width_slow->GetXaxis()->SetRangeUser(0, 50.0);
  h_width_fast->GetXaxis()->SetRangeUser(0, 50.0);
  h_width_slow->Scale(scale);
  h_width_slow->Draw("hist");
  h_width_fast->Draw("hist same");
  legend->Draw("same");
    normLabel->Draw("same");
  outFile = topDir + "WidthCompare.pdf";
  can->Print(outFile.c_str(), "pdf");

  can->Clear();
  can->Update();
  
  TProfile *prof_slow = h_timeVPulse_slow->ProfileX("profSlow");
  TProfile *prof_fast = h_timeVPulse_fast->ProfileX("profFast");
  prof_slow->SetMarkerColor(kBlack);
  prof_fast->SetMarkerColor(kRed+2);
  prof_slow->SetLineColor(kBlack);
  prof_fast->SetLineColor(kRed+2);
  prof_slow->SetMarkerStyle(kFullDotLarge);
  prof_fast->SetMarkerStyle(kFullDotLarge);
  prof_fast->SetMarkerSize(1.08);
  prof_slow->SetMarkerSize(1.08);
  
  prof_slow->GetXaxis()->SetTitle("Days Since Start of Data Taking");
  prof_fast->GetXaxis()->SetTitle("Days Since Start of Data Taking");
  prof_fast->GetYaxis()->SetTitle("Mean Number of 1 PE Pulses");
  prof_fast->GetYaxis()->SetRangeUser(0.0, 12.0);
  prof_fast->GetXaxis()->SetRangeUser(0.0, 12.0);
  
  prof_fast->SetMarkerSize(1.75);
  prof_slow->SetMarkerSize(1.75);
  prof_fast->Draw("");
  prof_slow->Draw("same");
  outFile = topDir +"SPENumberCompare.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  delete prof_slow;
  delete prof_fast;
  
  prof_slow = h_rateVPulse_slow->ProfileX("profSlow");
  prof_fast = h_rateVPulse_fast->ProfileX("profFast");
  prof_slow->SetMarkerColor(kBlack);
  prof_fast->SetMarkerColor(kRed+2);
  prof_slow->SetLineColor(kBlack);
  prof_fast->SetLineColor(kRed+2);
  prof_slow->SetMarkerStyle(kFullDotLarge);
  prof_fast->SetMarkerStyle(kFullDotLarge);
  prof_fast->SetMarkerSize(1.08);
  prof_slow->SetMarkerSize(1.08);
  
  prof_slow->GetXaxis()->SetTitle("Days Since Start of Data Taking");
  prof_fast->GetXaxis()->SetTitle("Days Since Start of Data Taking");
  prof_fast->GetYaxis()->SetTitle("Mean 1 PE Rate (kHz)");
  prof_fast->GetYaxis()->SetRangeUser(0.0, 40.0);
  prof_fast->GetXaxis()->SetRangeUser(0.0, 12.0);
  
  prof_fast->SetMarkerSize(1.75);
  prof_slow->SetMarkerSize(1.75);
  prof_fast->Draw("");
  prof_slow->Draw("same");
  outFile = topDir +"SPERateCompare.pdf";
  can->Print(outFile.c_str(), "pdf");
  
  can->Clear();
  can->Update();
  
  for(int day = 0; day < maxDays; ++day){
     double scaleSF_fast = 1 / h_speNumberPerDay_fast[day]->Integral();
     double scaleSF_slow = 1 / h_speNumberPerDay_slow[day]->Integral();
     //h_speNumberPerDay_fast[day]->Scale(scaleSF);
     
     //scaleSF = 1 /  h_speNumberPerDay_slow[day]->Integral();
     TF1 f1_slow("f1Slow", "TMath::Poisson(x, [0])", 0.0, 25.0);
     TF1 f1_fast("f1Fast", "TMath::Poisson(x, [0])", 0.0, 25.0);
     double parsOne[1] = {2.0};
     f1_fast.SetParameters(parsOne);
     f1_slow.SetParameters(parsOne);
     
     h_speNumberPerDay_slow[day]->Scale(scaleSF_slow);
     h_speNumberPerDay_fast[day]->Scale(scaleSF_fast);
   
     h_speNumberPerDay_fast[day]->SetMarkerStyle(kFullDotMedium);
     h_speNumberPerDay_slow[day]->SetMarkerStyle(kFullDotMedium);
     
     h_speNumberPerDay_fast[day]->SetMarkerColor(kRed+2);
     h_speNumberPerDay_slow[day]->SetMarkerColor(kBlack);
     
     h_speNumberPerDay_fast[day]->SetLineColor(kRed+2);
     h_speNumberPerDay_slow[day]->SetLineColor(kBlack);
     
     h_speNumberPerDay_fast[day]->SetMarkerSize(3);
     h_speNumberPerDay_slow[day]->SetMarkerSize(3);
     
     //h_speNumberPerDay_slow[day]->SetMaximum(0.22);
     
     h_speNumberPerDay_slow[day]->GetXaxis()->SetRangeUser(0, 25.0);
     h_speNumberPerDay_fast[day]->GetXaxis()->SetRangeUser(0, 25.0);
     
     h_speNumberPerDay_slow[day]->GetYaxis()->SetRangeUser(0, 0.3);
     h_speNumberPerDay_fast[day]->GetYaxis()->SetRangeUser(0, 0.3);
     
     h_speNumberPerDay_fast[day]->SetTitle(Form("Number of SPE Day %d", day));
     h_speNumberPerDay_fast[day]->GetXaxis()->SetTitle("Number of Single PE Pulses per Readout");
     //h_speNumberPerDay_fast[day]->GetYaxis()->SetTitle("Number Readouts");
     
     h_speNumberPerDay_fast[day]->Draw("PE");
     h_speNumberPerDay_slow[day]->Draw("P same");
     
     h_speNumberPerDay_fast[day]->Fit(&f1_fast, "QRN");
     h_speNumberPerDay_slow[day]->Fit(&f1_slow, "QRN");
     
     f1_fast.SetLineColor(kRed+2);
     f1_slow.SetLineColor(kBlack);
     
     f1_fast.Draw("same");
     f1_slow.Draw("same");
     
     
     legend->Draw("same");
     
     outFile = topDir + Form("NSPEDay%d.pdf", day);
     can->Print(outFile.c_str(), "pdf");
  
  
  } 
  


}//end of comparison plots


double RadonFinder::CalcLogError(double val, double err){
   const double lnTen = log(10);
   return (fabs(err/val*lnTen));

}

void RadonFinder::SetRedHeatPalette(){
  const int NRGBs = 7;
  static bool initialized=false;
  const unsigned int n_color_contours = 999;
  static int* colors=new int[n_color_contours];

  if(!initialized){
    // White -> red
    Double_t stops[NRGBs] = { 0.0, 0.16, 0.32, 0.48, 0.64, 0.8, 1.000};
    Double_t red[NRGBs]   = { 1.00, 1.00, 0.99, 0.84, 0.70, 0.55, 0.20 };
    Double_t green[NRGBs] = { 0.96, 0.88, 0.42, 0.23, 0.09, 0.06, 0.00 };
    Double_t blue[NRGBs]  = { 0.94, 0.82, 0.29, 0.17, 0.11, 0.08, 0.05 };
    int colmin=TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, n_color_contours);
    for(uint i=0; i<n_color_contours; ++i) colors[i]=colmin+i;

    initialized=true;
  }
  gStyle->SetNumberContours(n_color_contours);
  gStyle->SetPalette(n_color_contours, colors);
}

double RadonFinder::MinToHeight(double minutes){
   
   if(minutes > 1909.0)
     minutes = 1909.0;
   //third degree polynomial fit  
   return 2.54*(2.0e-8*pow(minutes, 3) - 5.0e-5*pow(minutes, 2) + 0.0431*minutes + 35.289);
   

}

double RadonFinder::MinToVolume(double minutes){
   const double tallBoRadius = (56.0 / 2.0); //cm
   const double pi = 3.14159;
   return pi*pow(tallBoRadius, 2)*MinToHeight(minutes);

}

/*
runTime_t RadonFinder::DecodeDate(){
    runTime_t rval;
    double year = floor((Date % 1000000000000) / 100000000 );
    //cout << "Year: " << year << endl;
    double month = floor((Date % 100000000) / 1000000 );
    //cout << "Month: " << month << endl;
    double day = floor((Date % 1000000) / 10000);
    //cout << "Day : " << day << endl;
    double hour = floor((Date % 10000) / 100);
    //cout << "Date: " << Date << endl;
    //cout << "Hour : " << hour << endl;
    double minute = floor((Date % 100) );
    //cout << "Minute: " << minute << endl;
    
    rval.year = year;
    rval.month = month;
    rval.day = day;
    rval.hour = hour
    rval.minute = minute;
    
    return rval;
    
    
    
    //rval.day = int(day);


}*/


/*void RadonFinder::GetSPERateHist(TProfile *prof, TH1F *returnHist, int holdOffTime){
    
    double binVal = 0.0;
    double binErr = 0.0;
    const double denom = (N_SAMPLES-holdOffTime)*TICK_TO_NS;  
    for(int bin=1; bin < prof->GetNbinsX()+2; ++bin){
       binVal = prof->GetBinContent(bin) / denom;
       returnHist->SetBinContent(bin, binVal);
    
    }


}*/
