#define LEDHistos_cxx
#include "LEDHistos.h"
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



void LEDHistos::MakeHists()
{
 
   const double tick_to_ns = 0.8;
   const double tick_to_us = 0.8*1e-3;
   const double tick_to_s = 0.8*1e-9;
   const int maxDepth = 20;
   const double minAmp = 0.01;
   const int maxWidth = 5;
   const double startDay = 13.0;
   const int maxTicks = 2500;
   
   
   TFile *histoFile = new TFile("HistogramsLED.root", "RECREATE");
   histoFile->cd();
   
   TH1F *h_widths   = new TH1F("width", "Pulse Widths", 100, 0 , 100);
   TH1F *h_pedestal   = new TH1F("pedestal", "Pedestal", 100, -0.0225 , 0.0225);
   //TH1F *h_areas    = new TH1F("area", "Pulse Area", 2000, 0 , 20.0);
   TH1F *h_areas    = new TH1F("area", "Pulse Area", 1000, 0 , 40.0);
   TH1F *h_amps     = new TH1F("amp", "Pulse Amplitude", 450, 0 ,0.225);
   TH1F *h_ampOArea = new TH1F("ampoverarea", "Pulse Area / Amplitude", 200, 0 , 0.2);
   TH1F *h_timeDiff = new TH1F("timediff", "Time Between Neighboring Pulses", 250*tick_to_ns, 0, 50000*tick_to_us);
   TH1F *h_triggerDiff = new TH1F("triggerdiff", "Time Between Trigger and Subsequent Pulses", 250*tick_to_ns, 0, 50000*tick_to_us);
   
   TH2F *h_areaVWidth = new TH2F("areavwidth", "Pulse Area vs Pulse Width", 30, 0, 30, 3000, 0, 30);
   TH2F *h_ampVArea   = new TH2F("ampvarea", "Pulse Amplitude vs Pulse Area", 250, 0 , 20.0, 50, 0, 0.5);
   TH2F *h_wavesVTime = new TH2F("wavesvtime", "Number of Pulses vs Day", 30, 0 , 30.0, 50, 0, 50.0);
   
   TH1F *h_areaPerWidth[maxDepth];
   TH1F *h_ampOverAreaPerWidth[maxDepth];
   
   for(int i=0; i < maxDepth; ++i){
      h_areaPerWidth[i] = new TH1F(Form("areaPerWidth%d", i),  Form("Pulse Area %d Ticks Wide", i), 250, 0 , 20.0);
      h_ampOverAreaPerWidth[i] = new TH1F(Form("ampOverAreaPerWidth%d", i),  Form("Pulse Area / Amplitude %d Ticks Wide", i), 200, 0 , 0.2);
   
   }
   
   int currentWave = -1;
   double currentTime = 0.0;
   double triggerTime = 0.0; 
   double currentAmp = 0.0;
   int pulsesPerForm = 0;
 
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      const int pulseWidth = (EndTick - StartTick);
      const double pulseArea = Area*tick_to_ns;
      const double ampOverArea = Amplitude / pulseArea;
      double day = floor((Date % 1000000) / 10000);
      day =( day - startDay);
      
      /*if(Date < 2.016e11){
        day = (day - startDay);
      
      }
      
      else{
        day = (31.0 - startDay) + day;
      }*/
      //cout << day << endl;
      ++pulsesPerForm;
      
      for(int i=0; i < maxTicks; ++i){
          if(Pedestal[i] != -10.0)
	    h_pedestal->Fill(-Pedestal[i]);
      
      }
      
      
      //indicates we've moved to a new waveform
      if(currentWave != WaveNo){
         
        currentWave = WaveNo;
	triggerTime = StartTick*tick_to_us;
	currentTime = StartTick*tick_to_us;
	currentAmp  = Amplitude;
	if(currentWave != -1) 
	  h_wavesVTime->Fill(day, pulsesPerForm);
	pulsesPerForm = 0;
	//cout << day <<endl;
      }
      
      else{
         if(fabs(Amplitude) >= minAmp && fabs(currentAmp) >= minAmp && pulseWidth < maxWidth){
	   h_timeDiff->Fill((StartTick*tick_to_us - currentTime));
	   h_triggerDiff->Fill((StartTick*tick_to_us - triggerTime));
	 }
      
      }
      
      if(fabs(Amplitude) > minAmp){
        h_amps->Fill(fabs(Amplitude));
        h_areas->Fill(fabs(pulseArea*tick_to_ns));
      
      }
      
      if(StartTick != -1 && EndTick != -1){
        	
        h_widths->Fill(pulseWidth);	  
        h_areaVWidth->Fill( pulseWidth, fabs(pulseArea) );
        
	
	if(pulseWidth < maxDepth && fabs(Amplitude) >= minAmp){
	  h_areaPerWidth[pulseWidth]->Fill(fabs(pulseArea) );
	  h_ampOverAreaPerWidth[pulseWidth]->Fill( ampOverArea );
	
	}
	    
	if(Area != 0 && fabs(Amplitude) >= minAmp){
	   
	   h_ampVArea->Fill(fabs(pulseArea), fabs(Amplitude));	
	   h_ampOArea->Fill( ampOverArea );
	}  
	
      }
      currentTime = StartTick*tick_to_us;
      currentAmp  = Amplitude;
   }//end of big loop
   histoFile->Write();
}//end of makeHists

void LEDHistos::MakePlots(){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.1);
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetLineWidth(3);
  gStyle->SetLineColor(1);
  
  const double loge = 0.434;
  const double lnTwo = 0.693;
  const int maxDepth = 20;
  
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
  
  //Get the Histogram File
  TFile *histoFile = new TFile("HistogramsLED.root", "READ");
  TH1F *h_timeDiff  = (TH1F*)histoFile->Get("timediff");
  TH1F *h_ampOArea  = (TH1F*)histoFile->Get("ampoverarea");
  TH1F *h_amp       = (TH1F*)histoFile->Get("amp");
  TH1F *h_area      = (TH1F*)histoFile->Get("area");
  TH1F *h_width     = (TH1F*)histoFile->Get("width");
  TH1F *h_ped     = (TH1F*)histoFile->Get("pedestal");
  TH2F *h_ampVArea  = (TH2F*)histoFile->Get("ampvarea");
  TH1F *h_ampOverAreaPerWidth[maxDepth];
  
  for(int i = 0; i < maxDepth; ++i){
     h_ampOverAreaPerWidth[i] = (TH1F*)histoFile->Get(Form("ampOverAreaPerWidth%d", i));
  
  }
  
  if(h_timeDiff == 0){
    std::cout << "Can't get Histogram" << std::endl;
    return;    
  }
  
  
  TCanvas *can = new TCanvas("can", "", 1500, 1500);
  
  
  can->Clear();
  can->Update();
  
  h_ampOArea->SetLineColor(kBlack);
  h_ampOArea->SetTitle("Pulse Amplitude / Area");
  h_ampOArea->GetXaxis()->SetTitle("ns^{-1}");
  h_ampOArea->GetXaxis()->SetRangeUser(0, 0.14);
  h_ampOArea->SetLineWidth(3);
  h_ampOArea->Draw("hist");
  can->Print("AreaOverAmp.pdf", "pdf");
  
  can->Clear();
  can->Update();
  
  h_amp->SetLineColor(kBlack);
  h_amp->SetTitle("Pulse Amplitude");
  h_amp->GetXaxis()->SetTitle("V");
  h_amp->SetLineWidth(3);
  h_amp->Draw("hist");
  gPad->SetLogy(0);
  can->Print("Amp.pdf", "pdf");
  
  can->Clear();
  can->Update();
  
  h_ped->SetLineColor(kRed);
  h_ped->SetTitle("Pulse Amplitude");
  h_ped->GetXaxis()->SetTitle("V");
  //h_ped->GetXaxis()->SetRangeUser(0, 0.12);
  h_ped->SetLineWidth(3);
  h_ped->Draw("hist");
  h_amp->Draw("hist same");
  gPad->SetLogy(1);
  can->Print("AmpAndPed.pdf", "pdf");
  
  can->Clear();
  can->Update();
  
  h_width->SetLineColor(kBlack);
  h_width->SetTitle("Pulse Width");
  h_width->GetXaxis()->SetTitle("Number of Ticks (1 tick = 8.0 ns)");
  h_width->GetXaxis()->SetRangeUser(0, 20.0);
  h_width->SetLineWidth(3);
  h_width->Draw("hist");
  can->Print("Width.pdf", "pdf");
  
  can->Clear();
  can->Update();
  
  h_area->SetLineColor(kBlack);
  h_area->SetTitle("Pulse Area ");
  h_area->GetXaxis()->SetRangeUser(0, 4.0);
  h_area->GetXaxis()->SetTitle("V*ns");
  h_area->SetLineWidth(3);
  h_area->Draw("hist");
  can->Print("Area.pdf", "pdf");
  
  can->Clear();
  can->Update();
  
  
  h_area->GetXaxis()->SetRangeUser(0, 4.0);
  h_area->Draw("hist");
  gPad->SetLogy(1);
  can->Print("AreaLog.pdf", "pdf");
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
  
  
  can->Print("AmpVArea.pdf", "pdf");
  
  can->Clear();
  can->Update();
  
  h_ampVArea->GetYaxis()->SetTitle("Pulse Amplitude");
  h_ampVArea->GetXaxis()->SetTitle("Pulse Area");
  h_ampVArea->GetXaxis()->SetRangeUser(0, 6.0);
  h_ampVArea->Draw("colrz");
  gPad->Update();
  gPad->SetLogz(1);
  TPaletteAxis *palette0 = (TPaletteAxis*)h_ampVArea->GetListOfFunctions()->FindObject("palette");
  palette0->SetX1NDC(0.85);
  palette0->SetX2NDC(0.9);
  
  
  can->Print("AmpVAreaZoom.pdf", "pdf");
  
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
  can->Print("AmpOverAreaBins.pdf", "pdf");
  
  
  
  
     

}//end of MakePlots

void LEDHistos::MakeComparisonPlots(){

  gStyle->SetOptStat(0);
  gStyle->SetPalette(1);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.1);
  gStyle->SetLabelSize(0.05, "x");
  gStyle->SetLabelSize(0.05, "y");
  gStyle->SetLineWidth(3);
  gStyle->SetLineColor(1);
  
  const double loge = 0.434;
  const double lnTwo = 0.693;
  const int maxDepth = 20;
  
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
  TFile *histoFileSlow = new TFile("HistogramsSlowFill.root", "READ");
  TFile *histoFileFast = new TFile("HistogramsFastFill.root", "READ");   
  
  TH1F *h_timeDiff_slow  = (TH1F*)histoFileSlow->Get("timediff");
  TH1F *h_ampOArea_slow  = (TH1F*)histoFileSlow->Get("ampoverarea");
  TH1F *h_amp_slow       = (TH1F*)histoFileSlow->Get("amp");
  TH1F *h_area_slow      = (TH1F*)histoFileSlow->Get("area");
  TH1F *h_width_slow     = (TH1F*)histoFileSlow->Get("width");
  TH2F *h_ampVArea_slow  = (TH2F*)histoFileSlow->Get("ampvarea");
  
  TH1F *h_timeDiff_fast  = (TH1F*)histoFileFast->Get("timediff");
  TH1F *h_ampOArea_fast  = (TH1F*)histoFileFast->Get("ampoverarea");
  TH1F *h_amp_fast       = (TH1F*)histoFileFast->Get("amp");
  TH1F *h_area_fast      = (TH1F*)histoFileFast->Get("area");
  TH1F *h_width_fast     = (TH1F*)histoFileFast->Get("width");
  TH2F *h_ampVArea_fast  = (TH2F*)histoFileFast->Get("ampvarea");
  
  //area normalize the two runs
  int startBin = h_amp_slow->GetXaxis()->FindBin(0.05);
  double scale = ( h_amp_fast->Integral(startBin, h_amp_fast->GetNbinsX()-1) / h_amp_slow->Integral() );
  cout << scale << endl;
  
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
  
  can->Print("AmpComparison.pdf", "pdf");
  
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
  
  h_area_slow->GetXaxis()->SetRangeUser(0, 20.0);
  h_area_fast->GetXaxis()->SetRangeUser(0, 20.0);  
  h_area_slow->Scale(scale);
  h_area_slow->Draw("hist");
  h_area_fast->Draw("hist same");
  legend->Draw("same");
  can->Print("AreaComparison.pdf", "pdf");
  
  can->Clear();
  can->Update();
  
  h_area_slow->GetXaxis()->SetRangeUser(0, 80.0);
  h_area_fast->GetXaxis()->SetRangeUser(0, 80.0);
  h_area_slow->Scale(scale);
  h_area_slow->Draw("hist");
  h_area_fast->Draw("hist same");
  gPad->SetLogy(1);
  legend->Draw("same");
  can->Print("AreaComparisonLog.pdf", "pdf");
  
  can->Clear();
  can->Update();
  gPad->SetLogy(0);
  
  h_ampOArea_slow->SetLineColor(kBlack);
  h_ampOArea_fast->SetLineColor(kRed+2);
  h_ampOArea_slow->GetXaxis()->SetTitle("ns^{-1]");
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
  can->Print("AmpOverAreaCompare.pdf", "pdf");
  
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
  h_width_slow->GetXaxis()->SetRangeUser(0, 20.0);
  h_width_fast->GetXaxis()->SetRangeUser(0, 20.0);
  h_width_slow->Scale(scale);
  h_width_fast->Draw("hist");
  h_width_slow->Draw("hist same");
  legend->Draw("same");
  can->Print("WidthCompare.pdf", "pdf");
  
  
    




}//end of comparison plots


int LEDHistos::GetOnePE(double percentage){
    //Get the Histogram File
    TFile *histoFile = new TFile("HistogramsLED.root", "READ");
    TH1F *h_area      = (TH1F*)histoFile->Get("area");
    int firstBin = 2; //igonore the zero biin
    int lastBin  = h_area->GetNbinsX();
    double totalArea = h_area->Integral(firstBin, lastBin);
    int rbin = -1;
    
    for(int bin = firstBin; bin < h_area->GetNbinsX(); ++bin){
        double area = h_area->Integral(firstBin, bin+1);
	if( 100*(area / totalArea) >= percentage ){
	   rbin = bin;
	   break;
	}        
    }
    cout << "Return Bin: " << rbin << endl;
    cout << "Retrun Area: " << rbin*h_area->GetBinWidth(rbin) << " V*ns" << endl;
    return rbin;
    
    
    

}

double LEDHistos::CalcLogError(double val, double err){
   const double lnTen = log(10);
   return (fabs(err/val*lnTen));

}

void LEDHistos::SetRedHeatPalette(){
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
