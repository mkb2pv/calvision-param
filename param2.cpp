#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <stdio.h>
#include <math.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TSystem.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TRandom2.h>
#include <THStack.h>
#include <TStyle.h>
#include <chrono>
#include "materials.hpp" // contains material information for PWO

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[]){

  auto t0 = high_resolution_clock::now();

  TApplication theApp("App", &argc, argv);
  //gStyle->SetOptStat(0);

  const int NBINSLAMBDA = 70;

  // convert scintillation energy spectrum to wavelength
  double PhotonWavelength_FAST[nEntries_FAST];
  for(int i=0; i<nEntries_FAST; i++){
      PhotonWavelength_FAST[i] = 1239.8/PhotonEnergy_FAST[i];
  }

  // put the scintillation spectrum into a histogram for use later
  TH1F *scint_spectrum = new TH1F("scint_spectrum","scintillation spectrum;nm",NBINSLAMBDA,300,1000);
  
  for(int i=0; i<nEntries_FAST; i++){
    scint_spectrum->Fill(PhotonWavelength_FAST[i],FastComponent[i]);
  }

  //normalize to 1 to use as pdf
  scint_spectrum->Scale(1/scint_spectrum->Integral());
  
  // set parameters for decay time distribution
  p_decay_time->SetParameters(5,15,0.3,0.7);
  TCanvas tcd = TCanvas();


  TChain *chain = new TChain("tree","combined photon files");
    
  for(int i=1;i<=11;i++){
    //chain->Add(("/project/HEP_EF/calvision/singlebar2/singleOP/1Mopticalphoton_"+to_string(i)+".root").c_str());
    //chain->Add(("/project/HEP_EF/calvision/singlebar2/singleOP_new_bh/1Mopticalphoton_"+to_string(i)+".root").c_str());
    //chain->Add(("/project/HEP_EF/calvision/singlebar2/singleOP_new_hh/1Mopticalphoton_"+to_string(i)+".root").c_str());
    chain->Add(("/project/HEP_EF/calvision/singlebar2/Hayden_results/photons/1Mopticalphoton_"+to_string(i)+".root").c_str());
  } 

  TTree *tree=chain;  // rename for convenience

  //string timeS = "SDStime_r_S";
  //string timeC = "SDCtime_r_S";

  //string timeS = "SiPMS_time_r_S";
  //string timeC = "SiPMC_time_r_S";

  string timeS = "SiPMS_time_r_S";

  string lambd = "1239.8/(inputMomentum[3]*1e9)";
  
  // bin all photon events by travel time, wavelength, and z pos
  // note photons that don't reach the crystal have their travel time
  // set to -1 and will therefore be included in the underflow bin
  TH3F *h_time_pdfs = new TH3F("h_time_pdfs","Time dist by wavelength and z pos;time ns;lambda nm;z pos mm",100,0,10,NBINSLAMBDA,300,1000,18,217.5,397.5);
  tree->Draw(("inputInitialPosition[2]:"+lambd+":"+timeS+">>h_time_pdfs").c_str());

  // create time pdfs integrated on wavelength
  TH2F *h_time_pdfs_int = new TH2F("h_time_pdfs_int","Time dist by z pos (integrated over wavelength);time ns;z pos mm",100,0,10,18,217.5,397.5);
  
  for(int i=0;i<=h_time_pdfs_int->GetNbinsX();i++){ // iterate over time
    for(int j=1;j<=h_time_pdfs_int->GetNbinsY();j++){ // iterate over z pos
      double w = 0;
      for(int k=1;k<=NBINSLAMBDA;k++){ // iterate over wavelength
	// sum of number of single photons detected weighted by scintillation spectrum (which is normalized to 1)
	double pspec = scint_spectrum->GetBinContent(k);
	if(pspec <=0) continue;
	w += h_time_pdfs->GetBinContent(i,k,j)*pspec;
      }
      h_time_pdfs_int->SetBinContent(i,j,w);
    }
  }

  // convolve timing distributions with decay time dist
  TH2F *h_time_pdfs_int_conv = new TH2F("h_time_pdfs_int_conv","Time dist by z-pos integrated over wavelength including decay time convolution;time ns;z pos mm",500,0,50,18,217.5,397.5);
  
  for(int i=0;i<=h_time_pdfs_int_conv->GetNbinsX();i++){ // iterate over possible sum travel plus decay times
    for(int j=1;j<=h_time_pdfs_int_conv->GetNbinsY();j++){ // iterate over z pos
      double w = 0;
      for(int k=1;k<=h_time_pdfs_int->GetNbinsX();k++){ // iterate over possible travel times
	// sum over all travel times of number of photons in this layer with that travel time
	// multiplied by probability of decay time equal to travel + decay time minus travel time
	double s = h_time_pdfs_int_conv->GetXaxis()->GetBinCenter(i) - h_time_pdfs_int->GetXaxis()->GetBinCenter(k);
	if(s<0) continue;
	w += h_time_pdfs_int->GetBinContent(k,j)*p_decay_time->Eval(s);
      }
      h_time_pdfs_int_conv->SetBinContent(i,j,w);
    }    
  }

  //normalize all pdfs to the total probability of detection for each layer
  for(int i=1;i<=h_time_pdfs_int_conv->GetNbinsY();i++){ // iterate over position
    TH1D *proj = h_time_pdfs_int->ProjectionX("h_proj",i,i);
    double Ndet = proj->Integral(); // number of single photons in this layer that were detected (after weighting by scintillation spectrum)
    double Nlost = proj->GetBinContent(0); // number of single photons in this layer that were lost (after weighting by scintillation spectrum)
    double norm = Ndet/(Ndet+Nlost)/h_time_pdfs_int_conv->ProjectionX("h_proj2",i,i)->Integral(); // want to normalize to fraction of photons that was detected
    for(int j=1;j<=h_time_pdfs_int_conv->GetNbinsX();j++){
      h_time_pdfs_int_conv->SetBinContent(j,i,h_time_pdfs_int_conv->GetBinContent(j,i)*norm);
    }
  }

  // faster to look up the distributions in an array than use the projection function each time in hit loop
  TH1D *arr_time_pdfs_int_conv[18];
  for(int i=0;i<18;i++){
    arr_time_pdfs_int_conv[i]=h_time_pdfs_int_conv->ProjectionX(("h_proj"+to_string(i+1)).c_str(),i+1,i+1);
  }


  TCanvas *test2 = new TCanvas();
  h_time_pdfs_int->Draw("colz");
  test2->Update();

  TCanvas *test4 = new TCanvas();
  h_time_pdfs_int_conv->Draw("colz");
  test4->Update();

  TCanvas *test5 = new TCanvas();
  scint_spectrum->Draw("hist");
  test5->Update();
  
  // setup for reading hits data
  TFile *hit_file = TFile::Open("/project/HEP_EF/calvision/singlebar2/Hayden_results/fastslow/mu_10G_3_withCounts.root");
  TTreeReader reader("tree",hit_file);
  TTreeReaderValue<vector<Float_t>> hit_energy(reader, "ECAL_r_hit_energy");
  TTreeReaderValue<vector<Float_t>> hit_time(reader, "ECAL_r_hit_time");
  TTreeReaderValue<vector<Float_t>> hit_pos(reader, "ECAL_r_hit_zPos");
  //TTreeReaderValue<vector<int>> hit_scin_geant(reader, "ECAL_r_hit_scin"); // geant value of scintillation photons detected due to this hit

  TRandom2 rand = TRandom2();
  
  TH1F *hist_E = new TH1F("hist_E","Hit Energies;MeV",100,0,1);
  TH1F *hist_Nphotons = new TH1F("hist_Nphotons","Number of Photons Generated by Hits;# photons",100,0,500);
  TH1F *hist_Nphotons_det = new TH1F("hist_Nphotons_det","Number of Photons Detected by Hits;# photons",10,0,10);
  TH1F *hist_travel = new TH1F("hist_travel","Travel Times of Detected Photons;ns",100,0,10);
  TH1F *hist_decay = new TH1F("hist_decay","Scintillation Decay Times of Detected Photons;ns",500,0,50);
  TH1F *hist_hit_times = new TH1F("hist_hit_times","Hit Times;ns",100,0.5,1.3);
  TH1F *hist_arrivals = new TH1F("hist_arrivals","Photon Arrival Times (Fast Parameterization);ns",500,0,50);
  //TH1F *hist_zpos_generated = new TH1F("hist_zpos_generated","Generated zpos;mm",50,217.5,397.5);
  //TH1F *hist_zpos_detected = new TH1F("hist_zpos_detected","Detected zpos;mm",50,217.5,397.5);
  TH1F *hist_energy_deposited = new TH1F("hist_energy_deposited","Energy Depositions in Crystal per 10 GeV Muon;mm;MeV",50,217.5,397.5);

  auto t1 = high_resolution_clock::now();

  // loop over all events in the tree
  int nhits;
  double hit_E, hit_z, hit_t, hit_Nphoton_avg;

  while(reader.Next()){
    nhits = (*hit_energy).size();
    
    // loop over all hits in this event
    for(int i =0; i<nhits; i++){
      hit_E = (*hit_energy)[i]*1000; // MeV
      hit_z = (*hit_pos)[i]; // mm
      hit_t = (*hit_time)[i]; // ns

      hist_energy_deposited->Fill(hit_z,hit_E);

      hit_Nphoton_avg = 450*hit_E; // 450 photons/MeV is used by Geant
      
      int ndx_pos = (int)floor((hit_z - 217.5)/10);
      if(ndx_pos > 17) ndx_pos = 17; // this seems to happen occasionally
      //cout << to_string(ndx_pos) << endl;
      // get photon time distribution for this layer
      //TH1D *proj = h_time_pdfs_int_conv->ProjectionX("h_proj",ndx_pos+1,ndx_pos+1);
      TH1D *proj = arr_time_pdfs_int_conv[ndx_pos]; // faster than ProjectionX()

      // theoretically if the photon has positive probability to be detected
      // then the distribution of travel times should be nonzero
      // but due to binning or some other effect there are a few that 
      // give a 0 integral error, so this avoids that
      if(proj->Integral() == 0){
	continue;
      }

      // shift distribution by hit time
      int shift = (int)round(hit_t/0.1);
      for(int i=hist_arrivals->GetNbinsX();i>=shift+1;i--){
	// add photon time distribution after shifting by hit time to final arrival time distribution
	//hist_arrivals->Fill(hist_arrivals->GetBinCenter(i),proj->GetBinContent(i-shift)*hit_Nphoton_avg);
	hist_arrivals->SetBinContent(i,hist_arrivals->GetBinContent(i)+proj->GetBinContent(i-shift)*hit_Nphoton_avg); // faster than the line above
      }
   
      hist_Nphotons_det->Fill(hit_Nphoton_avg);
      hist_Nphotons->Fill(hit_Nphoton_avg);
    }
  }

  auto t2 = high_resolution_clock::now();

  cout << "Total photons detected: " +to_string(hist_arrivals->Integral()) << endl;

  chrono::duration<double> time1 = t1 - t0;
  cout << "Time for single photon processing: " << time1.count() << " s" << endl;

  chrono::duration<double> time2 = t2 - t1;
  cout << "Time for hit processing: " << time2.count() << " s" << endl;
  cout << "Total time: " << time1.count() + time2.count() << " s" << endl;

  TCanvas *tc = new TCanvas();
  tc->Divide(3,2);
  tc->cd(1);
  hist_E->Draw("hist");
  tc->cd(2);
  hist_Nphotons->Draw("hist");
  tc->cd(3);
  hist_Nphotons_det->Draw("hist");
  tc->cd(4);
  hist_travel->Draw("hist");
  tc->cd(5);
  hist_hit_times->Draw("hist");
  tc->cd(6);
  hist_decay->Draw("hist");
  tc->Update();


  // now compare the same histogram of arrival times from Geant simulation
  TH1F *hist_arrivals_geant = (TH1F*)hit_file->Get("h_phot_time_SiPMS_Scin");
  hist_arrivals_geant->SetTitle("Photon arrival times (Geant)");

  TCanvas *tc3 = new TCanvas();
  tc3->Divide(2,1);
  tc3->cd(1);
  hist_arrivals->Draw("hist");
  tc3->cd(2);
  hist_arrivals_geant->Draw("hist");
  tc3->Update();

  TCanvas *tc4 = new TCanvas();
  THStack *stack = new THStack("stack","Comparison of Geant and Fast Photon Arrival Times;Arrival Time ns;Number of Photons");
  TH1F *hist_arrivals_geant_cp = (TH1F*)hist_arrivals_geant->Clone();
  hist_arrivals_geant_cp->SetTitle("Geant;");
  //hist_arrivals_geant_cp->Scale(1/hist_arrivals_geant_cp->Integral());
  TH1F *hist_arrivals_cp = (TH1F*)hist_arrivals->Clone();
  hist_arrivals_cp->SetTitle("Fast;ns");
  hist_arrivals_cp->SetLineColor(kRed);
  //hist_arrivals_cp->Scale(1/hist_arrivals_cp->Integral());
  stack->Add(hist_arrivals_cp,"hist");
  stack->Add(hist_arrivals_geant_cp,"hist");
  stack->Draw("nostack");
  gPad->BuildLegend(0.75,0.75,0.95,0.95);
  tc4->Update();

  // additional diagnostic plots

  //TCanvas *tc5 = new TCanvas();
  //tc5->Divide(2,1);
  //tc5->cd(1);
  //hist_wavelengths->Scale(1/hist_wavelengths->Integral());
  //hist_wavelengths->Draw("hist");
  //hist_zpos_generated->Draw("hist");
  //tc5->cd(2);
  //hist_zpos_detected->Draw("hist");
  // scint_spectrum->Draw("hist");
  //tc5->Update();

  //TCanvas *tc6 = new TCanvas();
  //TH1F *hist_zpos_fraction = (TH1F*)hist_zpos_detected->Clone();
  //hist_zpos_fraction->Divide(hist_zpos_generated);
  //hist_zpos_fraction->SetTitle("fraction of generated photons detected;mm");
  //hist_zpos_fraction->Draw();
  //tc6->Update();
  
  //TCanvas *tc7 = new TCanvas();
  //tc7->SetLogy();
  //TH1D *proj = h_time_pdfs->ProjectionX("h_proj",0,-1,11,11);
  //proj->SetTitle("Travel Time Probability Distribution for Initial Z-Position 317.5-327.5 mm;ns");
  //proj->SetStats(0);
  //proj->Scale(1/proj->Integral());
  //proj->Draw("hist");
  //h_time_pdfs->Draw("box");
  //tc7->Update();

  //TCanvas *tc8 = new TCanvas();
  // hist_energy_deposited->Scale(1/500.0);
  //hist_energy_deposited->Draw("hist min0");
  //hist_energy_deposited->SetStats(0);
  //tc8->Update();

  ///TCanvas *tc9 = new TCanvas();
  //tc9->Divide(2,1);
  //tc9->cd(1);
  //hist_E->Draw("hist");
  //tc9->cd(2);
  //hist_Nphotons->Draw("hist");
  //tc9->Update();

  //TCanvas *tc10 = new TCanvas();
  //tc10->Divide(2,1);
  //tc10->cd(1);
  //hist_Nphotons_det->Draw("hist");
  //tc10->cd(2);
  //hist_travel->Draw("hist");
  //tc10->Update();

  //TCanvas *tc11 = new TCanvas();
  //tc11->Divide(2,1);
  //tc11->cd(1);
  //hist_hit_times->Draw("hist");
  //tc11->cd(2);
  //hist_decay->Draw("hist");
  //tc11->Update();

  cout << "\nTo exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;
  theApp.SetIdleTimer(600,".q");  // set up a failsafe timer to end the program
  theApp.Run(true);

  return 0;
}
