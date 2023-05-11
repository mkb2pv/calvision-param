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

  // convert scintillation energy spectrum to wavelength
  double PhotonWavelength_FAST[nEntries_FAST];
  for(int i=0; i<nEntries_FAST; i++){
      PhotonWavelength_FAST[i] = 1239.8/PhotonEnergy_FAST[i];
  }

  // put the scintillation spectrum into a histogram for use later
  double lambda_step = PhotonWavelength_FAST[nEntries_FAST-2]-PhotonWavelength_FAST[nEntries_FAST-1];
  double lambda_min = PhotonWavelength_FAST[nEntries_FAST-1]-lambda_step/2;
  double lambda_max = PhotonWavelength_FAST[0]+lambda_step/2;

  TH1F *scint_spectrum = new TH1F("scint_spectrum","scintillation spectrum;nm",(int)round((lambda_max-lambda_min)/lambda_step),lambda_min,lambda_max);
  
  for(int i=0; i<nEntries_FAST; i++){
    scint_spectrum->Fill(PhotonWavelength_FAST[i],FastComponent[i]);
  }

  scint_spectrum->Scale(1/scint_spectrum->Integral());

  // fix empty bin in scintillation spectrum by putting in tgraph, interpolating, and restoring to hist
  TGraph *tg_scint_spectrum = new TGraph();
  tg_scint_spectrum->SetTitle("PbWO4 Scintillation Spectrum;Wavelength nm");
  for(int i=0; i<nEntries_FAST; i++){
    tg_scint_spectrum->SetPoint(tg_scint_spectrum->GetN(),PhotonWavelength_FAST[i],FastComponent[i]);
  }

  for(int i=1;i<scint_spectrum->GetNbinsX()+1;i++){
    double x = scint_spectrum->GetBinCenter(i);
    double y = tg_scint_spectrum->Eval(x);
    scint_spectrum->SetBinContent(i,y);
  }

  // store the scintillation yield in a tgraph for later use
  for(int i=0; i<nEntries_SCY; i++){
    ScintilYield[i] = 0.3 * MeV * ScintilYield[i] * ElectronEnergy_SCY[i]; // this formula is used in the original file, unclear why
  }
  TGraph *scint_yield = new TGraph(nEntries_SCY,ElectronEnergy_SCY,ScintilYield);
  scint_yield->SetTitle("Scintillation photon yield vs energy;MeV");
  
  // set parameters for decay time distribution
  p_decay_time->SetParameters(5,15,0.3,0.7);

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
  const int NBINSLAMBDA = 70;
  TCanvas* test = new TCanvas();

  // 2d hist of total number of photons simulated by layer and wavelength
  TH2F *h_tot = new TH2F("h_tot","total events by layer and lambda;pos mm;lamba nm",18,217.5,397.5,NBINSLAMBDA,300,1000);
  tree->Draw((lambd+":inputInitialPosition[2]>>h_tot").c_str(),"","colz");
  gSystem->ProcessEvents();

  // 2d hist of probability of detection of photons by layer and wavelength
  TH2F *h_pdetect = new TH2F("h_pdetect","Probability of Detection vs. Z-position and Wavelength;Position along crystal (z)  mm;Wavelength nm",18,217.5,397.5,NBINSLAMBDA,300,1000);
  tree->Draw((lambd+":inputInitialPosition[2]>>h_pdetect").c_str(),(timeS+">-1").c_str(),"colz");
  h_pdetect->Divide(h_tot);
  h_pdetect->SetStats(0);
  

  TH3F *h_time_pdfs = new TH3F("h_time_pdfs","Time dist by wavelength and z pos;time ns;lambda nm;z pos mm",100,0,10,NBINSLAMBDA,300,1000,18,217.5,397.5);
  tree->Draw(("inputInitialPosition[2]:"+lambd+":"+timeS+">>h_time_pdfs").c_str(),(timeS+">-1").c_str());

  // faster to look up distributions in array than use the projection function in the hit loop
  TH1D *arr_time_pdfs[18][NBINSLAMBDA];
  for(int i=0;i<18;i++){
    for(int j=0;j<NBINSLAMBDA;j++){
      arr_time_pdfs[i][j]=h_time_pdfs->ProjectionX(("h_proj"+to_string(i+1)+"_"+to_string(j+1)).c_str(),j+1,j+1,i+1,i+1);
    }
  }

  h_pdetect->Draw("colz");
  test->Update();

  TCanvas *test2 = new TCanvas();
  tg_scint_spectrum->Draw();
  test2->Update();

  TCanvas tcd = TCanvas();
  p_decay_time->SetTitle("Scintillation Decay Time Distribution;Decay Time ns");
  p_decay_time->Draw();
  tcd.Update();
  
  // setup for reading hits data
  TFile *hit_file = TFile::Open("/project/HEP_EF/calvision/singlebar2/Hayden_results/fastslow/mu_10G_3_withCounts.root");
  TTreeReader reader("tree",hit_file);
  TTreeReaderValue<vector<Float_t>> hit_energy(reader, "ECAL_r_hit_energy");
  TTreeReaderValue<vector<Float_t>> hit_time(reader, "ECAL_r_hit_time");
  TTreeReaderValue<vector<Float_t>> hit_pos(reader, "ECAL_r_hit_zPos");
  //TTreeReaderValue<vector<int>> hit_scin_geant(reader, "ECAL_r_hit_scin"); // geant value of scintillation photons detected due to this hit

  TRandom2 rand = TRandom2();
  
  TH1F *hist_E = new TH1F("hist_E","Hit Energies;Energy Deposited MeV;# Hits",100,0,1.5);
  TH1F *hist_Nphotons = new TH1F("hist_Nphotons","Number of Photons Generated by Hits;# Photons Generated;# Hits",100,0,500);
  TH1F *hist_Nphotons_det = new TH1F("hist_Nphotons_det","Number of Photons Detected by Hits;# Photons Detected;# Hits",6,0,6);
  TH1F *hist_travel = new TH1F("hist_travel","Travel Times of Detected Photons;Travel Time ns;# Photons",100,0,10);
  TH1F *hist_decay = new TH1F("hist_decay","Scintillation Decay Times of Detected Photons;ns",500,0,50);
  TH1F *hist_hit_times = new TH1F("hist_hit_times","Hit Times;ns",100,0.5,1.3);
  TH1F *hist_arrivals = new TH1F("hist_arrivals","Photon Arrival Times (Fast Parameterization);ns",500,0,50);
  TH1F *hist_wavelengths = new TH1F("hist_wavelengths","Generated wavelengths;nm",100,0,1200);
  TH1F *hist_zpos_generated = new TH1F("hist_zpos_generated","Generated zpos;mm",18,217.5,397.5);
  TH1F *hist_zpos_detected = new TH1F("hist_zpos_detected","Detected zpos;mm",18,217.5,397.5);
  TH1F *hist_energy_deposited = new TH1F("hist_energy_deposited","Energy Depositions in Crystal per 10 GeV Muon;mm;MeV",50,217.5,397.5);

  auto t1 = high_resolution_clock::now();

  // loop over all events in the tree
  int nhits, hit_Nphoton;
  double hit_E, hit_z, hit_t, hit_Nphoton_avg, lambda, pdet, travel_t, scint_t;
  int loss = 0;

  while(reader.Next()){
    nhits = (*hit_energy).size();
    
    // loop over all hits in this event
    for(int i =0; i<nhits; i++){
      hit_E = (*hit_energy)[i]*1000; // MeV
      hit_z = (*hit_pos)[i]; // mm
      hit_t = (*hit_time)[i]; // ns

      hist_energy_deposited->Fill(hit_z,hit_E);
      
      //hit_Nphoton_avg = scint_yield->Eval(hit_E *1000); // avg expected scintillation photons from this hit energy, some units issue requires multiplying by 1000 to get right order of magnitude

      hit_Nphoton_avg = 450*hit_E; // 450 photons/MeV is used by Geant

      if(hit_Nphoton_avg<100){
	hit_Nphoton = rand.Poisson(hit_Nphoton_avg);
      }
      else{
	hit_Nphoton = rand.Gaus(hit_Nphoton_avg,TMath::Sqrt(hit_Nphoton_avg));
      }
      
      hist_E->Fill(hit_E);
      hist_hit_times->Fill(hit_t);
      hist_Nphotons->Fill(hit_Nphoton);

      int ndx_pos = (int)floor((hit_z - 217.5)/10);
      int hit_Nphotons_det = 0;
      // treat all photons generated by this hit
      for(int j=0; j<hit_Nphoton; j++){
	lambda = scint_spectrum->GetRandom(); // generate wavelength randomly from scintillation spectrum

	//hist_wavelengths->Fill(lambda);
	hist_zpos_generated->Fill(hit_z);
	int ndx_lam = (int)floor((lambda - 300)/10);
        // probability of detection corresponding to this wavelength and z-position
	pdet = h_pdetect->GetBinContent(ndx_pos+1,ndx_lam+1);


	if(rand.Rndm() > pdet){ // if the photon is not detected continue to the next one
	  loss++;
	  continue;
	}
	
	// generate random travel time using pdfs generated from optical photons
	// ProjectionX() projects the 3d hist onto the time (x) axis
	// the ndx_lam arguments restrict the y-index used for the projection
	// to the corresponding wavelength slice
	// ndx_pos arguments restrict the z-index used to the corresponding
	// z-position slice
	// the +1 is needed because ROOT bin indices start at 1
	//TH1D *proj = h_time_pdfs->ProjectionX("h_proj",ndx_lam+1,ndx_lam+1,ndx_pos+1,ndx_pos+1);
	TH1D *proj = arr_time_pdfs[ndx_pos][ndx_lam]; // faster than ProjectionX()

	// theoretically if the photon has positive probability to be detected
	// then the distribution of travel times should be nonzero
	// but due to binning or some other effect there are a few that 
	// give a 0 integral error, so this avoids that
	if(proj->Integral() == 0){
	  continue;
	}

	travel_t = proj->GetRandom();
	hist_travel->Fill(travel_t);

	// draw a random scintillation decay time
	scint_t = p_decay_time->GetRandom();
	hist_decay->Fill(scint_t);

	// sum up hit time, scintillation decay time, and travel to to 
	// get final time of arrival

	hist_arrivals->Fill(hit_t + scint_t + travel_t);
	hist_wavelengths->Fill(lambda);
	hist_zpos_detected->Fill(hit_z);
	hit_Nphotons_det++;
      }

      hist_Nphotons_det->Fill(hit_Nphotons_det);
    }
  }

  auto t2 = high_resolution_clock::now();
  
  cout << "Total photons detected: " +to_string(hist_arrivals->GetEntries()) << endl;
  cout << "Total photons lost: " + to_string(loss) << endl;

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
  hist_arrivals_geant_cp->SetTitle("Geant;ns");
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

  TCanvas *tc6 = new TCanvas();
  TH1F *hist_zpos_fraction = (TH1F*)hist_zpos_detected->Clone();
  hist_zpos_fraction->Divide(hist_zpos_generated);
  hist_zpos_fraction->SetTitle("fraction of generated photons detected;mm");
  hist_zpos_fraction->Draw();
  tc6->Update();
  
  //TCanvas *tc7 = new TCanvas();
  //tc7->SetLogy();
  //TH1D *proj = h_time_pdfs->ProjectionX("h_proj",0,-1,11,11);
  //proj->SetTitle("Travel Time Probability Distribution for Initial Z-Position 317.5-327.5 mm;Travel time ns");
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


  //TCanvas tc12 = TCanvas();
  //hist_E->Draw("hist");
  //tc12.Update();
  //TCanvas tc13 = TCanvas();
  //hist_Nphotons->Draw("hist");
  //tc13.Update();
  //TCanvas tc14 = TCanvas();
  //hist_Nphotons_det->Draw("hist");
  //tc14.Update();
  //TCanvas tc15 = TCanvas();
  //hist_travel->Draw("hist");
  //tc15.Update();
  

  cout << "\nTo exit, quit ROOT from the File menu of the plot (or use control-C)" << endl;
  theApp.SetIdleTimer(600,".q");  // set up a failsafe timer to end the program
  theApp.Run(true);

  return 0;
}
