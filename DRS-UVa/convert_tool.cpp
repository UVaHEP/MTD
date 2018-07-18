#include <fstream>
#include <iostream>
#include <ios>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string.h>
#include <getopt.h>



#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"

using std::string; 
double findRMS(double ch[1024], float t[1024], float lowtime, float hightime);
double findmean(double ch[1024], float t[1024], float lowtime, float hightime);
double findmean(double ch[1024], float t[1024], float lowtime, float hightime)
{
  double sum = 0;
  int counter = 1;
  for (int i=0; i<1024; i++)
    {
      if (t[i] < lowtime || t[i] > hightime)
	{
	  sum += ch[i];
	  counter++;
	  
	}
      
    }
  double mean = sum/counter;
  return mean;
}

double findRMS(double ch[1024], float t[1024], float lowtime, float hightime)
{
  double sum = 0;
  int counter = 1;
  for (int i=0; i<1024; i++)
    {
      if (t[i] < lowtime || t[i] > hightime)
	{
	  sum += ch[i]*ch[i];
	  counter++;
	  
	}
      
    }
  double RMS = TMath::Sqrt(sum/counter);
  return RMS;
}


bool processFile(std::ifstream &file, TString &oName, int chCount = 1) {

  //Initializing the output ROOT file
  TFile *outfile = new TFile(oName, "recreate");
  //declaring variables to be used in ROOT file
  double chmV[chCount][1024]; //This is the voltage recorded on channel 1 in mV
  float time[chCount][1024];  //This is the time of a sample in ns with time[0] = 0
  int event = 0;
  double RMS[chCount];  //RMS of pedestal
  double amp[chCount];
  double mean[chCount];
  double integral[chCount];
  std::string statements[] = {
    "chmV[%d][1024]/D",
    "time[%d][1024]/F",
    "RMS[1]/D",
    "amp[1]/D",
    "mean[1]/D",
    "integral[%d]/D"
  }; 
  //Initializing the Tree and branches in the ROOT file
  TString Format; 
  TTree *pulse = new TTree("pulse", "This is a tree");
  Format.Form("chmV[%d][1024]/D", chCount); 
  TBranch *ch_b = pulse->Branch("chmV", &chmV, Format);
  Format.Form("time[%d][1024]/F", chCount);  
  TBranch *time_b = pulse->Branch("time", &time, Format); 
  TBranch *event_b = pulse->Branch("event", &event, "event/I");
  Format.Form("RMS[1]/D", chCount);  
  TBranch *RMS_b = pulse->Branch("RMS", &RMS, Format);
  Format.Form("amp[1]/D", chCount);  
  TBranch *amp_b = pulse->Branch("amp", &amp, Format);
  Format.Form("mean[1]/D", chCount);  
  TBranch *mean_b = pulse->Branch("mean", &mean, Format);
  Format.Form("integral[%d]/D", chCount);  
  TBranch *integral_b = pulse->Branch("integral", &integral,Format );


  //Declaring some dummy variables for reading headers we don't care about.
  char tmpHeader[4];
  char tmpsmallHeader[2];
  int SerialNumber;
  short Date[8];
  unsigned short ch[chCount][1024]; //This will be the raw value of the channel, which we will convert to mV later.
  bool endoffile = false;
  char tmpBoardNumber[5];

  //Reading bytes until we get to stuff we actually care about.  
  file.read((char *) &tmpHeader, 4);
  std::cout << "header 1 reads " << tmpHeader << std::endl;
  file.read((char *) &tmpHeader, 4);
  std::cout << "header 2 reads " << tmpHeader << std::endl; 
  file.read((char *) &tmpsmallHeader, 2);
  file.read((char *) &SerialNumber, 2);
  std::cout << "header 3 reads " << tmpsmallHeader << " " << SerialNumber << std::endl;
  //  file.read((char *) &tmpHeader, 4);
  //  std::cout << "header 4 reads " << tmpHeader << std::endl;
  //looping over each channel getting time
  for (int k=0; k<chCount; k++)
    {
      file.read((char *) &tmpHeader, 4);
      file.read((char *) &time[k], 4096);  //This gets the time values in nanoseconds for each sample in the buffer.
      for (int i=1; i<1024; i++)
	{
	  time[k][i] += time[k][i-1];
	}
    }
  std::cout << time[0][0] << "," << time[0][1] << "," << time[0][1020] << std::endl; 

  std::cout << "Time calculated" << std::endl;

  
  while(!file.eof())
    {
      if (event%500==0)
	std::cout << "Processing event #" << event << std::endl;

      //Reading more bytes that we don't care about...
      file.read((char *) &tmpHeader, 4);
      //std::cout << "beginning of loop: " << tmpHeader << std::endl;
      file.read((char *) &SerialNumber, 4);
      file.read((char *) &Date, 16);
      file.read((char *) &tmpBoardNumber, 4);
      file.read((char *) &tmpBoardNumber, 4);
      //file.read((char *) &tmpHeader, 4);
      //file.read((char *) &SerialNumber, 4);
      // std::cout << "Scaler # " << SerialNumber << std::endl;

      //Reading the values of each sample in the buffer for an event
      for (int m=0; m<chCount; m++)
	{
	  file.read((char *) &tmpHeader, 4);
	  file.read((char *) &SerialNumber, 4);
	  file.read((char *) &ch[m], 2048);
	  
      //Converting raw value to mV
	  for (int j=0; j<1024; j++)
	    {
	      chmV[m][j] = (ch[m][j]/65535.-0.5)*1000;

	    }
	  //Need to fix this TGraph to make it an array
	  TGraph *g1 = new TGraph();  //Making a TGraph for each event that we will use to fit each pulse with a function

	  RMS[m] = findRMS(chmV[m], time[m], 100, 140);
	  mean[m] = findmean(chmV[m], time[m], 100, 140);

      //Filling TGraph for each event (pulse)
	  for (int j=0; j<1024; j++)
	    {
	      chmV[m][j] *= -1;
	      //	      chmV[m][j] = mean[m] - chmV[m][j];
	      g1->SetPoint(j, time[m][j], chmV[m][j]);
	  
	    }
     
	  amp[m] = chmV[m][TMath::LocMax(1024, chmV[m])];  //Finding amplitude of the pulse
	  //Should the function be an array????
	  //	  fn1->SetRange(50, 110);  //Setting fit range the match the region in which we see most pulses
      
	  //	  g1->Fit(fn1, "QR"); //Fitting without printing out all of the statistics to the screen (Quiet mode)
      
      //Drawing a pulse to check that the fit is working
	  if (event==100)
	    {
	      TCanvas *c2 =  new TCanvas("c2", "c2", 500, 500);
	      g1->Draw();
	    }
      //Integrating the fit function to get the area under the pulse (related to net charge of the pulse)
	  //TF1 *fn1 = g1->GetFunction("ourfn");
	  //	  integral[m] = fn1->Integral(50, 110);
	}
      event++;
      pulse->Fill(); //Filling the tree's branches
      if (file.eof())
        {
	  std::cout << "End of file reached ......." << std::endl;
	  std::cout << "Processed " << event << " events ..." << std::endl;
          endoffile = true;
          break;
        }

    }
  std::cout << "Tree was filled" << std::endl;
  outfile->Write();
  outfile->Close();
  std::cout << "ROOT file written and closed" << std::endl;


  return true; 
}


int main(int argc, char **argv) {

  int nchan = 1;
  int opt;
  TString input; 
  TString output = "output.root";
  while ((opt = getopt(argc, argv, "f:c:o:")) != -1) {
    switch(opt) {
    case 'f':
      //input file
      input = optarg;
      break;
    case 'c':
      //Channel count
      nchan = atoi(optarg);
      break;
    case 'o':
      //output file
      output = optarg;
      break;
    default:
      std::cerr << "Usage:" << std::endl <<
	"-f <input .dat file>" << std::endl <<
	"-o <output root file>" << std::endl <<
	"-c channelCount" << std::endl;
      break;
    }
  }

  if (input.Length() == 0) {
    std::cout << "Please give me an input file." << std::endl;
    exit(-1); 
  }
  std::ifstream infile; 

  std::cout << "Opening:" << input << std::endl;

  std::cout << "Writing to:" << output << std::endl; 
  
  infile.open(input, std::ios::binary | std::ios::in);
  if (infile.good())  {
    processFile(infile, output);
    infile.close();
  }
  else
    std::cout << "Failed to open:" << input << std::endl; 
  

  return 0; 
}

