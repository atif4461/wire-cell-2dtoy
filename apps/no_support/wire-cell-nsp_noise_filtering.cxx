#include "WCPSst/GeomDataSource.h"
#include "WCPSst/DatauBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/BadTiling.h"

#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/SimpleBlobToyTiling.h"

#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMatrixIterate_SingleWire.h"
#include "WCP2dToy/ToyMatrixIterate_Only.h"


#include "WCP2dToy/ToyMatrixMarkov.h"
#include "WCP2dToy/ToyMetric.h"
#include "WCP2dToy/BlobMetric.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

#include "WCPData/GeomCluster.h"


#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"
#include "WCP2dToy/ToySignalSimu.h"
#include "WCP2dToy/ToySignalSimuTrue.h"

#include "WCP2dToy/DataSignalGaus.h"
#include "WCP2dToy/DataSignalWien_ROI.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include <iostream>

using namespace WCP;
using namespace std;




/**
 * @brief Main program entry point.
 *
 * This function is the main entry point of the program. It reads command line arguments,
 * initializes data sources, and performs various operations such as reading frames,
 * filling histograms, and writing output files.
 *
 * @param argc Number of command line arguments.
 * @param argv Array of command line argument strings.
 * @return Program exit status.
 
int main(int argc, char* argv[])
 * @brief Prints usage message and exits if insufficient command line arguments are provided.
 
if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num " << endl;
    return 1;
}
 * @brief Initializes geometry data source from file specified by command line argument.
 
WCPSst::GeomDataSource gds(argv[1]);
 * @brief Retrieves extent of geometry data source.
 
std::vector<double> ex = gds.extent();
 * @brief Prints extent of geometry data source to error stream.
 
cerr << "Extent: "
     << " x:" << ex[0]/units::mm << " mm"
     << " y:" << ex[1]/units::m << " m"
     << " z:" << ex[2]/units::m << " m"
     << endl;
 * @brief Prints pitch values for each wire plane type.
 
cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
     << " " << gds.pitch(WirePlaneType_t(1)) 
     << " " << gds.pitch(WirePlaneType_t(2))
     << endl;
 * @brief Prints angle values for each wire plane type.
 
cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
     << " " << gds.angle(WirePlaneType_t(1)) 
     << " " << gds.angle(WirePlaneType_t(2))
     << endl;
 * @brief Defines unit distance value.
 
float unit_dis = 1.119;  // 70 KV @ 273 V/cm
 * @brief Defines total time bin count.
 
int total_time_bin=9594;
 * @brief Defines frame length.
 
int frame_length = 3200;
 * @brief Defines maximum number of events.
 
int max_events = 100;
 * @brief Parses event number from command line argument.
 
int eve_num  = atoi(argv[3]);
 * @brief Defines time offset value.
 
int time_offset = -92.;
 * @brief Retrieves root file path from command line argument.
 
const char* root_file = argv[2];
 * @brief Declares variables to store run, subrun, and event numbers.
 
int run_no, subrun_no, event_no;
 * @brief Initializes frame data source from root file and geometry data source.
 
WCPSst::DatauBooNEFrameDataSource data_fds(root_file,gds,total_time_bin);
 * @brief Jumps to specified event in frame data source.
 
data_fds.jump(eve_num);
 * @brief Retrieves run, subrun, and event numbers from frame data source.
 
run_no = data_fds.get_run_no();
subrun_no = data_fds.get_subrun_no();
event_no = data_fds.get_event_no();
 * @brief Prints run, subrun, and event numbers to console.
 
cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;
 * @brief Retrieves u, v, and w plane maps from frame data source.
 
ChirpMap& uplane_map = data_fds.get_u_cmap();
ChirpMap& vplane_map = data_fds.get_v_cmap();
ChirpMap& wplane_map = data_fds.get_w_cmap();
 * @brief Opens output root file for writing.
 
TFile *file = new TFile(Form("nsp1_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
 * @brief Selects wires in each plane type.
 
GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));
 * @brief Retrieves sizes of wire selections.
 
Int_t nwire_u = wires_u.size();
Int_t nwire_v = wires_v.size();
Int_t nwire_w = wires_w.size();
 * @brief Creates raw histogram objects for each plane type.
 
TH2F *hu_raw = new TH2F("hu_raw","hu_raw",nwire_u,-0.5,nwire_u-0.5,total_time_bin,0,total_time_bin);
TH2F *hv_raw = new TH2F("hv_raw","hv_raw",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin,0,total_time_bin);
TH2F *hw_raw = new TH2F("hw_raw","hw_raw",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin,0,total_time_bin);
 * @brief Temporary histogram pointer.
 
TH2F *htemp;
 * @brief Retrieves current frame from frame data source.
 
const Frame& frame = data_fds.get();
 * @brief Iterates over traces in frame and fills histograms accordingly.
 
size_t ntraces = frame.traces.size();
for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp = hu_raw;
    }else if (plane == WirePlaneType_t(1)){
      htemp = hv_raw;
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp = hw_raw;
      chid -= nwire_u + nwire_v;
    }
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(chid+1,tt,trace.charge.at(i));
    }
}
 * @brief Creates tree object for storing run information.
 
TTree *Trun = new TTree("Trun","Trun");
 * @brief Sets directory for tree object.
 
Trun->SetDirectory(file);
 * @brief Branches for storing detector, event, run, and subrun numbers.
 
int detector = 0; // MicroBooNE
Trun->Branch("detector",&detector,"detector/I");
Trun->Branch("eventNo",&event_no,"eventNo/I");
Trun->Branch("runNo",&run_no,"runNo/I");
Trun->Branch("subRunNo",&subrun_no,"subRunNo/I");
 * @brief Branches for storing unit distance, total time bin, frame length, and maximum events.
 
Trun->Branch("unit_dis",&unit_dis,"unit_dis/F");
Trun->Branch("total_time_bin",&total_time_bin,"total_time_bin/I");
Trun->Branch("frame_length",&frame_length,"frame_length/I");
Trun->Branch("max_events",&max_events,"max_events/I");
Trun->Branch("eve_num",&eve_num,"eve_num/I");
Trun->Branch("time_offset",&time_offset,"time_offset/I");
 * @brief Fills tree object with current values.
 
Trun->Fill();
 * @brief Creates tree object for storing noisy channel information.
 
TTree *T_lf = new TTree("T_lf","T_lf");
 * @brief Sets directory for tree object.
 
T_lf->SetDirectory(file);
 * @brief Branch for storing channel number.
 
Int_t channel;
T_lf->Branch("channel",&channel,"channel/I");
 * @brief Retrieves set of noisy channels from frame data source.
 
std::set<int> lf_noisy_channels = data_fds.get_lf_noisy_channels();
 * @brief Iterates over noisy channels and fills tree object.
 
for (auto it = lf_noisy_channels.begin(); it!= lf_noisy_channels.end(); it++){
    channel = *it;
    T_lf->Fill();
}
 * @brief Creates tree object for storing bad channel information.
 
TTree *T_bad = new TTree("T_bad","T_bad");
 * @brief Sets directory for tree object.
 
T_bad->SetDirectory(file);
 * @brief Branches for storing channel ID, plane, start time, and end time.
 
Int_t chid, plane;
Int_t start_time,end_time;
T_bad->Branch("chid",&chid,"chid/I");
T_bad->Branch("plane",&plane,"plane/I");
T_bad->Branch("start_time",&start_time,"start_time/I");
T_bad->Branch("end_time",&end_time,"end_time/I");
 * @brief Iterates over u, v, and w plane maps and fills tree object.
 
for (auto it = uplane_map.begin(); it!=uplane_map.end();it++){
    chid = it->first;
    plane = 0;
    start_time = it->second.first;
    end_time = it->second.second;
    T_bad->Fill();
}
for (auto it = vplane_map.begin(); it!=vplane_map.end();it++){
    chid = it->first + nwire_u;
    plane = 1;
    start_time = it->second.first;
    end_time = it->second.second;
    T_bad->Fill();
}
for (auto it = wplane_map.begin(); it!=wplane_map.end();it++){
    chid = it->first + nwire_u + nwire_v;
    plane = 2;
    start_time = it->second.first;
    end_time = it->second.second;
    T_bad->Fill();
}
 * @brief Writes output file to disk.
 
file->Write();
 * @brief Closes output file.
 
file->Close();* This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.2.
*/ 
int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root eve_num " << endl;
    return 1;
  }

  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;
  cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
       << " " << gds.pitch(WirePlaneType_t(1)) 
       << " " << gds.pitch(WirePlaneType_t(2))
       << endl;
  cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
       << " " << gds.angle(WirePlaneType_t(1)) 
       << " " << gds.angle(WirePlaneType_t(2))
       << endl;
  
  float unit_dis = 1.119;  // 70 KV @ 273 V/cm

  int total_time_bin=9594;

  int frame_length = 3200;
  int max_events = 100;
  int eve_num  = atoi(argv[3]);

  int time_offset = -92.;
  
  const char* root_file = argv[2];
   
  int run_no, subrun_no, event_no;
  
  WCPSst::DatauBooNEFrameDataSource data_fds(root_file,gds,total_time_bin);
  data_fds.jump(eve_num);

  run_no = data_fds.get_run_no();
  subrun_no = data_fds.get_subrun_no();
  event_no = data_fds.get_event_no();
  
  

  cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;
  
  ChirpMap& uplane_map = data_fds.get_u_cmap();
  ChirpMap& vplane_map = data_fds.get_v_cmap();
  ChirpMap& wplane_map = data_fds.get_w_cmap();
  
  TFile *file = new TFile(Form("nsp1_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  Int_t nwire_u = wires_u.size();
  Int_t nwire_v = wires_v.size();
  Int_t nwire_w = wires_w.size();

  TH2F *hu_raw = new TH2F("hu_raw","hu_raw",nwire_u,-0.5,nwire_u-0.5,total_time_bin,0,total_time_bin);
  TH2F *hv_raw = new TH2F("hv_raw","hv_raw",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin,0,total_time_bin);
  TH2F *hw_raw = new TH2F("hw_raw","hw_raw",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin,0,total_time_bin);


  TH2F *htemp;
  
  const Frame& frame = data_fds.get();
  size_t ntraces = frame.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp = hu_raw;
    }else if (plane == WirePlaneType_t(1)){
      htemp = hv_raw;
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp = hw_raw;
      chid -= nwire_u + nwire_v;
    }
    for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp->SetBinContent(chid+1,tt,trace.charge.at(i));
    }
  }
  
  // finish saving

  TTree *Trun = new TTree("Trun","Trun");
  Trun->SetDirectory(file);

  int detector = 0; // MicroBooNE
  Trun->Branch("detector",&detector,"detector/I");

  Trun->Branch("eventNo",&event_no,"eventNo/I");
  Trun->Branch("runNo",&run_no,"runNo/I");
  Trun->Branch("subRunNo",&subrun_no,"subRunNo/I");  
  Trun->Branch("unit_dis",&unit_dis,"unit_dis/F");
  Trun->Branch("total_time_bin",&total_time_bin,"total_time_bin/I");
  Trun->Branch("frame_length",&frame_length,"frame_length/I");
  Trun->Branch("max_events",&max_events,"max_events/I");
  Trun->Branch("eve_num",&eve_num,"eve_num/I");
  Trun->Branch("time_offset",&time_offset,"time_offset/I");
  Trun->Fill();

  TTree *T_lf = new TTree("T_lf","T_lf");
  Int_t channel;
  T_lf->SetDirectory(file);
  T_lf->Branch("channel",&channel,"channel/I");
  std::set<int> lf_noisy_channels = data_fds.get_lf_noisy_channels();
  for (auto it = lf_noisy_channels.begin(); it!= lf_noisy_channels.end(); it++){
    channel = *it;
    T_lf->Fill();
  }

  TTree *T_bad = new TTree("T_bad","T_bad");
  Int_t chid, plane;
  Int_t start_time,end_time;
  T_bad->Branch("chid",&chid,"chid/I");
  T_bad->Branch("plane",&plane,"plane/I");
  T_bad->Branch("start_time",&start_time,"start_time/I");
  T_bad->Branch("end_time",&end_time,"end_time/I");
  T_bad->SetDirectory(file);
  for (auto it = uplane_map.begin(); it!=uplane_map.end();it++){
    chid = it->first;
    plane = 0;
    start_time = it->second.first;
    end_time = it->second.second;
    T_bad->Fill();
  }
  for (auto it = vplane_map.begin(); it!=vplane_map.end();it++){
    chid = it->first + nwire_u;
    plane = 1;
    start_time = it->second.first;
    end_time = it->second.second;
    T_bad->Fill();
  }
  for (auto it = wplane_map.begin(); it!=wplane_map.end();it++){
    chid = it->first + nwire_u + nwire_v;
    plane = 2;
    start_time = it->second.first;
    end_time = it->second.second;
    T_bad->Fill();
  }
  



  file->Write();
  file->Close();
  
  

}
