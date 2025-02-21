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

//#include "WCP2dToy/DataSignalGaus.h"
//#include "WCP2dToy/DataSignalWien_ROI.h"

#include "WCP2dToy/uBooNE_Data_2D_Deconvolution.h"

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
 * @brief Main program entry point
 * 
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 * @return Program exit status
 
int main(int argc, char* argv[]) 
 * @brief Prints usage message and exits if insufficient arguments are provided
 
if (argc < 3) {
  
}
 * @brief Creates a geometry data source object from a file
 * 
 * @param argv[1] Path to ChannelWireGeometry.txt file
 
WCPSst::GeomDataSource gds(argv[1]);
 * @brief Retrieves and prints the extent of the geometry data source
 
std::vector<double> ex = gds.extent();
 * @brief Prints pitch values for each wire plane type
 
cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
       << " " << gds.pitch(WirePlaneType_t(1)) 
       << " " << gds.pitch(WirePlaneType_t(2))
       << endl;
 * @brief Prints angle values for each wire plane type
 
cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
       << " " << gds.angle(WirePlaneType_t(1)) 
       << " " << gds.angle(WirePlaneType_t(2))
       << endl;
 * @brief Opens a ROOT file and retrieves trees
 * 
 * @param argv[2] Path to root file
 
TString filename = argv[2];
TFile *file = new TFile(filename);
TTree *Trun = (TTree*)file->Get("Trun");
TTree *T_bad = (TTree*)file->Get("T_bad");
TTree *T_lf = (TTree*)file->Get("T_lf");
 * @brief Clones trees
 
TTree *Trun1 = Trun->CloneTree();
TTree *T_bad1 = T_bad->CloneTree();
TTree *T_lf1 = T_lf->CloneTree();
 * @brief Retrieves histograms from the file
 
TH2F *hu_raw = (TH2F*)file->Get("hu_raw");
TH2F *hv_raw = (TH2F*)file->Get("hv_raw");
TH2F *hw_raw = (TH2F*)file->Get("hw_raw");
 * @brief Defines time offsets
 
float toffset_1=0.0;
float toffset_2=0.0;
float toffset_3=0.0;
 * @brief Creates a uBooNE frame data source object
 * 
 * @param hu_raw Horizontal raw histogram
 * @param hv_raw Vertical raw histogram
 * @param hw_raw Depth raw histogram
 * @param T_bad Bad tree
 * @param T_lf LF tree
 * @param Trun Run tree
 * @param gds Geometry data source
 
WCPSst::DatauBooNEFrameDataSource data_fds(hu_raw,hv_raw,hw_raw,T_bad,T_lf,Trun,gds);
 * @brief Retrieves run, subrun, and event numbers
 
int run_no = data_fds.get_run_no();
int subrun_no = data_fds.get_subrun_no();
int event_no = data_fds.get_event_no();
 * @brief Retrieves channel maps
 
ChirpMap& uplane_map = data_fds.get_u_cmap();
ChirpMap& vplane_map = data_fds.get_v_cmap();
ChirpMap& wplane_map = data_fds.get_w_cmap();
 * @brief Performs deconvolution using Wiener filter
 * 
 * @param data_fds Data source
 * @param gds Geometry data source
 * @param uplane_map U-plane channel map
 * @param vplane_map V-plane channel map
 * @param wplane_map W-plane channel map
 * @param bins_per_frame Bins per frame
 * @param toffset_1 Time offset 1
 * @param toffset_2 Time offset 2
 * @param toffset_3 Time offset 3
 
WCP2dToy::uBooNEData2DDeconvolutionFDS wien_fds(data_fds,gds,uplane_map, vplane_map, wplane_map,100,toffset_1,toffset_2,toffset_3);
 * @brief Jumps to the first event
 
wien_fds.jump(0);
 * @brief Retrieves wire selections for each plane
 
GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));
 * @brief Counts number of wires in each plane
 
Int_t nwire_u = wires_u.size();
Int_t nwire_v = wires_v.size();
Int_t nwire_w = wires_w.size();
 * @brief Creates a new ROOT file and sets its directory
 * 
 * @param run_no Run number
 * @param subrun_no Subrun number
 * @param event_no Event number
 
TFile *file1 = new TFile(Form("nsp2_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
Trun1->SetDirectory(file1);
T_bad1->SetDirectory(file1);
T_lf1->SetDirectory(file1);
 * @brief Writes trees to the new file
 
Trun1->Write();
T_bad1->Write();
T_lf1->Write();
 * @brief Retrieves the number of time bins per frame
 
int total_time_bin = wien_fds.Get_Bins_Per_Frame();
 * @brief Creates histograms for deconvolved data
 * 
 * @param nwire_u Number of wires in U-plane
 * @param nwire_v Number of wires in V-plane
 * @param nwire_w Number of wires in W-plane
 * @param total_time_bin Total time bins per frame
 
TH2I *hu_decon = new TH2I("hu_decon","hu_decon",nwire_u,-0.5,nwire_u-0.5,total_time_bin,0,total_time_bin);
TH2I *hv_decon = new TH2I("hv_decon","hv_decon",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin,0,total_time_bin);
TH2I *hw_decon = new TH2I("hw_decon","hw_decon",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin,0,total_time_bin);
 * @brief Fills histograms with deconvolved data
 * 
 * @param frame1 Deconvolved frame
 * @param ntraces Number of traces in the frame
 
const Frame& frame1 = wien_fds.get();
int ntraces = frame1.traces.size();
for (size_t ind=0; ind<ntraces; ++ind) {
  const Trace& trace = frame1.traces[ind];
  int tbin = trace.tbin;
  int chid = trace.chid;
  int nbins = trace.charge.size();
  WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
  if (plane == WirePlaneType_t(0)){
    htemp1 = hu_decon;
  }else if (plane == WirePlaneType_t(1)){
    htemp1 = hv_decon;
    chid -= nwire_u;
  }else if (plane == WirePlaneType_t(2)){
    htemp1 = hw_decon;
    chid -= nwire_u + nwire_v;
  }
   for (int i = tbin;i!=tbin+nbins;i++){
    int tt = i+1;
    htemp1->SetBinContent(chid+1,tt,trace.charge.at(i));
  }
}
 * @brief Writes the new file and closes it
 
file1->Write();
file1->Close();* This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.2.
*/ 
int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root " << endl;
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

  TString filename = argv[2];
  TFile *file = new TFile(filename);
  TTree *Trun = (TTree*)file->Get("Trun");
  TTree *T_bad = (TTree*)file->Get("T_bad");
  TTree *T_lf = (TTree*)file->Get("T_lf");

  TTree *Trun1 = Trun->CloneTree();
  TTree *T_bad1 = T_bad->CloneTree();
  TTree *T_lf1 = T_lf->CloneTree();

  TH2F *hu_raw = (TH2F*)file->Get("hu_raw");
  TH2F *hv_raw = (TH2F*)file->Get("hv_raw");
  TH2F *hw_raw = (TH2F*)file->Get("hw_raw");

  

 
  
  float toffset_1=0.0;//-0.787;
  float toffset_2=0.0;//-0.603;
  float toffset_3=0.0;

  // toffset_1 = 0.;
  // toffset_2 = 0.;
  
  
  WCPSst::DatauBooNEFrameDataSource data_fds(hu_raw,hv_raw,hw_raw,T_bad,T_lf,Trun,gds);


  //data_fds.jump(0);

  int run_no = data_fds.get_run_no();
  int subrun_no = data_fds.get_subrun_no();
  int event_no = data_fds.get_event_no();
  cout << "Run No: " << run_no << " " << subrun_no << " " << event_no << endl;
  
 
  


  ChirpMap& uplane_map = data_fds.get_u_cmap();
  ChirpMap& vplane_map = data_fds.get_v_cmap();
  ChirpMap& wplane_map = data_fds.get_w_cmap();
  cout << "Bad Channels: " << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() << endl;

  // const Frame& frame1 = data_fds.get();
  // size_t ntraces = frame1.traces.size();
  // cout << ntraces << endl;

  
  cout << "Deconvolution with Wiener filter" << endl; 
  WCP2dToy::uBooNEData2DDeconvolutionFDS wien_fds(data_fds,gds,uplane_map, vplane_map, wplane_map,100,toffset_1,toffset_2,toffset_3);
  wien_fds.jump(0);
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  Int_t nwire_u = wires_u.size();
  Int_t nwire_v = wires_v.size();
  Int_t nwire_w = wires_w.size();


  
  
  
  TFile *file1 = new TFile(Form("nsp2_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
  Trun1->SetDirectory(file1);
  T_bad1->SetDirectory(file1);
  T_lf1->SetDirectory(file1);

  Trun1->Write();
  T_bad1->Write();
  T_lf1->Write();

  int total_time_bin = wien_fds.Get_Bins_Per_Frame();

  
  TH2I *hu_decon = new TH2I("hu_decon","hu_decon",nwire_u,-0.5,nwire_u-0.5,total_time_bin,0,total_time_bin);
  TH2I *hv_decon = new TH2I("hv_decon","hv_decon",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin,0,total_time_bin);
  TH2I *hw_decon = new TH2I("hw_decon","hw_decon",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin,0,total_time_bin);

  //save histograms ... 
  TH2I *htemp1;
  const Frame& frame1 = wien_fds.get();
  int ntraces = frame1.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp1 = hu_decon;
    }else if (plane == WirePlaneType_t(1)){
      htemp1 = hv_decon;
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp1 = hw_decon;
      chid -= nwire_u + nwire_v;
    }
     for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp1->SetBinContent(chid+1,tt,trace.charge.at(i));
    }
  }

    
  file1->Write();
  file1->Close();


}
