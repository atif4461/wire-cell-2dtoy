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
#include "WCP2dToy/uBooNE_Data_ROI.h"
#include "WCP2dToy/uBooNE_Data_After_ROI.h"

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
 * @brief Check command line argument count
 *
 * @return Error message if insufficient arguments
 
if (argc < 4) {
     * @brief Print usage message
 
cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt decon.root raw.root" << endl;
    return 1;
}
 * @brief Create geometry data source object
 *
 * @param argv Command line argument array
 
WCPSst::GeomDataSource gds(argv[1]);
 * @brief Get extent of geometry data source
 *
 * @return Vector of doubles representing extent
 
std::vector<double> ex = gds.extent();
 * @brief Print extent values
 
cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;
 * @brief Print pitch values for each wire plane type
 
cout << "Pitch: " << gds.pitch(WirePlaneType_t(0)) 
       << " " << gds.pitch(WirePlaneType_t(1)) 
       << " " << gds.pitch(WirePlaneType_t(2))
       << endl;
 * @brief Print angle values for each wire plane type
 
cout << "Angle: " << gds.angle(WirePlaneType_t(0)) 
       << " " << gds.angle(WirePlaneType_t(1)) 
       << " " << gds.angle(WirePlaneType_t(2))
       << endl;
 * @brief Open root file
 *
 * @param argv Command line argument array
 
TString filename = argv[2];
TFile *file = new TFile(filename);
 * @brief Retrieve trees from root file
 
TTree *Trun = (TTree*)file->Get("Trun");
TTree *T_bad = (TTree*)file->Get("T_bad");
TTree *T_lf = (TTree*)file->Get("T_lf");
 * @brief Retrieve histograms from root file
 
TH2I *hu_decon = (TH2I*)file->Get("hu_decon");
TH2I *hv_decon = (TH2I*)file->Get("hv_decon");
TH2I *hw_decon = (TH2I*)file->Get("hw_decon");
 * @brief Open second root file
 *
 * @param argv Command line argument array
 
filename = argv[3];
TFile *file2 = new TFile(filename);
 * @brief Retrieve histograms from second root file
 
TH2F *hu_raw = (TH2F*)file2->Get("hu_raw");
TH2F *hv_raw = (TH2F*)file2->Get("hv_raw");
TH2F *hw_raw = (TH2F*)file2->Get("hw_raw");
 * @brief Create data frame data source object
 *
 * @param hu_raw Horizontal raw histogram
 * @param hv_raw Vertical raw histogram
 * @param hw_raw Depth raw histogram
 * @param T_bad Bad tree
 * @param T_lf LF tree
 * @param Trung Run tree
 * @param gds Geometry data source
 
WCPSst::DatauBooNEFrameDataSource data_fds(hu_raw,hv_raw,hw_raw,T_bad,T_lf,Trun,gds);
 * @brief Create uBooNE data 2D deconvolution FDs object
 *
 * @param hu_decon Horizontal deconvoluted histogram
 * @param hv_decon Vertical deconvoluted histogram
 * @param hw_decon Depth deconvoluted histogram
 * @param T_bad Bad tree
 * @param gds Geometry data source
 
WCP2dToy::uBooNEData2DDeconvolutionFDS wien_fds(hu_decon,hv_decon,hw_decon,T_bad, gds);
 * @brief Retrieve chirp maps
 
ChirpMap& uplane_map = wien_fds.get_u_cmap();
ChirpMap& vplane_map = wien_fds.get_v_cmap();
ChirpMap& wplane_map = wien_fds.get_w_cmap();
 * @brief Retrieve set of noisy channels
 
std::set<int>& lf_noisy_channels = data_fds.get_lf_noisy_channels();
 * @brief Declare variables for run number, sub-run number, and event number
 
int run_no;
int subrun_no;
int event_no;
 * @brief Set branch addresses for run tree
 
Trun->SetBranchAddress("eventNo",&event_no);
Trun->SetBranchAddress("runNo",&run_no);
Trun->SetBranchAddress("subRunNo",&subrun_no);
 * @brief Get first entry from run tree
 
Trun->GetEntry(0);
 * @brief Declare rebining factor
 
int rebin = 4;
 * @brief Create uBooNE ROI object
 *
 * @param data_fds Data frame data source
 * @param wien_fds uBooNE data 2D deconvolution FDS
 * @param gds Geometry data source
 * @param uplane_map U-plane chirp map
 * @param vplane_map V-plane chirp map
 * @param wplane_map W-plane chirp map
 * @param lf_noisy_channels Noisy channels
 
WCP2dToy::uBooNEDataROI uboone_rois(data_fds,wien_fds,gds,uplane_map,vplane_map,wplane_map,lf_noisy_channels);
 * @brief Create ROI FDS object
 *
 * @param wien_fds uBooNE data 2D deconvolution FDS
 * @param gds Geometry data source
 * @param uboone_rois uBooNE ROI
 * @param rebin Rebinning factor
 
WCP2dToy::uBooNEDataAfterROI roi_fds(wien_fds,gds,uboone_rois,rebin);
 * @brief Jump to first entry in ROI FDS
 
roi_fds.jump(0);
 * @brief Retrieve wire selections for each plane type
 
GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));
 * @brief Get sizes of wire selections
 
Int_t nwire_u = wires_u.size();
Int_t nwire_v = wires_v.size();
Int_t nwire_w = wires_w.size();
 * @brief Get total time bins
 
int  total_time_bin = hu_decon->GetNbinsY();
 * @brief Open output root file
 *
 * @param run_no Run number
 * @param subrun_no Sub-run number
 * @param event_no Event number
 
TFile *file1 = new TFile(Form("nsp3_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");
 * @brief Create histograms for ROI
 
TH2F *hu_roi = new TH2F("hu_roi","hu_roi",nwire_u,-0.5,nwire_u-0.5,total_time_bin/rebin,0,total_time_bin);
TH2F *hv_roi = new TH2F("hv_roi","hv_roi",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin/rebin,0,total_time_bin);
TH2F *hw_roi = new TH2F("hw_roi","hw_roi",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin/rebin,0,total_time_bin);
 * @brief Get current frame from ROI FDS
 
const Frame& frame1 = roi_fds.get();
 * @brief Loop through traces in frame
 
int ntraces = frame1.traces.size();
for (size_t ind=0; ind<ntraces; ++ind) {
         * @brief Get current trace
     
    const Trace& trace = frame1.traces[ind];

         * @brief Get time bin, channel ID, and charge size
     
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();

         * @brief Determine plane type and adjust channel ID accordingly
     
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    TH2F* htemp1;
    if (plane == WirePlaneType_t(0)){
        htemp1 = hu_roi;
    }else if (plane == WirePlaneType_t(1)){
        htemp1 = hv_roi;
        chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
        htemp1 = hw_roi;
        chid -= nwire_u + nwire_v;
    }

         * @brief Fill histogram with charge values
     
    for (int i = tbin;i!=tbin+nbins;i++){
        int tt = i+1;
        htemp1->SetBinContent(chid+1,tt,trace.charge.at(i));
    }
}
 * @brief Write and close output root file
 
file1->Write();
file1->Close();* This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.01.
*/ 
int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt decon.root raw.root" << endl;
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

  TH2I *hu_decon = (TH2I*)file->Get("hu_decon");
  TH2I *hv_decon = (TH2I*)file->Get("hv_decon");
  TH2I *hw_decon = (TH2I*)file->Get("hw_decon");
  
  filename = argv[3];
  TFile *file2 = new TFile(filename);
  TH2F *hu_raw = (TH2F*)file2->Get("hu_raw");
  TH2F *hv_raw = (TH2F*)file2->Get("hv_raw");
  TH2F *hw_raw = (TH2F*)file2->Get("hw_raw");
 
  WCPSst::DatauBooNEFrameDataSource data_fds(hu_raw,hv_raw,hw_raw,T_bad,T_lf,Trun,gds);

  WCP2dToy::uBooNEData2DDeconvolutionFDS wien_fds(hu_decon,hv_decon,hw_decon,T_bad, gds);
  ChirpMap& uplane_map = wien_fds.get_u_cmap();
  ChirpMap& vplane_map = wien_fds.get_v_cmap();
  ChirpMap& wplane_map = wien_fds.get_w_cmap();
  std::set<int>& lf_noisy_channels = data_fds.get_lf_noisy_channels();

  int run_no;
  int subrun_no;
  int event_no;
 
  Trun->SetBranchAddress("eventNo",&event_no);
  Trun->SetBranchAddress("runNo",&run_no);
  Trun->SetBranchAddress("subRunNo",&subrun_no);
  Trun->GetEntry(0);


  //std::cout << uplane_map.size() << " " << vplane_map.size() << " " << wplane_map.size() <<  " " << uplane_map[880].first << " " << uplane_map[880].second << std::endl;
  int rebin = 4;
  WCP2dToy::uBooNEDataROI uboone_rois(data_fds,wien_fds,gds,uplane_map,vplane_map,wplane_map,lf_noisy_channels);
  WCP2dToy::uBooNEDataAfterROI roi_fds(wien_fds,gds,uboone_rois,rebin);
  roi_fds.jump(0);
  
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));
  

  Int_t nwire_u = wires_u.size();
  Int_t nwire_v = wires_v.size();
  Int_t nwire_w = wires_w.size();

  int  total_time_bin = hu_decon->GetNbinsY();

  TFile *file1 = new TFile(Form("nsp3_%d_%d_%d.root",run_no,subrun_no,event_no),"RECREATE");

  TH2F *hu_roi = new TH2F("hu_roi","hu_roi",nwire_u,-0.5,nwire_u-0.5,total_time_bin/rebin,0,total_time_bin);
  TH2F *hv_roi = new TH2F("hv_roi","hv_roi",nwire_v,-0.5+nwire_u,nwire_v-0.5+nwire_u,total_time_bin/rebin,0,total_time_bin);
  TH2F *hw_roi = new TH2F("hw_roi","hw_roi",nwire_w,-0.5+nwire_u+nwire_v,nwire_w-0.5+nwire_u+nwire_v,total_time_bin/rebin,0,total_time_bin);
  TH2F *htemp1;
  
  const Frame& frame1 = roi_fds.get();
  int ntraces = frame1.traces.size();
  for (size_t ind=0; ind<ntraces; ++ind) {
    const Trace& trace = frame1.traces[ind];
    int tbin = trace.tbin;
    int chid = trace.chid;
    int nbins = trace.charge.size();
    WirePlaneType_t plane = gds.by_channel(chid).at(0)->plane();
    if (plane == WirePlaneType_t(0)){
      htemp1 = hu_roi;
    }else if (plane == WirePlaneType_t(1)){
      htemp1 = hv_roi;
      chid -= nwire_u;
    }else if (plane == WirePlaneType_t(2)){
      htemp1 = hw_roi;
      chid -= nwire_u + nwire_v;
    }
     for (int i = tbin;i!=tbin+nbins;i++){
      int tt = i+1;
      htemp1->SetBinContent(chid+1,tt,trace.charge.at(i));
    }
  }
  
  file1->Write();
  file1->Close();
  


  // std::vector<std::pair<int,int>>& rois = uboone_rois.get_self_rois(500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " S " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }
  // rois = uboone_rois.get_others_rois(500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " O " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }
  // rois = uboone_rois.get_combined_rois(500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " C " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }

  // rois = uboone_rois.get_self_rois(3500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " S " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }
  // rois = uboone_rois.get_others_rois(3500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " O " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }
  // rois = uboone_rois.get_combined_rois(3500);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " C " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }


  // rois = uboone_rois.get_combined_rois(7000);
  // std::cout << rois.size() << std::endl;
  // for (int i=0;i!=rois.size();i++){
  //   std::cout << i << " C " << rois.at(i).first << " " << rois.at(i).second << std::endl;
  // }


  
  
}
