#include "WCP2dToy/pd_Data_FDS.h"

using namespace WCP;

/**
 * Constructor for pdDataFDS class
 * @param gds GeomDataSource object containing geometry data
 * @param hu_decon pointer to TH2I histogram for u-plane deconvolution
 * @param hv_decon pointer to TH2I histogram for v-plane deconvolution
 * @param hw_decon pointer to TH2I histogram for w-plane deconvolution
 * @param eve_num event number
 */
// The above comment was written by an LLM. 
WCP2dToy::pdDataFDS::pdDataFDS(const WCP::GeomDataSource& gds, TH2I *hu_decon, TH2I *hv_decon, TH2I *hw_decon, int eve_num)
  : gds(gds)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  
  
  frame.clear();		// win or lose, we start anew

  frame.index =eve_num;
  bins_per_frame = hu_decon->GetNbinsY();
  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {
    WCP::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hu_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    WCP::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    WCP::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  //std::cout << frame.traces.size() << " " << bins_per_frame << std::endl;
}



/**
 * Constructor for pdDataFDS class 
 * @param gds Geom data source object
 * @param hu_decon pointer to histogram of deconvoluted u plane data
 * @param hv_decon pointer to histogram of deconvoluted v plane data
 * @param hw_decon pointer to histogram of deconvoluted w plane data
 * @param eve_num event number
 */
// The above comment was written by an LLM. 
WCP2dToy::pdDataFDS::pdDataFDS(const WCP::GeomDataSource& gds, TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, int eve_num)
  : gds(gds)
{
  GeomWireSelection wires_u = gds.wires_in_plane(WirePlaneType_t(0));
  GeomWireSelection wires_v = gds.wires_in_plane(WirePlaneType_t(1));
  GeomWireSelection wires_w = gds.wires_in_plane(WirePlaneType_t(2));

  nwire_u = wires_u.size();
  nwire_v = wires_v.size();
  nwire_w = wires_w.size();

  
  
  frame.clear();		// win or lose, we start anew

  frame.index =eve_num;
  bins_per_frame = hu_decon->GetNbinsY();
  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {
    WCP::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hu_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    WCP::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    WCP::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  //std::cout << frame.traces.size() << " " << bins_per_frame << std::endl;
}

/**
 * Refreshes the pdDataFDS object with new data from the given histograms
 * @param hu_decon histogram for U plane deconvolution
 * @param hv_decon histogram for V plane deconvolution
 * @param hw_decon histogram for W plane deconvolution
 * @param eve_num event number
 */
// The above comment was written by an LLM. 
void WCP2dToy::pdDataFDS::refresh(TH2F *hu_decon, TH2F *hv_decon, TH2F *hw_decon, int eve_num){
    
  frame.clear();		// win or lose, we start anew

  frame.index =eve_num;
  bins_per_frame = hu_decon->GetNbinsY();
  // U plane
  for (size_t ind=0; ind < hu_decon->GetNbinsX(); ++ind) {
    WCP::Trace trace;
    trace.chid = ind;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hu_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }
  
  // V plane
  for (size_t ind=0; ind < hv_decon->GetNbinsX(); ++ind) {
    WCP::Trace trace;
    trace.chid = ind + nwire_u;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hv_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }

  // W plane
  for (size_t ind=0; ind < hw_decon->GetNbinsX(); ++ind) {
    WCP::Trace trace;
    trace.chid = ind + nwire_u + nwire_v;
    trace.tbin = 0;		// full readout, if zero suppress this would be non-zero
    trace.charge.resize(bins_per_frame, 0.0);
    
    for (int ibin=0; ibin != bins_per_frame; ibin++) {
      trace.charge.at(ibin) = hw_decon->GetBinContent(ind+1,ibin+1);
    }
    frame.traces.push_back(trace);
  }
}


WCP2dToy::pdDataFDS::~pdDataFDS(){
}

int WCP2dToy::pdDataFDS::jump(int frame_number){
  return frame.index;
}

int WCP2dToy::pdDataFDS::size() const{
  return 1;
}
