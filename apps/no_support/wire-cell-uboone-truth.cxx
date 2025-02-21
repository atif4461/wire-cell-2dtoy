#include "WCPSst/GeomDataSource.h"
//#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCPData/MergeGeomCell.h"
#include "WCPData/GeomCluster.h"
//#include "WCPNav/SliceDataSource.h"

#include "WCPNav/FrameDataSource.h"
#include "WCPNav/SimDataSource.h"
#include "WCPNav/SliceDataSource.h"
#include "WCPSst/Util.h"
#include "WCPData/SimTruth.h"
#include "WCP2dToy/ToyDepositor.h"
#include "WCPNav/GenerativeFDS.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TColor.h"
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
 * @brief Check command line arguments
 *
 * @details Verify the number of command line arguments is sufficient
 
if (argc < 3) {
   * @brief Print usage message
 *
 * @details Display instructions for correct program invocation
 
  cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
  return 1;
}
 * @brief Initialize geometry data source
 *
 * @param argv[1] Path to ChannelWireGeometry.txt file
 
WCPSst::GeomDataSource gds(argv[1]);
 * @brief Get extent of geometry data
 *
 * @return Vector of doubles representing extent
 
std::vector<double> ex = gds.extent();
 * @brief Print extent of geometry data
 *
 * @details Display dimensions of geometry data
 
cerr << "Extent: "
     << " x:" << ex[0]/units::mm << " mm"
     << " y:" << ex[1]/units::m << " m"
     << " z:" << ex[2]/units::m << " m"
     << endl;
 * @brief Set path to root file
 *
 * @param argv[2] Path to celltree.root file
 
const char* root_file = argv[2];
 * @brief Set tree path
 *
 * @details Specify location of tree within root file
 
const char* tpath = "/Event/Sim";
 * @brief Create frame data source
 *
 * @param root_file Path to root file
 * @return Pointer to FrameDataSource object
 
WCP::FrameDataSource* fds = WCPSst::make_fds(root_file);
 * @brief Check if frame data source creation was successful
 *
 * @details Handle error case where frame data source could not be created
 
if (!fds) {
  cerr << "ERROR: failed to get FDS from " << root_file << endl;
  return 1;
}
 * @brief Define maximum number of events
 *
 * @details Specify upper limit for event processing
 
int max_events = 100;
 * @brief Define event number
 *
 * @details Identify specific event for processing
 
int eve_num = 4;
 * @brief Create ToyDepositor object
 *
 * @param fds Pointer to FrameDataSource object
 
WCP::ToyDepositor toydep(fds);
 * @brief Get depositions for specified event
 *
 * @param eve_num Event number
 * @return PointValueVector containing depositions
 
const PointValueVector pvv = toydep.depositions(eve_num);
 * @brief Create GenerativeFDS object
 *
 * @param toydep ToyDepositor object
 * @param gds Geometry data source
 * @param max_events Maximum number of events
 * @param threshold Threshold value
 
WCP::GenerativeFDS gfds(toydep,gds,2400,max_events,2.0*1.6*units::millimeter);
 * @brief Jump to specified event
 *
 * @param eve_num Event number
 
gfds.jump(eve_num);
 * @brief Create ToyuBoNESliceDataSource object
 *
 * @param gfds GenerativeFDS object
 * @param threshold Threshold value
 
WCPSst::ToyuBooNESliceDataSource sds(gfds,1500);
 * @brief Define arrays for storing coordinates
 *
 * @details Allocate memory for coordinate storage
 
const int N = 100000;
Double_t x[N],y[N],z[N];
Double_t xt[N],yt[N],zt[N];
 * @brief Initialize counters
 *
 * @details Reset counters for tracking purposes
 
int ncount = 0;
int ncount_t = 0;
 * @brief Create arrays for tiling objects
 *
 * @details Allocate memory for tiling objects
 
WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];
 * @brief Define cluster sets
 *
 * @details Initialize sets for clustering purposes
 
GeomClusterSet cluster_set, cluster_delset;
 * @brief Initialize counter for merged cells
 *
 * @details Track number of merged cells
 
int ncount_mcell = 0;
 * @brief Loop through events
 *
 * @details Process each event in sequence
 
for (int i=1191;i<2400;i++) {
 * @brief Jump to current event
 *
 * @param i Event index
 
  sds.jump(i);
 * @brief Get slice for current event
 *
 * @return Slice object
 
  WCP::Slice slice = sds.get();
 * @brief Check if slice has groups
 *
 * @details Verify presence of groups in slice
 
  if ( slice.group().size() >0){
 * @brief Print event index and group size
 *
 * @details Display event information
 
    cout << i << " " << slice.group().size() << endl;
 * @brief Create ToyTiling object
 *
 * @param slice Slice object
 * @param gds Geometry data source
 
    toytiling[i] = new WCP2dToy::ToyTiling(slice,gds);
 * @brief Get all cells for ToyTiling object
 *
 * @return GeomCellSelection object
 
    GeomCellSelection allcell = toytiling[i]->get_allcell();
 * @brief Print number of cells
 *
 * @details Display cell count
 
    cout << allcell.size() << endl;
 * @brief Create MergeToyTiling object
 *
 * @param toytiling ToyTiling object
 * @param i Event index
 * @param threshold Threshold value
 
    mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i,3,1);
 * @brief Create TruthToyTiling object
 *
 * @param toytiling ToyTiling object
 * @param pvv PointValueVector object
 * @param i Event index
 * @param gds Geometry data source
 
    truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds);
 * @brief Get all merged cells for MergeToyTiling object
 *
 * @return GeomCellSelection object
 
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
 * @brief Get all wires for MergeToyTiling object
 *
 * @return GeomWireSelection object
 
    GeomWireSelection allmwire = mergetiling[i]->get_allwire();
 * @brief Print number of merged cells and wires
 *
 * @details Display merged cell and wire counts
 
    cout << allmcell.size() << endl;
 * @brief Update cluster set
 *
 * @details Modify cluster set based on merged cells
 
    if (cluster_set.empty()){
      for (int j=0;j!=allmcell.size();j++){
        GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
        cluster_set.insert(cluster);
      }
    } else {
      for (int j=0;j!=allmcell.size();j++){
        int flag = 0;
        int flag_save = 0;
        GeomCluster *cluster_save = 0;
        cluster_delset.clear();

        for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
          flag += (*it)->AddCell(*((MergeGeomCell*)allmcell[j]));
          if (flag==1 && flag!= flag_save){
            cluster_save = *it;
          } else if (flag>1 && flag!= flag_save){
            cluster_save->MergeCluster(*(*it));
            cluster_delset.insert(*it);
          }
          flag_save = flag;
        }

        for (auto it = cluster_delset.begin();it!=cluster_delset.end();it++){
          cluster_set.erase(*it);
          delete (*it);
        }

        if (flag==0){
          GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
          cluster_set.insert(cluster);
        }
      }
    }
 * @brief Store cell coordinates
 *
 * @details Save cell center coordinates
 
    for (int j=0;j!=allcell.size();j++){
      Point p = allcell[j]->center();
      x[ncount] = i*0.32;
      y[ncount] = p.y/units::cm;
      z[ncount] = p.z/units::cm;
      ncount ++;
    }
 * @brief Calculate total number of merged cells across clusters
 *
 * @details Sum up merged cell counts from all clusters
 
    int ncount_mcell_cluster = 0;
    for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
      ncount_mcell_cluster += (*it)->get_allcell().size();
    }
    ncount_mcell += allmcell.size();
 * @brief Print summary statistics
 *
 * @details Display various counts and metrics
 
    cout << i << " " << allcell.size() << " " << allmcell.size() << " " << cluster_set.size()  << endl;
 * @brief Get charge map from TruthToyTiling object
 *
 * @return CellChargeMap object
 
    CellChargeMap ccmap = truthtiling[i]->ccmap();
 * @brief Find minimum and maximum charges
 *
 * @details Determine range of charge values
 
    Double_t charge_min = 10000;
    Double_t charge_max = 0;
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      double charge = it->second;
      if (charge > charge_max) charge_max = charge;
      if (charge < charge_min) charge_min = charge;
    }
 * @brief Store truth cell coordinates
 *
 * @details Save coordinates of truth cells
 
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      Point p = it->first->center();
      xt[ncount_t] = i*0.32;
      yt[ncount_t] = p.y/units::cm;
      zt[ncount_t] = p.z/units::cm;
      ncount_t ++;
    }
  }
}
 * @brief Create application object
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 
TApplication theApp("theApp",&argc,argv);
 * @brief Set application return status
 *
 * @details Configure application to return from Run method
 
theApp.SetReturnFromRun(true);
 * @brief Create canvas object
 *
 * @param title Canvas title
 * @param name Canvas name
 * @param width Canvas width
 * @param height Canvas height
 
TCanvas c1("ToyMC","ToyMC",800,600);
 * @brief Draw canvas
 *
 * @details Render the canvas
 
c1.Draw();
 * @brief Create ToyEventDisplay object
 *
 * @param c1 Canvas object
 * @param gds Geometry data source
 
WCP2dToy::ToyEventDisplay display(c1, gds);
 * @brief Set charge range for display
 *
 * @param charge_min Minimum charge
 * @param charge_max Maximum charge
 
display.charge_min = charge_min;
display.charge_max = charge_max;
 * @brief Configure graphics style
 *
 * @details Customize appearance of plots
 
gStyle->SetOptStat(0);
 * @brief Define color palette
 *
 * @param NRGBs Number of colors in palette
 * @param NCont Number of contours
 
const Int_t NRGBs = 5;
const Int_t NCont = 255;
Int_t MyPalette[NCont];
Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
Double_t red[NRGBs] = {0.0, 0.0, 0.87,1.0, 0.51};
Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2,0.0};
Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
 * @brief Create gradient color table
 *
 * @param NRGBs Number of colors in palette
 * @param stops Color stop points
 * @param red Red component values
 * @param green Green component values
 * @param blue Blue component values
 * @param NCont Number of contours
 
Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
 * @brief Set number of contours
 *
 * @param NCont Number of contours
 
gStyle->SetNumberContours(NCont);
 * @brief Set palette
 *
 * @param NCont Number of contours
 * @param MyPalette Palette array
 
for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
gStyle->SetPalette(NCont,MyPalette);
 * @brief Initialize display
 *
 * @param arg1 Initial value 1
 * @param arg2 Initial value 2
 * @param arg3 Initial value 3
 * @param arg4 Initial value 4
 
display.init(0,10.3698,-2.33/2.,2.33/2.);
 * @brief Draw Monte Carlo data
 *
 * @param arg1 Data identifier
 * @param arg2 Data type
 
display.draw_mc(1,WCP::PointValueVector(),"colz");
 * @brief Draw slice data
 *
 * @param slice Slice object
 
display.draw_slice(slice,"");
 * @brief Draw cells
 *
 * @param allcell GeomCellSelection object
 
display.draw_cells(toytiling[i]->get_allcell(),"*same");
 * @brief Draw merged cells
 *
 * @param allmcell GeomCellSelection object
 * @param flag Flag value
 
display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0);
 * @brief Draw truth cells
 *
 * @param ccmap CellChargeMap object
 
display.draw_truthcells(ccmap,"*same");
 * @brief Run application
 *
 * @details Execute the application's main loop
 
theApp.Run();
}
 * @brief Calculate final merged cell count across clusters
 *
 * @details Sum up merged cell counts from all clusters
 
int ncount_mcell_cluster = 0;
for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  ncount_mcell_cluster += (*it)->get_allcell().size();
}
 * @brief Print summary statistics
 *
 * @details Display final counts and metrics
 
cout << "Summary: " << ncount << " " << ncount_mcell << " " << ncount_mcell_cluster << endl;
 * @brief Return program exit status
 *
 * @details Indicate successful execution
 
return 0;* This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.2.
*/ 
int main(int argc, char* argv[])
{
  if (argc < 3) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
    return 1;
  }
  
  
  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;
  
  
  const char* root_file = argv[2];
  const char* tpath = "/Event/Sim";
  
  WCP::FrameDataSource* fds = 0;
  fds = WCPSst::make_fds(root_file);
  if (!fds) {
    cerr << "ERROR: failed to get FDS from " << root_file << endl;
    return 1;
  }
  
  int max_events = 100;
  int eve_num = 4;

  
  WCP::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(eve_num);
  

  WCP::GenerativeFDS gfds(toydep,gds,2400,max_events,2.0*1.6*units::millimeter);
  gfds.jump(eve_num);

  WCPSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons

  

  const int N = 100000;
  Double_t x[N],y[N],z[N];
  Double_t xt[N],yt[N],zt[N];
  int ncount = 0;
  int ncount_t = 0;
  

  WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
  WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
  WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];
  
  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;
  

  //int i=178;{
  int i=1191;{
    //int i=441;{
    // for (int i=0;i!=sds.size();i++){
    //for (int i=365;i!=378;i++){
 
    sds.jump(i);
    WCP::Slice slice = sds.get();
    if ( slice.group().size() >0){
      cout << i << " " << slice.group().size() << endl;

      toytiling[i] = new WCP2dToy::ToyTiling(slice,gds);
      GeomCellSelection allcell = toytiling[i]->get_allcell();

      cout << allcell.size() << endl;
      
      mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i,3,1);
      truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds);
      
      

      // for (int j=0;j!=allcell.size();j++){
      // 	std::cout << toytiling[i]->wires(*allcell.at(j)).size() << std::endl;
      // }
      

      GeomCellSelection allmcell = mergetiling[i]->get_allcell();
      GeomWireSelection allmwire = mergetiling[i]->get_allwire();
      
      cout << allmcell.size() << endl;
     
      if (cluster_set.empty()){
  	// if cluster is empty, just insert all the mcell, each as a cluster
       	for (int j=0;j!=allmcell.size();j++){
  	  GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
  	  cluster_set.insert(cluster);
  	}
      }else{
  	for (int j=0;j!=allmcell.size();j++){
  	  int flag = 0;
  	  int flag_save = 0;
  	  GeomCluster *cluster_save = 0;
	  
  	  cluster_delset.clear();

  	  // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	  //   if (i==318)
  	  //     cout << "b " << (*it)->get_allcell().size() << endl;
  	  // } 
	  

  	  // loop through merged cell
  	  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	    //loop through clusters
	   
  	    flag += (*it)->AddCell(*((MergeGeomCell*)allmcell[j]));
  	    if (flag==1 && flag != flag_save){
  	      cluster_save = *it;
  	    }else if (flag>1 && flag != flag_save){
  	      cluster_save->MergeCluster(*(*it));
  	      cluster_delset.insert(*it);
  	    }
  	    flag_save = flag;
  	    // if (i==318)
  	    //   cout << "c " << flag << endl;
  	  }

  	  for (auto it = cluster_delset.begin();it!=cluster_delset.end();it++){
  	    cluster_set.erase(*it);
  	    delete (*it);
  	  }
	  
	  

  	  // if (i==318)
  	  //   cout << j << " " << flag << endl;
  	  if (flag==0){
  	    GeomCluster *cluster = new GeomCluster(*((MergeGeomCell*)allmcell[j]));
  	    cluster_set.insert(cluster);
  	  }

  	  // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	  //   if (i==318)
  	  //     cout << (*it)->get_allcell().size() << endl;
  	  // }
	  
  	}
      }
      

      for (int j=0;j!=allcell.size();j++){
  	Point p = allcell[j]->center();
  	x[ncount] = i*0.32;
  	y[ncount] = p.y/units::cm;
  	z[ncount] = p.z/units::cm;
  	ncount ++;
      }


      int ncount_mcell_cluster = 0;
      for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  	ncount_mcell_cluster += (*it)->get_allcell().size();
      }
      ncount_mcell += allmcell.size();
      
      
      int ncells_qx = 0;
      for (int j=0;j!=allmcell.size();j++){
	MergeGeomCell* mcell = (MergeGeomCell*)allmcell.at(j);
	ncells_qx += mcell->get_allcell().size();
      }

      cout << i << " " << allcell.size() << " " << ncells_qx << " " << allmcell.size() << " " << allmwire.size() << " " << cluster_set.size()  << endl;
      

      // for (int j=0;j!=allmwire.size();j++){
      // 	cout << mergetiling.cells(*allmwire[j]).size() << endl;
      // }


      CellChargeMap ccmap = truthtiling[i]->ccmap();

      Double_t charge_min = 10000;
      Double_t charge_max = 0;

      for (auto it = ccmap.begin();it!=ccmap.end(); it++){
  	Point p = it->first->center();
      	xt[ncount_t] = i*0.32;
      	yt[ncount_t] = p.y/units::cm;
      	zt[ncount_t] = p.z/units::cm;
      	ncount_t ++;

  	double charge = it->second;
  	if (charge > charge_max) charge_max = charge;
  	if (charge < charge_min) charge_min = charge;
       	// cout << it->second << endl;
      }


      //loop through merged cell and compare with truth cells
      for (int j=0;j!=allmcell.size();j++){
  	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
  	mcell->CheckContainTruthCell(ccmap);
  	// 	cout << mergetiling.wires(*allmcell[j]).size() << endl;
      }


      
      // WireChargeMap wcmap = toytiling.wcmap();
      // for (auto it = wcmap.begin();it!=wcmap.end(); it++){
      // 	double charge = it->second;
      // 	if (charge > charge_max) charge_max = charge;
      // 	if (charge < charge_min) charge_min = charge;
      // }
    // int sum = 0;
    // for (int j=0;j!=allmcell.size();j++){
    //   sum += ((WCP::MergeGeomCell*)allmcell[j])->get_allcell().size() ;
    // }
    // cout << allcell.size() << " " << allmcell.size() << " "  << sum << endl;

     // for (int j=0;j!=allmcell.size();j++){
     //   cout << j << " " << mergetiling.wires(*allmcell[j]).size() << endl;
     // }

    // for (int j=0;j!=allwire.size();j++){
    //   const GeomWire *wire = allwire[j];
    //   const GeomCellSelection targetcells =  mergetiling.cells(*wire);
    //   cout << j << " " << targetcells.size() << endl;
    // }

    

    // GeomCellSelection allcell = toytiling.get_allcell();
    // GeomWireSelection allwire = toytiling.get_allwire();
    // cout << i << " " << allcell.size() << " " << allwire.size() << endl;
    //
    //}
  

  // cout << toytiling.wiremap[allwire.at(0)].size() << endl;
  // cout << toytiling.cellmap[allcell.at(0)].size() << endl;
  // //  

      
    TApplication theApp("theApp",&argc,argv);
    theApp.SetReturnFromRun(true);
    
    TCanvas c1("ToyMC","ToyMC",800,600);
    c1.Draw();
    
    WCP2dToy::ToyEventDisplay display(c1, gds);
    display.charge_min = charge_min;
    display.charge_max = charge_max;


    gStyle->SetOptStat(0);
    
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Int_t MyPalette[NCont];
    Double_t stops[NRGBs] = {0.0, 0.34, 0.61, 0.84, 1.0};
    Double_t red[NRGBs] = {0.0, 0.0, 0.87 ,1.0, 0.51};
    Double_t green[NRGBs] = {0.0, 0.81, 1.0, 0.2 ,0.0};
    Double_t blue[NRGBs] = {0.51, 1.0, 0.12, 0.0, 0.0};
    Int_t FI = TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);
    for (int kk=0;kk!=NCont;kk++) MyPalette[kk] = FI+kk;
    gStyle->SetPalette(NCont,MyPalette);

    

    display.init(0,10.3698,-2.33/2.,2.33/2.);
    //display.init(1.1,1.8,0.7,1.0);
    
    display.draw_mc(1,WCP::PointValueVector(),"colz");
    
    

    display.draw_slice(slice,"");
    display.draw_cells(toytiling[i]->get_allcell(),"*same");
    display.draw_mergecells(mergetiling[i]->get_allcell(),"*same",0); //0 is normal, 1 is only draw the ones containt the truth cell
    display.draw_truthcells(ccmap,"*same");
    
    // display.draw_wires_charge(wcmap,"Fsame",FI);
    // display.draw_cells_charge(toytiling.get_allcell(),"Fsame");
    // display.draw_truthcells_charge(ccmap,"lFsame",FI);
    
    
    theApp.Run();
    }
  }

  int ncount_mcell_cluster = 0;
  for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
    ncount_mcell_cluster += (*it)->get_allcell().size();
  }


  cout << "Summary: " << ncount << " " << ncount_mcell << " " << ncount_mcell_cluster << endl;
  // TGraph2D *g = new TGraph2D(ncount,x,y,z);
  // TGraph2D *gt = new TGraph2D(ncount_t,xt,yt,zt);
  // TFile *file = new TFile("shower3D.root","RECREATE");
  // g->Write("shower3D");
  // gt->Write("shower3D_truth");

  // //save cluster
  // int ncluster = 0;
  // for (auto it = cluster_set.begin();it!=cluster_set.end();it++){
  //   ncount = 0;
  //   for (int i=0; i!=(*it)->get_allcell().size();i++){
  //     const MergeGeomCell *mcell = (const MergeGeomCell*)((*it)->get_allcell().at(i));
  //     for (int j=0; j!=mcell->get_allcell().size();j++){
  // 	Point p = mcell->get_allcell().at(j)->center();
  // 	x[ncount] = mcell->GetTimeSlice()*0.32;
  // 	y[ncount] = p.y/units::cm;
  // 	z[ncount] = p.z/units::cm;
  // 	ncount ++;
  //     }
  //   }
  //   // cout << ncount << endl;
  //   TGraph2D *g1 = new TGraph2D(ncount,x,y,z);
  //   g1->Write(Form("cluster_%d",ncluster));
  //   ncluster ++;
  // }

  // file->Write();
  // file->Close();

  return 0;
  
} // main()
