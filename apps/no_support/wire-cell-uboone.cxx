#include "WCPSst/GeomDataSource.h"
#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"

#include "WCPData/MergeGeomCell.h"
//#include "WCPNav/SliceDataSource.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
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
 * @param argc Number of command line arguments
 * @return True if valid, false otherwise
 
if (argc < 3) 
 * @brief Print usage message
 *
 * @param None
 * @return None
 
cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
 * @brief Return non-zero exit status
 *
 * @param None
 * @return Non-zero integer value
 
return 1;
 * @brief Create geometry data source object
 *
 * @param filename Geometry file path
 * @return Geometry data source object
 
WCPSst::GeomDataSource gds(argv[1]);
 * @brief Get extent of geometry data source
 *
 * @param gds Geometry data source object
 * @return Vector of doubles representing extent
 
std::vector<double> ex = gds.extent();
 * @brief Print extent values
 *
 * @param ex Vector of doubles representing extent
 * @return None
 
cerr << "Extent: "
     << " x:" << ex[0]/units::mm << " mm"
     << " y:" << ex[1]/units::m << " m"
     << " z:" << ex[2]/units::m << " m"
     << endl;
 * @brief Open ROOT file
 *
 * @param filename ROOT file path
 * @param mode File open mode
 * @return TFile object
 
TFile tfile(root_file,"read");
 * @brief Get tree from ROOT file
 *
 * @param tfile TFile object
 * @param path Tree path
 * @return TTree object
 
TTree* sst = dynamic_cast<TTree*>(tfile.Get(tpath));
 * @brief Create frame data source object
 *
 * @param sst TTree object
 * @param gds Geometry data source object
 * @return Frame data source object
 
WCPSst::ToyuBooNEFrameDataSource fds(*sst,gds);
 * @brief Print number of frames
 *
 * @param fds Frame data source object
 * @return None
 
std::cerr << "Got " << fds.size() 
          << " frames from " << tpath 
          << " in " << root_file << std::endl;
 * @brief Jump to frame index
 *
 * @param fds Frame data source object
 * @param index Frame index
 * @return None
 
fds.jump(1);
 * @brief Get frame object
 *
 * @param fds Frame data source object
 * @return Frame object
 
WCP::Frame frame = fds.get();
 * @brief Create slice data source object
 *
 * @param fds Frame data source object
 * @param index Slice index
 * @return Slice data source object
 
WCPSst::ToyuBooNESliceDataSource sds(fds,1);
 * @brief Loop over slices
 *
 * @param sds Slice data source object
 * @return None
 
//for (int i=0;i!=sds.size();i++){
 * @brief Jump to slice index
 *
 * @param sds Slice data source object
 * @param index Slice index
 * @return None
 
sds.jump(i);
 * @brief Get slice object
 *
 * @param sds Slice data source object
 * @return Slice object
 
WCP::Slice slice = sds.get();
 * @brief Check if group size is greater than zero
 *
 * @param slice Slice object
 * @return True if valid, false otherwise
 
if ( slice.group().size() >0){
 * @brief Create ToyTiling object
 *
 * @param slice Slice object
 * @param gds Geometry data source object
 * @return ToyTiling object
 
WCP2dToy::ToyTiling toytiling(slice,gds);
 * @brief Get all cell selection
 *
 * @param toytiling ToyTiling object
 * @return Cell selection object
 
GeomCellSelection allcell = toytiling.get_allcell();
 * @brief Create MergeToyTiling object
 *
 * @param toytiling ToyTiling object
 * @return MergeToyTiling object
 
WCP2dToy::MergeToyTiling mergetiling(toytiling);    
 * @brief Get all merge cell selection
 *
 * @param mergetiling MergeToyTiling object
 * @return Cell selection object
 
GeomCellSelection allmcell = mergetiling.get_allcell();
 * @brief Get all wire selection
 *
 * @param mergetiling MergeToyTiling object
 * @return Wire selection object
 
GeomWireSelection allwire = mergetiling.get_allwire();
 * @brief Print slice index and selection sizes
 *
 * @param i Slice index
 * @param allmcell Cell selection object
 * @param allwire Wire selection object
 * @return None
 
if (allcell.size()>0){
  cout << i << " " << allmcell.size() << " "  << allwire.size() << endl;
}
 * @brief Create application object
 *
 * @param argc Argument count
 * @param argv Argument vector
 * @return Application object
 
TApplication theApp("theApp",&argc,argv);
 * @brief Set application return flag
 *
 * @param theApp Application object
 * @return None
 
theApp.SetReturnFromRun(true);
 * @brief Create canvas object
 *
 * @param title Canvas title
 * @param name Canvas name
 * @param width Canvas width
 * @param height Canvas height
 * @return Canvas object
 
TCanvas c1("ToyMC","ToyMC",800,600);
 * @brief Draw canvas
 *
 * @param c1 Canvas object
 * @return None
 
c1.Draw();
 * @brief Create event display object
 *
 * @param c1 Canvas object
 * @param gds Geometry data source object
 * @return Event display object
 
WCP2dToy::ToyEventDisplay display(c1, gds);
 * @brief Set graphics style option
 *
 * @param None
 * @return None
 
gStyle->SetOptStat(0);
 * @brief Initialize event display
 *
 * @param display Event display object
 * @param parameters Initialization parameters
 * @return None
 
display.init(0,10.3698,-2.33/2.,2.33/2.);
 * @brief Draw Monte Carlo truth
 *
 * @param display Event display object
 * @param parameter Parameter
 * @param label Label string
 * @return None
 
display.draw_mc(1,WCP::PointValueVector(),"");
 * @brief Draw slice
 *
 * @param display Event display object
 * @param slice Slice object
 * @param label Label string
 * @return None
 
display.draw_slice(slice,"");
 * @brief Draw cells
 *
 * @param display Event display object
 * @param cells Cell selection object
 * @param label Label string
 * @return None
 
display.draw_cells(toytiling.get_allcell(),"*same");
 * @brief Draw merged cells
 *
 * @param display Event display object
 * @param cells Cell selection object
 * @param label Label string
 * @return None
 
display.draw_mergecells(mergetiling.get_allcell(),"*same");
 * @brief Run application
 *
 * @param theApp Application object
 * @return None
 
theApp.Run();
 * @brief Return zero exit status
 *
 * @param None
 * @return Zero integer value
 
return 0;* This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.01.
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
  
  TFile tfile(root_file,"read");
  TTree* sst = dynamic_cast<TTree*>(tfile.Get(tpath));
  WCPSst::ToyuBooNEFrameDataSource fds(*sst,gds);
  std::cerr << "Got " << fds.size() 
	    << " frames from " << tpath 
	    << " in " << root_file << std::endl;
  
  fds.jump(1);
  WCP::Frame frame = fds.get();
  
  WCPSst::ToyuBooNESliceDataSource sds(fds,1);
  
  // int i=1129;{
  int i=331;{
  //for (int i=0;i!=sds.size();i++){
  //for (int i=1143;i!=1145;i++){
  sds.jump(i);
  WCP::Slice slice = sds.get();
  if ( slice.group().size() >0){
    WCP2dToy::ToyTiling toytiling(slice,gds);
    GeomCellSelection allcell = toytiling.get_allcell();

   
    

    WCP2dToy::MergeToyTiling mergetiling(toytiling);    
    GeomCellSelection allmcell = mergetiling.get_allcell();
    GeomWireSelection allwire = mergetiling.get_allwire();
    //cout << i << endl;
    // int sum = 0;
    // for (int j=0;j!=allmcell.size();j++){
    //   GeomCellSelection cells = ((WCP::MergeGeomCell*)allmcell[j])->get_allcell();
    //   for (int k=0;k!=cells.size();k++){
    // 	GeomWireSelection wires = toytiling.wires(*cells[k]);
    // 	if( wires[0]->ident()==0 || wires[1]->ident()==0 || wires[2]->ident()==0){
    // 	  cout << i << " Wrong!!" << endl;
    // 	}
    //   }
    //   //   sum += ((WCP::MergeGeomCell*)allmcell[j])->get_allcell().size() ;
    // }

    //  //debug the toytiling itself. 
    // for (int j=0;j!=allcell.size();j++){
    //   GeomWireSelection wires = toytiling.wires(*allcell[j]);
    //   if( wires[0]->ident()==0 || wires[1]->ident()==0 || wires[2]->ident()==0){
    // 	cout << i << " Wrong!!" << endl;
    //   }  
    // }



    if (allcell.size()>0){
      cout << i << " " << allmcell.size() << " "  << allwire.size() << endl;
    }

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
    
    gStyle->SetOptStat(0);
    
    display.init(0,10.3698,-2.33/2.,2.33/2.);
    //display.init(0.6,1.0,0.0,0.3);
    //display.init(0.6,0.7,0.07,0.12);
    //display.init();
    display.draw_mc(1,WCP::PointValueVector(),"");
    //display.draw_mc(1,fds.mctruth,"");
    //display.draw_mc(2,fds.mctruth,"TEXT");
    
    
    display.draw_slice(slice,"");
    
    display.draw_cells(toytiling.get_allcell(),"*same");
    display.draw_mergecells(mergetiling.get_allcell(),"*same");
    //display.draw_mc(3,fds.mctruth,"*same");
    
    theApp.Run();
  }
  }

  return 0;
  
} // main()
