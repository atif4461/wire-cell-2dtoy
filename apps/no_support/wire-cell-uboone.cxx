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
 
if (argc < 3) 

 * @brief Print usage message
 * 
 
cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl

 * @brief Return non-zero exit status
 * 
 
return 1

 * @brief Create geometry data source object
 * 
 * @param argv Command line argument array
 
WCPSst::GeomDataSource gds(argv[1])

 * @brief Get extent of geometry data source
 * 
 * @return Vector of doubles representing extent
 
std::vector<double> ex = gds.extent()

 * @brief Print extent to error stream
 * 
 
cerr << "Extent: "
<< " x:" << ex[0]/units::mm << " mm"
<< " y:" << ex[1]/units::m << " m"
<< " z:" << ex[2]/units::m << " m"
<< endl

 * @brief Set root file path from command line argument
 * 
 
const char* root_file = argv[2]

 * @brief Set tree path
 * 
 
const char* tpath = "/Event/Sim"

 * @brief Open root file
 * 
 * @param root_file Root file path
 * @param mode File open mode
 
TFile tfile(root_file,"read")

 * @brief Get tree from root file
 * 
 * @param tfile Root file object
 * @param tpath Tree path
 * @return Pointer to tree object
 
TTree* sst = dynamic_cast<TTree*>(tfile.Get(tpath))

 * @brief Create frame data source object
 * 
 * @param sst Tree object pointer
 * @param gds Geometry data source object
 
WCPSst::ToyuBooNEFrameDataSource fds(*sst,gds)

 * @brief Print number of frames to error stream
 * 
 
std::cerr << "Got " << fds.size() 
<< " frames from " << tpath 
<< " in " << root_file << std::endl

 * @brief Jump to frame index
 * 
 * @param index Frame index
 
fds.jump(1)

 * @brief Get frame object
 * 
 * @return Frame object
 
WCP::Frame frame = fds.get()

 * @brief Create slice data source object
 * 
 * @param fds Frame data source object
 * @param index Slice index
 
WCPSst::ToyuBooNESliceDataSource sds(fds,1)

 * @brief Jump to slice index
 * 
 * @param i Slice index
 
sds.jump(i)

 * @brief Get slice object
 * 
 * @return Slice object
 
WCP::Slice slice = sds.get()

 * @brief Check if group size is greater than zero
 * 
 
if ( slice.group().size() >0)

 * @brief Create ToyTiling object
 * 
 * @param slice Slice object
 * @param gds Geometry data source object
 
WCP2dToy::ToyTiling toytiling(slice,gds)

 * @brief Get all cell selection
 * 
 * @return Cell selection object
 
GeomCellSelection allcell = toytiling.get_allcell()

 * @brief Create MergeToyTiling object
 * 
 * @param toytiling ToyTiling object
 
WCP2dToy::MergeToyTiling mergetiling(toytiling)

 * @brief Get all merge cell selection
 * 
 * @return Merge cell selection object
 
GeomCellSelection allmcell = mergetiling.get_allcell()

 * @brief Get all wire selection
 * 
 * @return Wire selection object
 
GeomWireSelection allwire = mergetiling.get_allwire()

 * @brief Print slice index and selection sizes
 * 
 
if (allcell.size()>0)
cout << i << " " << allmcell.size() << " "  << allwire.size() << endl

 * @brief Create application object
 * 
 * @param argc Argument count
 * @param argv Argument vector
 
TApplication theApp("theApp",&argc,argv)

 * @brief Set return from run flag
 * 
 
theApp.SetReturnFromRun(true)

 * @brief Create canvas object
 * 
 * @param title Canvas title
 * @param name Canvas name
 * @param width Canvas width
 * @param height Canvas height
 
TCanvas c1("ToyMC","ToyMC",800,600)

 * @brief Draw canvas
 * 
 
c1.Draw()

 * @brief Create event display object
 * 
 * @param c1 Canvas object
 * @param gds Geometry data source object
 
WCP2dToy::ToyEventDisplay display(c1, gds)

 * @brief Set opt stat
 * 
 * @param value Opt stat value
 
gStyle->SetOptStat(0)

 * @brief Initialize display
 * 
 * @param arg1 Initialization argument 1
 * @param arg2 Initialization argument 2
 * @param arg3 Initialization argument 3
 * @param arg4 Initialization argument 4
 
display.init(0,10.3698,-2.33/2.,2.33/2.)

 * @brief Draw MC
 * 
 * @param arg1 Drawing argument 1
 * @param arg2 Drawing argument 2
 * @param option Drawing option
 
display.draw_mc(1,WCP::PointValueVector(),"")

 * @brief Draw slice
 * 
 * @param slice Slice object
 * @param option Drawing option
 
display.draw_slice(slice,"")

 * @brief Draw cells
 * 
 * @param cells Cell selection object
 * @param option Drawing option
 
display.draw_cells(toytiling.get_allcell(),"*same")

 * @brief Draw merge cells
 * 
 * @param cells Merge cell selection object
 * @param option Drawing option
 
display.draw_mergecells(mergetiling.get_allcell(),"*same")

 * @brief Run application
 * 
 
theApp.Run()

 * @brief Return program exit status
 * 
 
return 0* This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.2.
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
