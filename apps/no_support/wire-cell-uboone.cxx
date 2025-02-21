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
 * @brief Main program entry point.
 *
 * This function is responsible for executing the main logic of the program.
 * It takes command line arguments, initializes necessary objects, and performs
 * actions based on those inputs.
 *
 * @param argc Number of command line arguments.
 * @param argv Array of command line argument strings.
 * @return Program exit status.
  
 * @brief Checks command line argument count and prints usage message if insufficient.
 *
 * Exits the program with an error code if fewer than 3 arguments are provided.
 *
 * @return None
  
if (argc < 3) {
}
 * @brief Initializes geometry data source object from file specified by first command line argument.
 *
 * Creates a WCPSst::GeomDataSource object using the first command line argument.
 *
 * @param argv Command line argument array
  
WCPSst::GeomDataSource gds(argv[1]);
 * @brief Retrieves extent of geometry data source.
 *
 * Gets the extent of the geometry data source as a vector of doubles.
 *
 * @return Vector of doubles representing extent
  
std::vector<double> ex = gds.extent();
 * @brief Prints extent of geometry data source to standard error stream.
 *
 * Outputs the extent of the geometry data source to the standard error stream.
 *
 * @return None
  
cerr << "Extent: "
     << " x:" << ex[0]/units::mm << " mm"
     << " y:" << ex[1]/units::m << " m"
     << " z:" << ex[2]/units::m << " m"
     << endl;
 * @brief Sets path to root file.
 *
 * Assigns the second command line argument to the root_file variable.
 *
 * @param argv Command line argument array
  
const char* root_file = argv[2];
 * @brief Opens root file and retrieves tree.
 *
 * Opens the root file in read mode and gets the tree from the file.
 *
 * @param root_file Path to root file
 * @return Pointer to TTree object
  
TFile tfile(root_file,"read");
TTree* sst = dynamic_cast<TTree*>(tfile.Get("/Event/Sim"));
 * @brief Creates ToYuBooNE Frame Data Source object.
 *
 * Constructs a WCPSst::ToyuBooNEFrameDataSource object using the tree and geometry data source.
 *
 * @param sst Pointer to TTree object
 * @param gds Geometry data source object
  
WCPSst::ToyuBooNEFrameDataSource fds(*sst,gds);
 * @brief Prints number of frames retrieved from tree.
 *
 * Outputs the number of frames retrieved from the tree to the standard error stream.
 *
 * @return None
  
std::cerr << "Got " << fds.size() 
       << " frames from " << "/Event/Sim" 
       << " in " << root_file << std::endl;
 * @brief Jumps to specific frame in frame data source.
 *
 * Moves to the frame at index 1 in the frame data source.
 *
 * @return None
  
fds.jump(1);
 * @brief Retrieves current frame from frame data source.
 *
 * Gets the current frame from the frame data source.
 *
 * @return WCP::Frame object
  
WCP::Frame frame = fds.get();
 * @brief Creates ToYuBooNE Slice Data Source object.
 *
 * Constructs a WCPSst::ToyuBooNESliceDataSource object using the frame data source and index 1.
 *
 * @param fds Frame data source object
 * @param index Index of slice
  
WCPSst::ToYuBooNESliceDataSource sds(fds,1);
 * @brief Loops through slices in slice data source.
 *
 * Iterates over the slices in the slice data source, starting from index 331.
 *
 * @return None
  
//for (int i=0;i!=sds.size();i++){
//for (int i=1143;i!=1145;i++){
int i=331;
{
 * @brief Jumps to specific slice in slice data source.
 *
 * Moves to the slice at index i in the slice data source.
 *
 * @return None
  
sds.jump(i);
 * @brief Retrieves current slice from slice data source.
 *
 * Gets the current slice from the slice data source.
 *
 * @return WCP::Slice object
  
WCP::Slice slice = sds.get();
 * @brief Checks if slice has groups.
 *
 * Verifies if the slice has any groups.
 *
 * @return Boolean value indicating presence of groups
  
if ( slice.group().size() >0){
 * @brief Creates Toy Tiling object.
 *
 * Constructs a WCP2dToy::ToyTilling object using the slice and geometry data source.
 *
 * @param slice Slice object
 * @param gds Geometry data source object
  
WCP2dToy::ToyTiling toytiling(slice,gds);
 * @brief Retrieves all cells from Toy Tiling object.
 *
 * Gets all cells from the Toy Tilling object.
 *
 * @return GeomCellSelection object
  
GeomCellSelection allcell = toytiling.get_allcell();
 * @brief Creates Merge Toy Tiling object.
 *
 * Constructs a WCP2dToy::MergeToyTiling object using the Toy Tilling object.
 *
 * @param toytiling Toy Tilling object
  
WCP2dToy::MergeToyTiling mergetiling(toytiling);
 * @brief Retrieves all merge cells from Merge Toy Tiling object.
 *
 * Gets all merge cells from the Merge Toy Tiling object.
 *
 * @return GeomCellSelection object
  
GeomCellSelection allmcell = mergetiling.get_allcell();
 * @brief Retrieves all wires from Merge Toy Tiling object.
 *
 * Gets all wires from the Merge Toy Tiling object.
 *
 * @return GeomWireSelection object
  
GeomWireSelection allwire = mergetiling.get_allwire();
 * @brief Prints size of all cell and wire selections.
 *
 * Outputs the sizes of the all cell and wire selections to the console.
 *
 * @return None
  
if (allcell.size()>0){
  cout << i << " " << allmcell.size() << " "  << allwire.size() << endl;
}
 * @brief Creates TApplication object.
 *
 * Constructs a TApplication object with the name "theApp".
 *
 * @return None
  
TApplication theApp("theApp",&argc,argv);
 * @brief Sets return from run flag.
 *
 * Enables the return from run flag for the TApplication object.
 *
 * @return None
  
theApp.SetReturnFromRun(true);
 * @brief Creates TCanvas object.
 *
 * Constructs a TCanvas object with the name "ToyMC" and dimensions 800x600.
 *
 * @return None
  
TCanvas c1("ToyMC","ToyMC",800,600);
 * @brief Draws TCanvas object.
 *
 * Draws the TCanvas object.
 *
 * @return None
  
c1.Draw();
 * @brief Creates Toy Event Display object.
 *
 * Constructs a WCP2dToy::ToyEventDisplay object using the TCanvas and geometry data source.
 *
 * @param c1 TCanvas object
 * @param gds Geometry data source object
  
WCP2dToy::ToyEventDisplay display(c1, gds);
 * @brief Sets opt stat flag.
 *
 * Disables the opt stat flag using the gStyle object.
 *
 * @return None
  
gStyle->SetOptStat(0);
 * @brief Initializes Toy Event Display object.
 *
 * Initializes the Toy Event Display object with parameters 0, 10.3698, -2.33/2., and 2.33/2..
 *
 * @return None
  
display.init(0,10.3698,-2.33/2.,2.33/2.);
 * @brief Draws MC event.
 *
 * Draws the MC event using the Toy Event Display object.
 *
 * @return None
  
display.draw_mc(1,WCP::PointValueVector(),"");
 * @brief Draws slice.
 *
 * Draws the slice using the Toy Event Display object.
 *
 * @return None
  
display.draw_slice(slice,"");
 * @brief Draws cells.
 *
 * Draws the cells using the Toy Event Display object.
 *
 * @return None
  
display.draw_cells(toytiling.get_allcell(),"*same");
 * @brief Draws merge cells.
 *
 * Draws the merge cells using the Toy Event Display object.
 *
 * @return None
  
display.draw_mergecells(mergetiling.get_allcell(),"*same");
 * @brief Runs TApplication object.
 *
 * Runs the TApplication object.
 *
 * @return None
  
theApp.Run();
}
}
 * @brief Returns program exit status.
 *
 * Returns the program exit status as an integer.
 *
 * @return Integer value representing exit status
  
return 0;* This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.5.
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
