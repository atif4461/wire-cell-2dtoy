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
 
 
 * @brief Prints usage message and exits if insufficient arguments are provided
 
 
 * @brief Creates a geometry data source from a file path
 *
 * @param argv Command line argument array containing file path
 * @return Geometry data source object
 
 
 * @brief Retrieves extent of geometry data source
 *
 * @return Vector of doubles representing extent
 
 
 * @brief Prints extent of geometry data source to error stream
 *
 * @param ex Vector of doubles representing extent
 
 
 * @brief Sets root file path from command line argument
 *
 * @param argv Command line argument array containing file path
 * @return Constant string pointer to root file path
 
 
 * @brief Creates a frame data source from a root file
 *
 * @param root_file Constant string pointer to root file path
 * @return Pointer to frame data source object
 
 
 * @brief Checks if frame data source creation was successful
 *
 * @param fds Pointer to frame data source object
 * @return Boolean indicating success or failure
 
 
 * @brief Creates a toy depositor object
 *
 * @param fds Pointer to frame data source object
 * @return Toy depositor object
 
 
 * @brief Retrieves depositions from toy depositor
 *
 * @param toydep Toy depositor object
 * @return Constant PointValueVector object
 
 
 * @brief Creates a generative frame data source object
 *
 * @param toydep Toy depositor object
 * @param gds Geometry data source object
 * @param parameters Integer and floating point values
 * @return Generative frame data source object
 
 
 * @brief Jumps to a specific index in generative frame data source
 *
 * @param gfds Generative frame data source object
 * @param index Integer value
 
 
 * @brief Creates a slice data source object
 *
 * @param gfds Generative frame data source object
 * @param parameter Integer value
 * @return Slice data source object
 
 
 * @brief Loops through slices in slice data source
 *
 * @param sds Slice data source object
 
 
 * @brief Jumps to a specific index in slice data source
 *
 * @param sds Slice data source object
 * @param i Integer value
 
 
 * @brief Retrieves a slice object
 *
 * @param sds Slice data source object
 * @return Slice object
 
 
 * @brief Checks if slice group has elements
 *
 * @param slice Slice object
 * @return Boolean indicating presence of elements
 
 
 * @brief Creates a toy tiling object
 *
 * @param slice Slice object
 * @param gds Geometry data source object
 * @return Toy tiling object
 
 
 * @brief Retrieves all cells from toy tiling object
 *
 * @param toytiling Toy tiling object
 * @return GeomCellSelection object
 
 
 * @brief Loops through cells in GeomCellSelection object
 *
 * @param allcell GeomCellSelection object
 
 
 * @brief Retrieves center point of cell
 *
 * @param allcell GeomCellSelection object
 * @param j Integer value
 * @return Point object
 
 
 * @brief Stores coordinates in arrays
 *
 * @param x Array of Doubles
 * @param y Array of Doubles
 * @param z Array of Doubles
 * @param ncount Integer value
 * @param i Integer value
 * @param p Point object
 
 
 * @brief Creates a TApplication object
 *
 * @param argc Number of command line arguments
 * @param argv Array of command line argument strings
 * @return TApplication object
 
 
 * @brief Sets return from run flag
 *
 * @param theApp TApplication object
 
 
 * @brief Creates a TCanvas object
 *
 * @param title String title
 * @param name String name
 * @param width Integer width
 * @param height Integer height
 * @return TCanvas object
 
 
 * @brief Draws TCanvas object
 *
 * @param c1 TCanvas object
 
 
 * @brief Sets option for statistics display
 *
 * @param gStyle Global style object
 
 
 * @brief Creates a TGraph2D object
 *
 * @param ncount Integer value
 * @param x Array of Doubles
 * @param y Array of Doubles
 * @param z Array of Doubles
 * @return TGraph2D object
 
 
 * @brief Draws TGraph2D object
 *
 * @param g TGraph2D object
 
 
 * @brief Opens a TFile object
 *
 * @param filename String filename
 * @param option String option
 * @return TFile object
 
 
 * @brief Writes TGraph2D object to TFile
 *
 * @param g TGraph2D object
 * @param name String name
 
 
 * @brief Writes TFile object
 *
 * @param file TFile object
 
 
 * @brief Closes TFile object
 *
 * @param file TFile object
 * This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.01.
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
  

  
  WCP::ToyDepositor toydep(fds);
  const PointValueVector pvv = toydep.depositions(1);
  

  WCP::GenerativeFDS gfds(toydep,gds,2400,5,2.0*1.6*units::millimeter);
  gfds.jump(1);

  WCPSst::ToyuBooNESliceDataSource sds(gfds,1500); //set threshold at 2000 electrons
  
  const int N = 100000;
  Double_t x[N],y[N],z[N];
  
  int ncount = 0;
  //int i=1140;{
  
  for (int i=0; i!=sds.size();i++){
    sds.jump(i);
    WCP::Slice slice = sds.get();
    if ( slice.group().size() >0){
      WCP2dToy::ToyTiling toytiling(slice,gds);
      //WCP2dToy::MergeToyTiling mergetiling(toytiling);
      
      GeomCellSelection allcell = toytiling.get_allcell();
      
      for (int j=0;j!=allcell.size();j++){
	Point p = allcell[j]->center();
	x[ncount] = i*0.32;
	y[ncount] = p.y/units::cm;
	z[ncount] = p.z/units::cm;
	ncount ++;
      }
      
      // GeomCellSelection allmcell = mergetiling.get_allcell();
      // GeomWireSelection allwire = mergetiling.get_allwire();

      // cout << i << " " << allmcell.size() << " " << allwire.size() << endl;

      // int sum = 0;
      // for (int j=0;j!=allmcell.size();j++){
      //   sum += ((WCP::MergeGeomCell*)allmcell[j])->get_allcell().size() ;
      // }
      // cout << allcell.size() << " " << allmcell.size() << " "  << sum << endl;
      
      
    }
  }


  TApplication theApp("theApp",&argc,argv);
  theApp.SetReturnFromRun(true);
  
  TCanvas c1("ToyMC","ToyMC",800,600);
  c1.Draw();
  gStyle->SetOptStat(0);

  TGraph2D *g = new TGraph2D(ncount,x,y,z);
  g->Draw("p");

  TFile *file = new TFile("shower3D.root","RECREATE");
  g->Write("shower3D");
  file->Write();
  file->Close();

  //theApp.Run();
  return 0;
  
} // main()
