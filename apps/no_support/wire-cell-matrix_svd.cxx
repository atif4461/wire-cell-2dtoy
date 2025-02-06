#include "WCPSst/GeomDataSource.h"
//#include "WCPSst/ToyuBooNEFrameDataSource.h"
#include "WCPSst/ToyuBooNESliceDataSource.h"
#include "WCP2dToy/ToyEventDisplay.h"
#include "WCP2dToy/ToyTiling.h"
#include "WCP2dToy/MergeToyTiling.h"
#include "WCP2dToy/TruthToyTiling.h"
#include "WCP2dToy/ToyMatrix.h"
#include "WCP2dToy/ToyMatrixExclusive.h"
#include "WCP2dToy/ToyMatrixKalman.h"
#include "WCP2dToy/ToyMatrixIterate.h"
#include "WCP2dToy/ToyMetric.h"

#include "WCPData/MergeGeomCell.h"
#include "WCPData/MergeGeomWire.h"

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
 
if (argc < 3) {
   * @brief Print usage message
 
cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt /path/to/celltree.root" << endl;
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
 * @brief Set root file path and tree path
 
const char* root_file = argv[2];
const char* tpath = "/Event/Sim";
 * @brief Create frame data source object
 *
 * @param root_file Root file path
 * @return Frame data source pointer
 
WCP::FrameDataSource* fds = 0;
fds = WCPSst::make_fds(root_file);
 * @brief Check if frame data source creation was successful
 *
 * @return Error message if creation failed
 
if (!fds) {
cerr << "ERROR: failed to get FDS from " << root_file << endl;
return 1;
}
 * @brief Create toy depositor object
 *
 * @param fds Frame data source pointer
 
WCP::ToyDepositor toydep(fds);
 * @brief Get depositions from toy depositor
 *
 * @param index Index value
 * @return Point value vector
 
const PointValueVector pvv = toydep.depositions(1);
 * @brief Create generative frame data source object
 *
 * @param toydep Toy depositor object
 * @param gds Geometry data source object
 * @param args Additional arguments
 
WCP::GenerativeFDS gfds(toydep,gds,2400,5,2.0*1.6*units::millimeter);
 * @brief Jump to specified index in generative frame data source
 *
 * @param index Index value
 
gfds.jump(1);
 * @brief Create slice data source object
 *
 * @param gfds Generative frame data source object
 * @param arg Argument value
 
WCPSst::ToyuBooNESliceDataSource sds(gfds,1500);
 * @brief Initialize counters
 
int ncount = 0;
int ncount1 = 0;
int ncount_t = 0;
 * @brief Allocate memory for tiling objects
 
WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];
WCP2dToy::ToyMatrix **toymatrix = new WCP2dToy::ToyMatrix*[2400];
WCP2dToy::ToyMatrixIterate **toymatrix_it = new WCP2dToy::ToyMatrixIterate*[2400];
 * @brief Create metric object
 
WCP2dToy::ToyMetric toymetric;
 * @brief Loop through slices
 *
 * @param start_num Starting index
 * @param end_num Ending index
 
for (int i=start_num;i!=end_num+1;i++){
 * @brief Jump to current index in slice data source
 
sds.jump(i);
 * @brief Get slice object
 *
 * @return Slice object
 
WCP::Slice slice = sds.get();
 * @brief Create tiling objects
 *
 * @param slice Slice object
 * @param gds Geometry data source object
 
toytiling[i] = new WCP2dToy::ToyTiling(slice,gds);
mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i);
 * @brief Get cell and wire selections
 *
 * @return Cell and wire selection objects
 
GeomCellSelection allcell = toytiling[i]->get_allcell();
GeomWireSelection allwire = toytiling[i]->get_allwire();
GeomCellSelection allmcell = mergetiling[i]->get_allcell();
GeomWireSelection allmwire = mergetiling[i]->get_allwire();
 * @brief Print cell and wire counts
 
cout << i << " " << allmcell.size() << " " << allmwire.size() << endl;
 * @brief Create truth tiling object
 *
 * @param toytiling Tiling object
 * @param pvv Point value vector
 * @param index Index value
 * @param gds Geometry data source object
 
truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds);
 * @brief Create matrix object
 *
 * @param toytiling Tiling object
 * @param mergetiling Merge tiling object
 * @param args Additional arguments
 
toymatrix[i] = new WCP2dToy::ToyMatrix(*toytiling[i],*mergetiling[i],1,2000);
 * @brief Check solve flag and create iterate object if necessary
 
if (toymatrix[i]->Get_Solve_Flag()==0)
toymatrix_it[i] = new WCP2dToy::ToyMatrixIterate(*toymatrix[i],toymatrix[i]->Get_svd_removed());
 * @brief Print chi2 and NDF values
 
cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
 * @brief Add cells to metric object
 *
 * @param allmcell Cell selection object
 * @param toymatrix Matrix object
 * @param ccmap Cell charge map
 
CellChargeMap ccmap = truthtiling[i]->ccmap();
if (toymatrix[i]->Get_Solve_Flag()!=0)
toymetric.Add(allmcell,*toymatrix[i],ccmap);
 * @brief Add solve flag to metric object
 *
 * @param toymatrix Matrix object
 
toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());

}
 * @brief Print metric object
 
toymetric.Print();
 * @brief Open output file
 *
 * @param filename File name
 * @param mode Mode string
 * @return File pointer
 
TFile *file = new TFile("shower3D.root","RECREATE");
 * @brief Create trees
 *
 * @param file File pointer
 * @param treename Tree name
 * @param title Title string
 * @return Tree pointer
 
TTree *t_true = new TTree("T_true","T_true");
TTree *t_rec = new TTree("T_rec","T_rec");
TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
 * @brief Declare variables for branch creation
 
Double_t x_save, y_save, z_save;
Double_t charge_save;
Double_t chi2_save;
Double_t ndf_save;
Double_t ncell_save;
 * @brief Set directory for trees
 *
 * @param file File pointer
 
t_true->SetDirectory(file);
t_rec->SetDirectory(file);
t_rec_charge->SetDirectory(file);
 * @brief Create branches
 *
 * @param tree Tree pointer
 * @param var Variable address
 * @param varname Variable name
 * @param vartype Variable type
 
t_true->Branch("x",&x_save,"x/D");
t_true->Branch("y",&y_save,"y/D");
t_true->Branch("z",&z_save,"z/D");
t_true->Branch("q",&charge_save,"q/D");
t_rec->Branch("x",&x_save,"x/D");
t_rec->Branch("y",&y_save,"y/D");
t_rec->Branch("z",&z_save,"z/D");
t_rec_charge->Branch("x",&x_save,"x/D");
t_rec_charge->Branch("y",&y_save,"y/D");
t_rec_charge->Branch("z",&z_save,"z/D");
t_rec_charge->Branch("q",&charge_save,"q/D");
t_rec_charge->Branch("ncell",&ncell_save,"ncell/D");
t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");
 * @brief Create graphs
 *
 * @return Graph pointer
 
TGraph2D *g = new TGraph2D();
TGraph2D *gt = new TGraph2D();
TGraph2D *g_rec = new TGraph2D();
 * @brief Loop through slices again
 *
 * @param start_num Starting index
 * @param end_num Ending index
 
for (int i=start_num;i!=end_num+1;i++){
 * @brief Get cell charge map
 *
 * @param truthtiling Truth tiling object
 * @return Cell charge map
 
CellChargeMap ccmap = truthtiling[i]->ccmap();
 * @brief Loop through cell charge map
 *
 * @param it Iterator
 
for (auto it = ccmap.begin();it!=ccmap.end(); it++){
Point p = it->first->center();
x_save = i*0.32;
y_save = p.y/units::cm;
z_save = p.z/units::cm;
charge_save = it->second;
 * @brief Set point in graph
 *
 * @param graph Graph pointer
 * @param ncount Counter variable
 * @param x X-coordinate
 * @param y Y-coordinate
 * @param z Z-coordinate
 
gt->SetPoint(ncount_t,x_save,y_save,z_save);
t_true->Fill();
 * @brief Increment counter
 
ncount_t ++;
}
 * @brief Get cell selection
 *
 * @param toytiling Tiling object
 * @return Cell selection
 
GeomCellSelection allcell = toytiling[i]->get_allcell();
 * @brief Loop through cell selection
 *
 * @param j Index variable
 
for (int j=0;j!=allcell.size();j++){
Point p = allcell[j]->center();
x_save = i*0.32;
y_save = p.y/units::cm;
z_save = p.z/units::cm;
 * @brief Set point in graph
 *
 * @param graph Graph pointer
 * @param ncount Counter variable
 * @param x X-coordinate
 * @param y Y-coordinate
 * @param z Z-coordinate
 
g->SetPoint(ncount,x_save,y_save,z_save);
t_rec->Fill();
 * @brief Increment counter
 
ncount ++;
}
 * @brief Get merge cell selection
 *
 * @param mergetiling Merge tiling object
 * @return Merge cell selection
 
GeomCellSelection allmcell = mergetiling[i]->get_allcell();
 * @brief Loop through merge cell selection
 *
 * @param j Index variable
 
for (int j=0;j!=allmcell.size();j++){
MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
 * @brief Check charge condition
 
if (charge>2000){
 * @brief Loop through sub-cells
 *
 * @param k Index variable
 
for (int k=0;k!=mcell->get_allcell().size();k++){
Point p = mcell->get_allcell().at(k)->center();
x_save = i*0.32;
y_save = p.y/units::cm;
z_save = p.z/units::cm;
charge_save = charge/mcell->get_allcell().size();
ncell_save = mcell->get_allcell().size();
chi2_save = toymatrix[i]->Get_Chi2();
ndf_save = toymatrix[i]->Get_ndf();
 * @brief Set point in graph
 *
 * @param graph Graph pointer
 * @param ncount Counter variable
 * @param x X-coordinate
 * @param y Y-coordinate
 * @param z Z-coordinate
 
g_rec->SetPoint(ncount1,x_save,y_save,z_save);
t_rec_charge->Fill();
 * @brief Increment counter
 
ncount1 ++;
}
}
}

}
 * @brief Write graphs to file
 *
 * @param graph Graph pointer
 * @param filename File name
 
g->Write("shower3D");
gt->Write("shower3D_truth");
g_rec->Write("shower3D_charge");
 * @brief Close file
 
file->Write();
file->Close();
 * @brief Print metric object
 
toymetric.Print();
 * @brief Return program exit status
 
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

  

  
  int ncount = 0;
  int ncount1 = 0;
  int ncount_t = 0;
  

  WCP2dToy::ToyTiling **toytiling = new WCP2dToy::ToyTiling*[2400];
  WCP2dToy::MergeToyTiling **mergetiling = new WCP2dToy::MergeToyTiling*[2400];
  WCP2dToy::TruthToyTiling **truthtiling = new WCP2dToy::TruthToyTiling*[2400];
  WCP2dToy::ToyMatrix **toymatrix = new WCP2dToy::ToyMatrix*[2400];
  WCP2dToy::ToyMatrixIterate **toymatrix_it = new WCP2dToy::ToyMatrixIterate*[2400];
  
  WCP2dToy::ToyMetric toymetric;


  //add in cluster
  GeomClusterSet cluster_set, cluster_delset;
  
  int ncount_mcell = 0;

  // int start_num = 0 ;
  // int end_num = sds.size()-1;

  int start_num =0;
  int end_num = sds.size()-1;

  // int start_num = 462;
  // int end_num = 465;

 
  for (int i=start_num;i!=end_num+1;i++){
 
    sds.jump(i);
    WCP::Slice slice = sds.get();
    //if ( slice.group().size() >0){
      
      toytiling[i] = new WCP2dToy::ToyTiling(slice,gds);
      mergetiling[i] = new WCP2dToy::MergeToyTiling(*toytiling[i],i);

      GeomCellSelection allcell = toytiling[i]->get_allcell();
      GeomWireSelection allwire = toytiling[i]->get_allwire();
      GeomCellSelection allmcell = mergetiling[i]->get_allcell();
      GeomWireSelection allmwire = mergetiling[i]->get_allwire();
      
      cout << i << " " << allmcell.size() << " " << allmwire.size() << endl;

      truthtiling[i] = new WCP2dToy::TruthToyTiling(*toytiling[i],pvv,i,gds);
      toymatrix[i] = new WCP2dToy::ToyMatrix(*toytiling[i],*mergetiling[i],1,2000);
      if (toymatrix[i]->Get_Solve_Flag()==0)
	toymatrix_it[i] = new WCP2dToy::ToyMatrixIterate(*toymatrix[i],toymatrix[i]->Get_svd_removed());
      
      
      cout << "chi2: " << toymatrix[i]->Get_Chi2() << endl;
      cout << "NDF: " << toymatrix[i]->Get_ndf() << endl;
      
      


      GeomCellSelection calmcell;
      for (int j=0;j!=allmcell.size();j++){
	MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
	double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
	double charge_err = toymatrix[i]->Get_Cell_Charge(mcell,2);
	
	//	cout << "Recon: " << j << " " << charge << " " << charge_err << endl;

	if (charge > 2000) calmcell.push_back(mcell);
      }
      

      


      CellChargeMap ccmap = truthtiling[i]->ccmap();
      if (toymatrix[i]->Get_Solve_Flag()!=0)
	toymetric.Add(allmcell,*toymatrix[i],ccmap);

      toymetric.AddSolve(toymatrix[i]->Get_Solve_Flag());

      Double_t charge_min = 10000;
      Double_t charge_max = 0;

     
      
      
   
    
  }

  toymetric.Print();
  


  TFile *file = new TFile("shower3D.root","RECREATE");
  TTree *t_true = new TTree("T_true","T_true");
  TTree *t_rec = new TTree("T_rec","T_rec");
  TTree *t_rec_charge = new TTree("T_rec_charge","T_rec_charge");
  
  Double_t x_save, y_save, z_save;
  Double_t charge_save;
  Double_t chi2_save;
  Double_t ndf_save;
  Double_t ncell_save;

  t_true->SetDirectory(file);
  t_true->Branch("x",&x_save,"x/D");
  t_true->Branch("y",&y_save,"y/D");
  t_true->Branch("z",&z_save,"z/D");
  t_true->Branch("q",&charge_save,"q/D");
  
  t_rec->SetDirectory(file);
  t_rec->Branch("x",&x_save,"x/D");
  t_rec->Branch("y",&y_save,"y/D");
  t_rec->Branch("z",&z_save,"z/D");
  
  t_rec_charge->SetDirectory(file);
  t_rec_charge->Branch("x",&x_save,"x/D");
  t_rec_charge->Branch("y",&y_save,"y/D");
  t_rec_charge->Branch("z",&z_save,"z/D");
  t_rec_charge->Branch("q",&charge_save,"q/D");
  t_rec_charge->Branch("ncell",&ncell_save,"ncell/D");
  t_rec_charge->Branch("chi2",&chi2_save,"chi2/D");
  t_rec_charge->Branch("ndf",&ndf_save,"ndf/D");

  TGraph2D *g = new TGraph2D();
  TGraph2D *gt = new TGraph2D();
  TGraph2D *g_rec = new TGraph2D();


  //save results 
  for (int i=start_num;i!=end_num+1;i++){
    //truth
    CellChargeMap ccmap = truthtiling[i]->ccmap();
    for (auto it = ccmap.begin();it!=ccmap.end(); it++){
      Point p = it->first->center();
      x_save = i*0.32;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      charge_save = it->second;
      
      gt->SetPoint(ncount_t,x_save,y_save,z_save);
      t_true->Fill();
            
      ncount_t ++;
    }
    
    //recon 1
    GeomCellSelection allcell = toytiling[i]->get_allcell();
    for (int j=0;j!=allcell.size();j++){
      Point p = allcell[j]->center();
      x_save = i*0.32;
      y_save = p.y/units::cm;
      z_save = p.z/units::cm;
      

      g->SetPoint(ncount,x_save,y_save,z_save);
      t_rec->Fill();

      ncount ++;
    }

    //recon 2 with charge
    GeomCellSelection allmcell = mergetiling[i]->get_allcell();
    for (int j=0;j!=allmcell.size();j++){
      MergeGeomCell *mcell = (MergeGeomCell*)allmcell[j];
      double charge = toymatrix[i]->Get_Cell_Charge(mcell,1);
      if (charge>2000){
	//truth
	for (int k=0;k!=mcell->get_allcell().size();k++){
	  Point p = mcell->get_allcell().at(k)->center();
	  x_save = i*0.32;
	  y_save = p.y/units::cm;
	  z_save = p.z/units::cm;
	  charge_save = charge/mcell->get_allcell().size();
	  ncell_save = mcell->get_allcell().size();
	  chi2_save = toymatrix[i]->Get_Chi2();
	  ndf_save = toymatrix[i]->Get_ndf();

	  g_rec->SetPoint(ncount1,x_save,y_save,z_save);
	  t_rec_charge->Fill();
	  
	  ncount1 ++;
	}
      }
    }

  }
 
  

  g->Write("shower3D");
  gt->Write("shower3D_truth");
  g_rec->Write("shower3D_charge");
  file->Write();
  file->Close();

  toymetric.Print();

  return 0;
  
} // main()
