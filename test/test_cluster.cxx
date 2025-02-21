#include "WCPData/Slim3DCluster.h"
#include "WCPSst/GeomDataSource.h"

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
  
int main(int argc, char* argv[]) 
 * @brief Checks command line arguments and prints usage message if invalid.
 *
 * Verifies that the required number of command line arguments is provided.
 * If not, it displays the correct usage and exits the program.
 
if (argc < 2) {
     * @brief Prints usage message to standard error stream.
   *
   * Displays the expected format for invoking the program.
   
  cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt" << std::endl;
  return 1;
}
 * @brief Creates a geometry data source object.
 *
 * Instantiates a WCPSst::GeomDataSource object using the first command line argument.
 
WCPSst::GeomDataSource gds(argv[1]);
 * @brief Retrieves the extent of the geometry data source.
 *
 * Obtains the spatial boundaries of the geometry data.
 *
 * @return Vector of doubles representing the extent.
 
std::vector<double> ex = gds.extent();
 * @brief Prints the extent of the geometry data source.
 *
 * Outputs the dimensions of the geometry to the standard error stream.
 
cerr << "Extent: "
     << " x:" << ex[0]/units::mm << " mm"
     << " y:" << ex[1]/units::m << " m"
     << " z:" << ex[2]/units::m << " m"
     << endl;
 * @brief Defines collections of wires in specific planes.
 *
 * Specifies sets of wires within particular planes for later reference.
 
const GeomWireSelection& uwires = gds.wires_in_plane(WirePlaneType_t(0));
const GeomWireSelection& vwires = gds.wires_in_plane(WirePlaneType_t(0));
const GeomWireSelection& wwires = gds.wires_in_plane(WirePlaneType_t(0));
 * @brief Constructs merge geometry cells with specified properties.
 *
 * Creates instances of SlimMergeGeomCell with distinct characteristics.
 
SlimMergeGeomCell *mcell1 = new SlimMergeGeomCell(1);
mcell1->SetTimeSlice(1);
 * @brief Populates the first merge cell with wires from various planes.
 *
 * Adds wires from different planes to the first merge cell.
 
for (int i=0;i!=100;i++){
  mcell1->AddWire(uwires.at(i),WirePlaneType_t(0));
  mcell1->AddWire(vwires.at(i),WirePlaneType_t(1));
  mcell1->AddWire(wwires.at(i),WirePlaneType_t(2));
}
 * @brief Creates additional merge cells and populates them with wires.
 *
 * Instantiates more SlimMergeGeomCell objects and adds wires to them.
 
SlimMergeGeomCell *mcell2 = new SlimMergeGeomCell(2);
mcell2->SetTimeSlice(2);
for (int i=25; i!=35;i++){
  mcell2->AddWire(uwires.at(i),WirePlaneType_t(0));
  mcell2->AddWire(vwires.at(i),WirePlaneType_t(1));
  mcell2->AddWire(wwires.at(i),WirePlaneType_t(2));
}
 * @brief Continues creating and populating merge cells.
 *
 * More SlimMergeGeomCell instances are created and populated with wires.
 
SlimMergeGeomCell *mcell3 = new SlimMergeGeomCell(2);
mcell3->SetTimeSlice(2);
for (int i=0; i!=10;i++){
  mcell3->AddWire(uwires.at(i),WirePlaneType_t(0));
  mcell3->AddWire(vwires.at(i),WirePlaneType_t(1));
  mcell3->AddWire(wwires.at(i),WirePlaneType_t(2));
}
 * @brief Creates another merge cell with specific settings.
 *
 * A new SlimMergeGeomCell is instantiated with certain attributes set.
 
SlimMergeGeomCell *mcell4 = new SlimMergeGeomCell(2);
mcell4->SetTimeSlice(3);
mcell4->add_bad_planes(WirePlaneType_t(2));
for (int i=8; i!=23;i++){
  mcell4->AddWire(uwires.at(i),WirePlaneType_t(0));
  mcell4->AddWire(vwires.at(i),WirePlaneType_t(1));
  mcell4->AddWire(wwires.at(i),WirePlaneType_t(2));
}
 * @brief Creates yet another merge cell similar to the previous one.
 *
 * Another instance of SlimMergeGeomCell is created with similar settings.
 
SlimMergeGeomCell *mcell4p = new SlimMergeGeomCell(2);
mcell4p->SetTimeSlice(3);
mcell4p->add_bad_planes(WirePlaneType_t(2));
for (int i=6; i!=25;i++){
  mcell4p->AddWire(uwires.at(i),WirePlaneType_t(0));
  mcell4p->AddWire(vwires.at(i),WirePlaneType_t(1));
  mcell4p->AddWire(wwires.at(i),WirePlaneType_t(2));
}
 * @brief Creates a final merge cell.
 *
 * The last SlimMergeGeomCell instance is created.
 
SlimMergeGeomCell *mcell5 = new SlimMergeGeomCell(2);
mcell5->SetTimeSlice(4);
for (int i=8; i!=23;i++){
  mcell5->AddWire(uwires.at(i),WirePlaneType_t(0));
  mcell5->AddWire(vwires.at(i),WirePlaneType_t(1));
  mcell5->AddWire(wwires.at(i),WirePlaneType_t(2));
}
 * @brief Assembles clusters from merge cells.
 *
 * Instances of Slim3DCluster are created and populated with merge cells.
 
Slim3DCluster *cluster = new Slim3DCluster(0,*mcell1);
cluster->AddCell(*mcell2);
cluster->AddCell(*mcell3);
cluster->AddCell(*mcell4);
cluster->AddCell(*mcell5);
cluster->Calc_Projection();
 * @brief Creates another cluster with slightly different composition.
 *
 * Another Slim3DCluster instance is created with some differences.
 
Slim3DCluster *cluster1 = new Slim3DCluster(1,*mcell1);
cluster1->AddCell(*mcell2);
// cluster1->AddCell(*mcell3);
cluster1->AddCell(*mcell4p);
cluster1->AddCell(*mcell5);
cluster1->Calc_Projection();
 * @brief Projects clusters onto planes.
 *
 * Projections of clusters onto specific planes are obtained.
 
Projected2DCluster *u_proj = cluster->get_projection(WirePlaneType_t(0));
Projected2DCluster *v_proj = cluster->get_projection(WirePlaneType_t(1));
Projected2DCluster *w_proj = cluster->get_projection(WirePlaneType_t(2));
 * @brief Projects the second cluster onto planes as well.
 *
 * Similar projections are obtained for the second cluster.
 
Projected2DCluster *u1_proj = cluster1->get_projection(WirePlaneType_t(0));
Projected2DCluster *v1_proj = cluster1->get_projection(WirePlaneType_t(1));
Projected2DCluster *w1_proj = cluster1->get_projection(WirePlaneType_t(2));
 * @brief Compares coverage between projected clusters.
 *
 * Evaluates how well one projection covers another.
 
std::cout << "Status: " << u_proj->judge_coverage(u1_proj) << std::endl;
std::cout << "Status: " << u1_proj->judge_coverage(u_proj) << std::endl;
 * @brief Calculates detailed comparison results between projections.
 *
 * Performs a more detailed analysis of the coverage relationship.
 
std::vector<int> comp_results = u_proj->calc_coverage(u1_proj);
std::cout << "Xin: " << " " << comp_results.at(0) << " " << comp_results.at(1) << " " <<
	comp_results.at(2) << " " << comp_results.at(3) << std::endl;* This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.5.
*/ 
int main(int argc, char* argv[])
{
  if (argc < 2) {
    cerr << "usage: wire-cell-uboone /path/to/ChannelWireGeometry.txt" << std::endl;
    return 1;
  }
  WCPSst::GeomDataSource gds(argv[1]);
  std::vector<double> ex = gds.extent();
  cerr << "Extent: "
       << " x:" << ex[0]/units::mm << " mm"
       << " y:" << ex[1]/units::m << " m"
       << " z:" << ex[2]/units::m << " m"
       << endl;

  const GeomWireSelection& uwires = gds.wires_in_plane(WirePlaneType_t(0));
  const GeomWireSelection& vwires = gds.wires_in_plane(WirePlaneType_t(0));
  const GeomWireSelection& wwires = gds.wires_in_plane(WirePlaneType_t(0));
  
  
  SlimMergeGeomCell *mcell1 = new SlimMergeGeomCell(1);
  mcell1->SetTimeSlice(1);
  for (int i=0;i!=100;i++){
    mcell1->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell1->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell1->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
  SlimMergeGeomCell *mcell2 = new SlimMergeGeomCell(2);
  mcell2->SetTimeSlice(2);
  for (int i=25; i!=35;i++){
    mcell2->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell2->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell2->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
  SlimMergeGeomCell *mcell3 = new SlimMergeGeomCell(2);
  mcell3->SetTimeSlice(2);
  for (int i=0; i!=10;i++){
    mcell3->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell3->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell3->AddWire(wwires.at(i),WirePlaneType_t(2));
  }

  SlimMergeGeomCell *mcell4 = new SlimMergeGeomCell(2);
  mcell4->SetTimeSlice(3);
  mcell4->add_bad_planes(WirePlaneType_t(2));
  for (int i=8; i!=23;i++){
    mcell4->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell4->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell4->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
  
  SlimMergeGeomCell *mcell4p = new SlimMergeGeomCell(2);
  mcell4p->SetTimeSlice(3);
  mcell4p->add_bad_planes(WirePlaneType_t(2));
  for (int i=6; i!=25;i++){
    mcell4p->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell4p->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell4p->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
  

  SlimMergeGeomCell *mcell5 = new SlimMergeGeomCell(2);
  mcell5->SetTimeSlice(4);
  for (int i=8; i!=23;i++){
    mcell5->AddWire(uwires.at(i),WirePlaneType_t(0));
    mcell5->AddWire(vwires.at(i),WirePlaneType_t(1));
    mcell5->AddWire(wwires.at(i),WirePlaneType_t(2));
  }
 
  
  Slim3DCluster *cluster = new Slim3DCluster(0,*mcell1);
  cluster->AddCell(*mcell2);
  cluster->AddCell(*mcell3);
  cluster->AddCell(*mcell4);
  cluster->AddCell(*mcell5);
  cluster->Calc_Projection();

  
  Slim3DCluster *cluster1 = new Slim3DCluster(1,*mcell1);
  cluster1->AddCell(*mcell2);
  // cluster1->AddCell(*mcell3);
  cluster1->AddCell(*mcell4p);
  cluster1->AddCell(*mcell5);
  cluster1->Calc_Projection();

  
  Projected2DCluster *u_proj = cluster->get_projection(WirePlaneType_t(0));
  Projected2DCluster *v_proj = cluster->get_projection(WirePlaneType_t(1));
  Projected2DCluster *w_proj = cluster->get_projection(WirePlaneType_t(2));

  Projected2DCluster *u1_proj = cluster1->get_projection(WirePlaneType_t(0));
  Projected2DCluster *v1_proj = cluster1->get_projection(WirePlaneType_t(1));
  Projected2DCluster *w1_proj = cluster1->get_projection(WirePlaneType_t(2));
  
  std::cout << "Status: " << u_proj->judge_coverage(u1_proj) << std::endl;
  std::cout << "Status: " << u1_proj->judge_coverage(u_proj) << std::endl;

  std::vector<int> comp_results = u_proj->calc_coverage(u1_proj);
  std::cout << "Xin: " << " " << comp_results.at(0) << " " << comp_results.at(1) << " " <<
	  comp_results.at(2) << " " << comp_results.at(3) << std::endl;
	

  
}
