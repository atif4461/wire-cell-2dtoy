#ifndef WIRECELL2DTOY_TOYEVENTDISPLAY_H
#define WIRECELL2DTOY_TOYEVENTDISPLAY_H
#include "WCPData/Point.h"
#include "WCPData/GeomCell.h"
#include "WCPData/MergeGeomCell.h"


#include "WCP2dToy/ToyMatrix.h"

#include "WCPNav/SliceDataSource.h"
#include "WCPNav/GeomDataSource.h"
#include "WCPNav/DetectorGDS.h"

#include "TPad.h"
#include "TH2F.h"
#include "TString.h"
#include "TGraph.h"
#include "TPolyLine.h"

namespace WCP2dToy {

/**
 * @brief Class representing a toy event display.
 *
 * This class is responsible for creating a graphical representation of an event,
 * including drawing cell charges, slices, wires, and points. It also provides
 * functionality for setting thresholds and clearing the display data.
 
  * @brief Constructor for the ToyEventDisplay class.
 *
 * Creates a new instance of the ToyEventDisplay class, initializing it with the
 * given TPad and GeomDataSource or DetectorGDS.
 *
 * @param pad The TPad to draw into.
 * @param gds The GeomDataSource or DetectorGDS to use.
 
  * @brief Destructor for the ToyEventDisplay class.
 *
 * Releases any resources held by the ToyEventDisplay instance.
 
  * @brief Initializes the display with the given extents.
 *
 * Sets up the display with the specified minimum and maximum x and y values.
 *
 * @param x_min The minimum x value.
 * @param x_max The maximum x value.
 * @param y_min The minimum y value.
 * @param y_max The maximum y value.
 * @return An integer indicating success or failure.
 
  * @brief Draws the cell charges.
 *
 * Displays the cell charges on the screen, using the given flag and cell charges.
 *
 * @param flag A flag controlling the drawing behavior.
 * @param cellcharges The vector of cell charges to draw.
 * @param option A string specifying additional drawing options.
 * @return An integer indicating success or failure.
 
  * @brief Draws a slice.
 *
 * Displays the given slice on the screen, using the specified option.
 *
 * @param slice The slice to draw.
 * @param option A string specifying additional drawing options.
 * @return An integer indicating success or failure.
 
  * @brief Draws a set of wires.
 *
 * Displays the given wires on the screen, using the specified option.
 *
 * @param wires The wires to draw.
 * @param option A string specifying additional drawing options.
 * @return An integer indicating success or failure.
 
  * @brief Draws a merged set of wires.
 *
 * Displays the given merged wires on the screen, using the specified option and color.
 *
 * @param wires The merged wires to draw.
 * @param option A string specifying additional drawing options.
 * @param color The color to use for drawing.
 * @return An integer indicating success or failure.
 
  * @brief Draws a set of points.
 *
 * Displays the given points on the screen, using the specified option and color.
 *
 * @param pcells The points to draw.
 * @param option A string specifying additional drawing options.
 * @param color The color to use for drawing.
 * @return An integer indicating success or failure.
 
  * @brief Draws a selection of cells.
 *
 * Displays the given cells on the screen, using the specified option and color.
 *
 * @param cellall The cells to draw.
 * @param option A string specifying additional drawing options.
 * @param color The color to use for drawing.
 * @return An integer indicating success or failure.
 
  * @brief Draws a merged selection of cells.
 *
 * Displays the given merged cells on the screen, using the specified option and flag.
 *
 * @param cellall The merged cells to draw.
 * @param option A string specifying additional drawing options.
 * @param flag A flag controlling the drawing behavior.
 * @return An integer indicating success or failure.
 
  * @brief Draws truth cells.
 *
 * Displays the given truth cells on the screen, using the specified option.
 *
 * @param ccmap The cell charge map to use.
 * @param option A string specifying additional drawing options.
 * @return An integer indicating success or failure.
 
  * @brief Draws truth cells with charge information.
 *
 * Displays the given truth cells with charge information on the screen, using the specified option and FI.
 *
 * @param ccmap The cell charge map to use.
 * @param option A string specifying additional drawing options.
 * @param FI Additional information to use for drawing.
 * @return An integer indicating success or failure.
 
  * @brief Draws wire charges.
 *
 * Displays the given wire charges on the screen, using the specified option and FI.
 *
 * @param wcmap The wire charge map to use.
 * @param option A string specifying additional drawing options.
 * @param FI Additional information to use for drawing.
 * @return An integer indicating success or failure.
 
  * @brief Draws cell charges.
 *
 * Displays the given cell charges on the screen, using the specified option.
 *
 * @param cellall The cells to draw.
 * @param option A string specifying additional drawing options.
 * @return An integer indicating success or failure.
 
  * @brief Draws reconstructed cells.
 *
 * Displays the given reconstructed cells on the screen, using the specified option and color.
 *
 * @param cellall The reconstructed cells to draw.
 * @param toymatrix The toy matrix to use.
 * @param option A string specifying additional drawing options.
 * @param color The color to use for drawing.
 * @return An integer indicating success or failure.
 
  * @brief Sets the reconstruction threshold.
 *
 * Sets the threshold value for reconstructing cells.
 *
 * @param abc The new threshold value.
 
  * @brief Sets the truth threshold.
 *
 * Sets the threshold value for truth cells.
 *
 * @param abc The new threshold value.
 
  * @brief Draws a bad region.
 *
 * Displays the given bad region on the screen, using the specified time, scale, plane, and option.
 *
 * @param chirpmap The chirp map to use.
 * @param time The time to use for drawing.
 * @param scale The scale to use for drawing.
 * @param plane The plane to use for drawing.
 * @param option A string specifying additional drawing options.
 
  * @brief Draws a bad cell.
 *
 * Displays the given bad cell on the screen.
 *
 * @param cells The bad cells to draw.
 
  * @brief Clears the event display data.
 *
 * Removes all data from the display.
 * This comment was generated by meta-llama/Llama-3.3-70B-Instruct:None at temperature 0.01.
*/ 
    class ToyEventDisplay {
    private:
    
	TPad &pad;
	int gds_flag;
	const WCP::GeomDataSource* gds;
	const WCP::DetectorGDS* dgds;

	TH2F *h1;
	TH2F *h2;
	TGraph *g1;
	TGraph *g2;
	TGraph *g2b;
	TPolyLine *g3;
    
    public:
	/// Create a ToyEventDisplay drawing into the given TPad using the given GeomDataSource.
	ToyEventDisplay(TPad& pad, const WCP::GeomDataSource& gds);
	ToyEventDisplay(TPad& pad, const WCP::DetectorGDS& gds);
	virtual ~ToyEventDisplay();
    
	/// Initialize with extents
	virtual int init(float x_min=4.9, float x_max=6.1, float y_min=-1.1, float y_max=1.1);
    
	/// Draw the cell charges
	virtual int draw_mc(int flag, const WCP::PointValueVector& cellcharges, TString option);
    
	/// Draw a slice
	virtual int draw_slice(const WCP::Slice& slice, TString option);
	virtual int draw_wires(WCP::GeomWireSelection& wires, TString option);
	virtual int draw_merged_wires(WCP::GeomWireSelection wires, TString option, int color = 2);

	virtual int draw_points(WCP::PointVector pcells, TString option, int color = 2);

	/// Draw a selection of cells
	virtual int draw_cells(const WCP::GeomCellSelection& cellall ,TString option,int color=4);
	virtual int draw_mergecells(const WCP::GeomCellSelection& cellall ,TString option, int flag=0);
	virtual int draw_truthcells(const WCP::CellChargeMap& ccmap,TString option);
	virtual int draw_truthcells_charge(const WCP::CellChargeMap& ccmap,TString option, int FI);
	virtual int draw_wires_charge(const WCP::WireChargeMap& wcmap,TString option, int FI);
	virtual int draw_cells_charge(const WCP::GeomCellSelection& cellall ,TString option);

	virtual int draw_reconcells(const WCP::GeomCellSelection& cellall, WCP2dToy::ToyMatrix *toymatrix ,TString option, int color = 4);

	void Set_ReconThreshold(int abc){recon_threshold = abc;};
	void Set_TruthThreshold(int abc){truth_threshold = abc;};

	void draw_bad_region(WCP::ChirpMap& chirpmap, int time, int scale, int plane, TString option);
	void draw_bad_cell(WCP::GeomCellSelection& cells);

	/// Clear visual the event display data.
	void clear();

	Double_t charge_min, charge_max;

    private:
	int recon_threshold;
	int truth_threshold;
    };

}

#endif
