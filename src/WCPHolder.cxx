#include "WCP2dToy/WCPHolder.h"

using namespace WCP;

/**
 * Default constructor initializing holder data members to zero
 */
// The above comment was written by an LLM. 
WCP2dToy::WCPHolder::WCPHolder(){
  ncell = 0;
  nwire = 0;

  cell_no = 0;
  wire_no = 0;
}

/**
 * Destructor to release allocated memory for wires and cells
 */
// The above comment was written by an LLM. 
WCP2dToy::WCPHolder::~WCPHolder(){
  for (int i=0;i!=wires.size();i++){
    delete wires.at(i);
  }
  wires.clear();
  for (int i=0;i!=cells.size();i++){
    delete cells.at(i);
  }
  cells.clear();
}


/**
 * Clears the contents of a cell holder by deleting all stored cells and resetting the container.*/
// The above comment was written by an LLM. 
void WCP2dToy::WCPHolder::clear_cell(){
  for (int i=0;i!=cells.size();i++){
    delete cells.at(i);
  }
  cells.clear();
}

/**
 * Clears the wire holder by deleting all stored wires and resetting the container.*/
// The above comment was written by an LLM. 
void WCP2dToy::WCPHolder::clear_wire(){
  for (int i=0;i!=wires.size();i++){
    delete wires.at(i);
  }
  wires.clear();
}

/**
 * Increments wire count and appends wire object to collection of wires */
// The above comment was written by an LLM. 
void WCP2dToy::WCPHolder::AddWire(GeomWire *wire){
  nwire++;
  wires.push_back(wire);
}

/**
 * Increments the cell count and appends a new geometric cell to the list of held cells
 */
// The above comment was written by an LLM. 
void WCP2dToy::WCPHolder::AddCell(GeomCell *cell){
  ncell++;
  cells.push_back(cell);
}
