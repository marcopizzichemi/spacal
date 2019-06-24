#ifndef Absorber_h
#define Absorber_h 1

#include <iostream>
#include <string>
#include <fstream>
#include <utility>

#include "ConfigFile.hh"
#include "MyMaterials.hh"
#include "LedFiberTiming.hh"
#include "DetectorParameterisation.hh"

#include "G4Material.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4VisAttributes.hh"
#include "G4UniformMagField.hh"

#include "Cell.hh"


// class for Module
class Absorber
{
public:

  //! ctor
  Absorber ();

  //! dtor
  ~Absorber ();

  //! construct method
  // G4VPhysicalVolume* Construct () ;
  //
  // //! other methods
  void  SetID(int aNum){id = aNum;};
  G4int GetID(){return id;}
  //
  void     SetName(G4String aString){name = aString;};
  G4String GetName(){return name;};

  void AddCell(Cell cell);
  Cell GetCell(int id);
  G4int GetNumberOfCells();
  //
  void     SetDimensions(G4double x, G4double y ,G4double z){size_x = x;size_y = y; size_z = z;};
  G4double GetSizeX(){return size_x;};
  G4double GetSizeY(){return size_y;};
  G4double GetSizeZ(){return size_z;};
  //
  void       SetPosition(G4double x, G4double y ,G4double z){pos_x = x;pos_y = y; pos_z = z;};
  G4double   GetPositionX(){return pos_x;};
  G4double   GetPositionY(){return pos_y;};
  G4double   GetPositionZ(){return pos_z;};
  //
  void          SetMaterial(G4int aMat,G4double W_fraction);
  G4Material*   GetMaterial(){return AbMaterial;};

  void     SetEsrThickness(G4double x){esr_thickness = x;};
  G4double GetEsrThickness(){return esr_thickness;};



private:

  //-------------------------------------------//
  // identifiers                               //
  //-------------------------------------------//
  G4int id;
  G4String name;

  std::vector<Cell> cells;

  //-------------------------------------------//
  // dimensions, material etc                  //
  //-------------------------------------------//
  // dimensions of the cell == dimension of the bulk absorber material
  G4double size_x;             // abs dimensions in x [mm]
  G4double size_y;             // abs dimensions in y [mm]
  G4double size_z;             // abs dimensions in z [mm]
  // // position in mother volume
  G4double pos_x;             // position in x [mm]
  G4double pos_y;             // position in y [mm]
  G4double pos_z;             // position in z [mm]
  // material
  G4Material* AbMaterial;
  G4double esr_thickness;



} ;

#endif /*Absorber_h*/
