#include "Absorber.hh"

#include "CreateTree.hh"

#include <algorithm>
#include <string>
#include <sstream>

#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4PVParameterised.hh"
#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include <G4Cons.hh>

using namespace CLHEP;

Absorber::Absorber()
{
  // constructor
  // doesn't actually do anything, all will be set by the methods, called in DetectorConstruction

}

Absorber::~Absorber(){}


void Absorber::AddCell(Cell cell)
{
  cells.push_back(cell);
}

Cell Absorber::GetCell(int id)
{
  return cells[id];
}

G4int Absorber::GetNumberOfCells()
{
  return cells.size();
}

void Absorber::SetMaterial(G4int mat,G4double W_fraction)
{
  AbMaterial = NULL ;
  if      ( mat == 1 ) AbMaterial = MyMaterials::Brass () ;
  else if ( mat == 2 ) AbMaterial = MyMaterials::Tungsten () ;
  else if ( mat == 3 ) AbMaterial = MyMaterials::Lead () ;
  else if ( mat == 4 ) AbMaterial = MyMaterials::Iron () ;
  else if ( mat == 5 ) AbMaterial = MyMaterials::Aluminium () ;
  else if ( mat == 6 ) AbMaterial = MyMaterials::CopperTungstenAlloy(W_fraction) ;
  else
  {
    G4cerr << "<Absorber>: Invalid absorber material specifier " << mat << G4endl ;
    exit (-1) ;
  }
  // G4cout << "Absorber " << name << " material: "<< AbMaterial << G4endl ;
}
