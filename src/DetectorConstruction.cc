//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes, nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#include "DetectorConstruction.hh"

#include "DetectorParameterisation.hh"
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

#include <Cell.hh>
#include <Absorber.hh>

using namespace CLHEP;



DetectorConstruction::DetectorConstruction (const string& configFileName)
{
  //-----------------------------------------------------
  //------------- Define colors --------------
  //-----------------------------------------------------
  white   = G4Colour(1.00, 1.00, 1.00) ;  // white
  grey    = G4Colour(0.50, 0.50, 0.50) ;  // grey
  black   = G4Colour(0.00, 0.00, 0.00) ;  // black
  red     = G4Colour(1.00, 0.00, 0.00) ;  // red
  green   = G4Colour(0.00, 1.00, 0.00) ;  // green
  blue    = G4Colour(0.00, 0.00, 1.00) ;  // blue
  cyan    = G4Colour(0.00, 1.00, 1.00) ;  // cyan
  air     = G4Colour(0.00, 1.00, 1.00) ;  // cyan
  magenta = G4Colour(1.00, 0.00, 1.00) ;  // magenta
  yellow  = G4Colour(1.00, 1.00, 0.00) ;  // yellow
  brass   = G4Colour(0.80, 0.60, 0.40) ;  // brass
  brown   = G4Colour(0.70, 0.40, 0.10) ;  // brown
  orange  = G4Colour(1.00, 0.33, 0.00) ;  // orange



  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------

  ConfigFile config (configFileName) ;
  config.readInto (checkOverlaps, "checkOverlaps") ;
  config.readInto (world_material, "world_material") ;
  config.readInto (W_fraction, "W_fraction") ;
  config.readInto (surface_lg,"surface_lg");
  config.readInto (glue_interface,"glue_interface");
  config.readInto (cone_material,"cone_material");
  B_field_intensity = config.read<double>("B_field_intensity") * tesla ;

  //---------------------------------------//
  // CALORIMETER                           //
  //---------------------------------------//
  config.readInto(modules_nx,"modules_nx");
  config.readInto(modules_ny,"modules_ny");
  //visibility
  config.readInto(visibility,"crystalsVisibility");
  config.readInto(wireFrame,"wireFrame");

  //---------------------------------------//
  // ABSORBER                              //
  //---------------------------------------//
  // read from config fil
  config.readInto(AbsName ,    "absorber_name");
  config.readInto(AbsSizeX,    "absorber_size_x");
  config.readInto(AbsSizeY,    "absorber_size_y");
  config.readInto(AbsSizeZ,    "absorber_size_z");
  config.readInto(AbsPositionX,"absorber_pos_x");
  config.readInto(AbsPositionY,"absorber_pos_y");
  config.readInto(AbsPositionZ,"absorber_pos_z");
  config.readInto(AbsMaterial, "absorber_material");

  cell_separation_type = config.read<int>("cell_separation_type",0);
  cell_separator_position = config.read<int>("cell_separator_position",0);
  esr_thickness = config.read<double>("esr_thickness",0);
  if(cell_separation_type != 2)
  {
    esr_thickness = 0;
  }

  // rescale AbsSizeZ, adding esr_thickness if it's there
  AbsSizeZ = AbsSizeZ + esr_thickness;
  absorber.SetName(AbsName);
  absorber.SetDimensions( AbsSizeX,
                          AbsSizeY,
                          AbsSizeZ);
  absorber.SetPosition(AbsPositionX,
                            AbsPositionY,
                            AbsPositionZ);
  absorber.SetMaterial(AbsMaterial,W_fraction);
  absorber.SetEsrThickness(esr_thickness);


  //---------------------------------------//
  // CELLS                                 //
  //---------------------------------------//
  // read from config file
  std::vector<G4String>   CellNames;
  std::vector<G4double>   CellPositionX;
  std::vector<G4double>   CellPositionY;
  std::vector<G4double>   CellPositionZ;
  std::vector<G4int>      CellMaterial;
  std::vector<G4int>      CellXelements;
  std::vector<G4int>      CellYelements;
  std::vector<G4double>   CellAirLayer  ;
  std::vector<G4int>      CellExtGapMaterial;
  std::vector<G4int>      CellIntGapMaterial;
  std::vector<G4int>      CellStaggering;
  std::vector<G4int>      CellStaggeringAxis;
  std::vector<G4double>   CellStaggeringSize  ;
  std::vector<G4int>      CellStaggeringParity;
  std::vector<G4int>      CellStaggeringRemove;
  std::vector<G4double>   CellCrystalSizeX;
  std::vector<G4double>   CellCrystalSizeY;
  std::vector<G4double>   CellCrystalSizeZ;
  std::vector<G4double>   CellCrystalPitchX;
  std::vector<G4double>   CellCrystalPitchY;
  config.readIntoVect(CellNames,        "cell_name");
  config.readIntoVect(CellPositionX,    "cell_pos_x");
  config.readIntoVect(CellPositionY,    "cell_pos_y");
  config.readIntoVect(CellPositionZ,    "cell_pos_z");
  config.readIntoVect(CellXelements,    "cell_x_elements");
  config.readIntoVect(CellYelements,    "cell_y_elements");
  config.readIntoVect(CellMaterial,     "cell_crystal_material");
  config.readIntoVect(CellCrystalSizeX ,"cell_crystal_size_x");
  config.readIntoVect(CellCrystalSizeY ,"cell_crystal_size_y");
  config.readIntoVect(CellCrystalSizeZ ,"cell_crystal_size_z");
  config.readIntoVect(CellCrystalPitchX,"cell_crystal_pitch_x");
  config.readIntoVect(CellCrystalPitchY,"cell_crystal_pitch_y");
  // optional keys
  // air layer. if nothing is given, no air layer (so build 3 meaningless arrays)
  bool airLayerGiven = config.read<bool>("cell_air_layer",0);
  G4cout << "airLayerGiven = " << airLayerGiven << G4endl;
  if(airLayerGiven)
  {
    config.readIntoVect(CellAirLayer,     "cell_air_layer");
    config.readIntoVect(CellExtGapMaterial, "cell_ext_gap_material");
    //internal air layer? if it's not given, it means both layers are external
    bool intAirLayerGiven = config.read<bool>("cell_int_gap_material",0);
    G4cout << "intAirLayerGiven = " << intAirLayerGiven << G4endl;

    if(intAirLayerGiven) // /internal air layer? if it's not given, it means both layers are external
    {
      config.readIntoVect(CellIntGapMaterial,  "cell_int_gap_material");
    }
    else //if it's not given, it means both layers are external
    {
      config.readIntoVect(CellIntGapMaterial,  "cell_ext_gap_material");
    }
  }
  else
  {
    for(int i = 0; i < CellNames.size(); i++)
    {
      CellAirLayer      .push_back(0);
      CellExtGapMaterial.push_back(1);
      CellIntGapMaterial.push_back(1);
    }
  }
  // staggering. if nothing is given, no staggering
  bool stagGiven = config.read<bool>("cell_staggering",0);
  if(stagGiven)
  {
    config.readIntoVect(CellStaggering,         "cell_staggering");
    config.readIntoVect(CellStaggeringAxis,     "cell_staggering_axis");
    config.readIntoVect(CellStaggeringSize  ,   "cell_staggering_size");
    config.readIntoVect(CellStaggeringParity,   "cell_staggering_parity");
    config.readIntoVect(CellStaggeringRemove,   "cell_staggering_remove");
  }
  else
  {
    for(int i = 0; i < CellNames.size(); i++)
    {
      CellStaggering.push_back(0);
      CellStaggeringAxis.push_back(0);
      CellStaggeringSize.push_back(0);
      CellStaggeringParity.push_back(0);
      CellStaggeringRemove.push_back(0);
    }
  }
  // check that all vectors have the same length
  bool cell_allEquals = true;
  if(CellNames.size() != CellPositionX.size())        cell_allEquals = false;
  if(CellNames.size() != CellPositionY.size())        cell_allEquals = false;
  if(CellNames.size() != CellPositionZ.size())        cell_allEquals = false;
  if(CellNames.size() != CellMaterial.size())         cell_allEquals = false;
  if(CellNames.size() != CellXelements.size())        cell_allEquals = false;
  if(CellNames.size() != CellYelements.size())        cell_allEquals = false;
  if(CellNames.size() != CellAirLayer.size())         cell_allEquals = false;
  if(CellNames.size() != CellExtGapMaterial.size())   cell_allEquals = false;
  if(CellNames.size() != CellIntGapMaterial.size())   cell_allEquals = false;
  if(CellNames.size() != CellStaggering.size())       cell_allEquals = false;
  if(CellNames.size() != CellStaggeringAxis.size())   cell_allEquals = false;
  if(CellNames.size() != CellStaggeringSize.size())   cell_allEquals = false;
  if(CellNames.size() != CellStaggeringParity.size()) cell_allEquals = false;
  if(CellNames.size() != CellStaggeringRemove.size()) cell_allEquals = false;
  if(CellNames.size() != CellCrystalSizeX .size())    cell_allEquals = false;
  if(CellNames.size() != CellCrystalSizeY .size())    cell_allEquals = false;
  if(CellNames.size() != CellCrystalSizeZ .size())    cell_allEquals = false;
  if(CellNames.size() != CellCrystalPitchX.size())    cell_allEquals = false;
  if(CellNames.size() != CellCrystalPitchY.size())    cell_allEquals = false;
  if(!cell_allEquals)
  {
    G4cerr << "<DetectorConstruction>: Size of cell arrays in template file do not match!" << G4endl ;
    exit (-1) ;
  }
  // surface conditions
  crystal_lateral_depolishing = config.read<double>("crystal_lateral_depolishing",0);
  crystal_exit_depolishing    = config.read<double>("crystal_exit_depolishing",0);
  // create cells
  // for(unsigned int iAbs = 0 ; iAbs < AbsNames.size(); iAbs++)
  // {
    for(unsigned int iCell = 0 ; iCell < CellNames.size(); iCell++)
    {
      Cell cell;
      cell.SetID(iCell);
      cell.SetName(CellNames[iCell]);
      G4double modPositionZ;
      if(CellPositionZ[iCell] > 0)
      {
        modPositionZ = CellPositionZ[iCell] + (absorber.GetEsrThickness())/2.0;
      }
      if(CellPositionZ[iCell] < 0)
      {
        modPositionZ = CellPositionZ[iCell] - (absorber.GetEsrThickness())/2.0;
      }
      if(CellPositionZ[iCell] == 0)
      {
        modPositionZ = 0;
      }
      cell.SetPosition          (CellPositionX[iCell],
                                 CellPositionY[iCell],
                                 modPositionZ);
      cell.SetCrystalMaterial   (CellMaterial[iCell]);
      cell.SetXelements         (CellXelements[iCell]);
      cell.SetYelements         (CellYelements[iCell]);
      cell.SetNominalCrystalDimensions (CellCrystalSizeX[iCell],
                                 CellCrystalSizeY[iCell],
                                 CellCrystalSizeZ[iCell]);
      cell.SetCrystalPitch      (CellCrystalPitchX[iCell],
                                 CellCrystalPitchY[iCell]);
      cell.SetAirLayer          (CellAirLayer[iCell]);
      cell.SetExtGapMaterial    (CellExtGapMaterial[iCell]);
      cell.SetIntGapMaterial    (CellIntGapMaterial[iCell]);
      cell.SetStaggering        (CellStaggering[iCell]);
      cell.SetStaggeringAxis    (CellStaggeringAxis[iCell]);
      cell.SetStaggeringSize    (CellStaggeringSize[iCell]);
      cell.SetStaggeringParity  (CellStaggeringParity[iCell]);
      cell.SetStaggeringRemove  (CellStaggeringRemove[iCell]);
      cell.MakeCellStruture();
      absorber.AddCell(cell);
    }

    G4cout << "Absorber " << absorber.GetName()
          <<  " cells = " << absorber.GetNumberOfCells() << G4endl;

  gap_abs_interface_material = config.read<int>("gap_abs_interface_material",0);
  gap_interface_readout_material = config.read<int>("gap_interface_readout_material",0);


  gapSize = config.read<double>("gap_size",0.01); // default air gap is 10 micron

  InterfaceSizeX = AbsSizeX;
  InterfaceSizeY = AbsSizeY;
  InterfaceSizeZ = config.read<double>("interfaceLength",30);
  if(gapSize > 0) // rescale if there is a gap
  {
    InterfaceSizeZ = InterfaceSizeZ+2.0*gapSize;
  }
  ReadoutSizeX = AbsSizeX;
  ReadoutSizeY = AbsSizeY;
  ReadoutSizeZ = config.read<double>("readout_length",50);

  module_size_x = (AbsSizeX);
  module_size_y = (AbsSizeY);
  module_size_z = (AbsSizeZ + 2.0*InterfaceSizeZ + 2.0*ReadoutSizeZ);

  module_pos_x = new G4double*[modules_nx];
  module_pos_y = new G4double*[modules_nx];
  module_pos_z = new G4double*[modules_nx];
  for(int jMod = 0; jMod < modules_ny; jMod++)
  {
    module_pos_x[jMod] = new G4double[modules_ny];
    module_pos_y[jMod] = new G4double[modules_ny];
    module_pos_z[jMod] = new G4double[modules_ny];
  }

  for(int iMod = 0; iMod < modules_nx; iMod++)
  {
    for(int jMod = 0; jMod < modules_ny; jMod++)
    {
      module_pos_x[iMod][jMod] = 0.5*module_size_x + iMod*(module_size_x) - (module_size_x * (modules_nx) * 0.5);
      module_pos_y[iMod][jMod] = 0.5*module_size_y + jMod*(module_size_y) - (module_size_y * (modules_ny) * 0.5);
      module_pos_z[iMod][jMod] = 0.0;
    }
  }

  calorimeter_x = modules_nx*module_size_x;
  calorimeter_y = modules_ny*module_size_y;
  calorimeter_z = module_size_z;

  expHall_x = 2.0*calorimeter_x;
  expHall_y = 2.0*calorimeter_y;
  expHall_z = 2.0*calorimeter_z;

  B_field_IsInitialized = false ;

  initializeMaterials () ;
}


//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


DetectorConstruction::~DetectorConstruction ()
{}


  //---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

G4VPhysicalVolume* DetectorConstruction::Construct ()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl ;





  //------------------------------------
  //------------- Geometry -------------
  //------------------------------------

  // The experimental Hall
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4cout << ">>>>>>>>>>>>> Build Experimental Hall "<< G4endl;
  G4VSolid * worldS = new G4Box ("worldS", 0.5 * expHall_x, 0.5 * expHall_y, 0.5 * expHall_z) ;
  G4LogicalVolume * worldLV = new G4LogicalVolume (worldS, WoMaterial, "worldLV", 0, 0, 0) ;
  G4VPhysicalVolume * worldPV = new G4PVPlacement (0, G4ThreeVector (), worldLV, "World", 0, false, 0, checkOverlaps) ;
  G4VisAttributes* VisAttWorld = new G4VisAttributes (black) ;
  VisAttWorld->SetVisibility (true) ;
  VisAttWorld->SetForceWireframe (true) ;
  worldLV->SetVisAttributes (VisAttWorld) ;

  // The calorimeter
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  // the calorimeter is made of 1 or more identical modules. They are arranged in a grid, with
  // the possibility to "remove" some module (to make space for a beam pipe)
  G4cout << ">>>>>>>>>>>>> Build calorimeter "<< G4endl;
  G4VSolid * calorimeterS = new G4Box ("calorimeterS", 0.5 * calorimeter_x, 0.5 * calorimeter_y, 0.5 * calorimeter_z ) ;
  G4LogicalVolume * calorimeterLV = new G4LogicalVolume (calorimeterS, WoMaterial, "calorimeterLV") ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.), calorimeterLV, "calorimeterPV", worldLV, false, 0, checkOverlaps) ;
  G4VisAttributes* Calo_VisAtt = new G4VisAttributes(orange);
  Calo_VisAtt->SetVisibility(false);
  Calo_VisAtt->SetForceWireframe(wireFrame);
  calorimeterLV->SetVisAttributes(Calo_VisAtt);


  //-----------------------------------------//
  // Modules
  //-----------------------------------------//

  int iForbid = 1;
  int jForbid = 1;
  G4cout << ">>>>>>>>>>>>> Build modules "<< G4endl;
  for(int iMod = 0; iMod < modules_nx ; iMod++)
  {
    for(int jMod = 0; jMod < modules_ny ; jMod++)
    {

      //avoid the module is it's forbidden
      if((iMod == iForbid) && (jMod == jForbid))
      {
        continue;
      }
      // each module is made of 5 block elements:
      // 1x "absorber"
      // 2x "interfaces"
      // 2x "readouts"

      // The "absorber" block is the core of the module. Physically, it consists of the absorbing material
      // plus the holes to host the crystals, plus the crystals. Logically, it also includes the surfaces
      // defined between the various elements of the block. The holes can be grouped in 1 or more cells.
      // Cells are collections of hole/crystals with omogeneous type of crystal. The absorber can be divided
      // in two by a layer of reflector (ESR).

      // On both sides of the absober, there is an "interface" volume. It is composed of 1 or more light guides
      // and two separation gaps. Each gap is separating the light guide(s) from the adjacent volumes (crystals on
      // one side, pmts on the other) and it cam be made of air or other materials

      // Finally, two readout volumes are placed at the two ends of the module. they are made of the readout
      // detectors, pmts for now


      //-----------------------------------//
      //-----------------------------------//
      // The module box
      //-----------------------------------//
      //-----------------------------------//
      //declear a stringstream for composition of names
      std::stringstream sname;
      std::stringstream smodule;
      int module_number = iMod*modules_ny + jMod;
      smodule << "module_" << module_number;
      std::string moduleName = smodule.str();

      sname << moduleName << "_S";
      G4VSolid * moduleS = new G4Box (sname.str().c_str(), 0.5 * module_size_x, 0.5 * module_size_y, 0.5 * module_size_z ) ;
      sname.str("");
      sname << moduleName << "_LV";
      G4LogicalVolume * moduleLV = new G4LogicalVolume (moduleS, WoMaterial, sname.str().c_str()) ;
      sname.str("");
      sname << moduleName << "_PV";
      new G4PVPlacement (0, G4ThreeVector (module_pos_x[iMod][jMod],
                                           module_pos_y[iMod][jMod],
                                           module_pos_z[iMod][jMod]),
                                           moduleLV, sname.str().c_str(),
                                           calorimeterLV, false, 0, checkOverlaps) ;
      sname.str("");
      G4VisAttributes* Module_VisAtt = new G4VisAttributes(black);  // color
      Module_VisAtt->SetVisibility(false);
      Module_VisAtt->SetForceWireframe(true);
      moduleLV->SetVisAttributes(Module_VisAtt);
      //-----------------------------------//


      //-----------------------------------//
      //---------------------------------- //
      // The interface volumes
      //-----------------------------------//
      //-----------------------------------//

      //-----------------------------------//
      // positive z inteface
      //-----------------------------------//
      // prepare pointers
      G4VSolid        *Interface_Positive_Z_S;
      G4LogicalVolume *Interface_Positive_Z_LV;
      G4PVPlacement   *Interface_Positive_Z_PV;
      sname << "Interface_Positive_Z_"
            << moduleName
            << "_S";
      Interface_Positive_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * InterfaceSizeZ ) ;
      sname.str("");
      sname << "Interface_Positive_Z_"
            << moduleName
            << "_LV";
      Interface_Positive_Z_LV = new G4LogicalVolume (Interface_Positive_Z_S, MyMaterials::AirKiller(), sname.str().c_str()) ;
      sname.str("");

      sname << "Interface_Positive_Z_"
            << moduleName
            << "_PV";
      Interface_Positive_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                     0,
                                                                     +(AbsSizeZ+InterfaceSizeZ)/2.0),
                                                                     Interface_Positive_Z_LV, sname.str().c_str(),
                                                                     moduleLV, false, 0, checkOverlaps) ;
      sname.str("");
      //-----------------------------------//
      // elements of positive z interface
      //-----------------------------------//
      // place an air/grease layer at the beginning ad at the end of the interfaces
      // prepare pointers
      G4VSolid        *Gap_Abs_Interface_Positive_Z_S;
      G4LogicalVolume *Gap_Abs_Interface_Positive_Z_LV;
      G4PVPlacement   *Gap_Abs_Interface_Positive_Z_PV;
      sname << "Gap_Abs_Interface_Positive_Z_"
            << moduleName
            << "_S";
      Gap_Abs_Interface_Positive_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * gapSize ) ;
      sname.str("");
      sname << "Gap_Abs_Interface_Positive_Z_"
            << moduleName
            << "_LV";
      Gap_Abs_Interface_Positive_Z_LV = new G4LogicalVolume (Gap_Abs_Interface_Positive_Z_S, GapAbsToInterfaceMaterial, sname.str().c_str()) ;
      sname.str("");

      sname << "Gap_Abs_Interface_Positive_Z_"
            << moduleName
            << "_PV";
      Gap_Abs_Interface_Positive_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                     0,
                                                                     -(InterfaceSizeZ-gapSize)/2.0),
                                                                     Gap_Abs_Interface_Positive_Z_LV, sname.str().c_str(),
                                                                     Interface_Positive_Z_LV, false, 0, checkOverlaps) ;
      sname.str("");
      G4VSolid        *Gap_Interface_Readout_Positive_Z_S;
      G4LogicalVolume *Gap_Interface_Readout_Positive_Z_LV;
      G4PVPlacement   *Gap_Interface_Readout_Positive_Z_PV;
      sname << "Gap_Interface_Readout_Positive_Z_"
            << moduleName
            << "_S";
      Gap_Interface_Readout_Positive_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * gapSize ) ;
      sname.str("");
      sname << "Gap_Interface_Readout_Positive_Z_"
            << moduleName
            << "_LV";
      Gap_Interface_Readout_Positive_Z_LV = new G4LogicalVolume (Gap_Interface_Readout_Positive_Z_S, GapInterfaceToReadoutMaterial, sname.str().c_str()) ;
      sname.str("");

      sname << "Gap_Interface_Readout_Positive_Z_"
            << moduleName
            << "_PV";
      Gap_Interface_Readout_Positive_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                     0,
                                                                     +(InterfaceSizeZ-gapSize)/2.0),
                                                                     Gap_Interface_Readout_Positive_Z_LV, sname.str().c_str(),
                                                                     Interface_Positive_Z_LV, false, 0, checkOverlaps) ;
      sname.str("");
      // end of positive z
      //-----------------------------------//


      //-----------------------------------//
      // negative z
      //-----------------------------------//
      G4VSolid        *Interface_Negative_Z_S;
      G4LogicalVolume *Interface_Negative_Z_LV;
      G4PVPlacement   *Interface_Negative_Z_PV;
      sname << "Interface_Negative_Z_"
            << moduleName
            << "_S";
      Interface_Negative_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * InterfaceSizeZ ) ;
      sname.str("");

      sname << "Interface_Negative_Z_"
            << moduleName
            << "_LV";
      Interface_Negative_Z_LV = new G4LogicalVolume (Interface_Negative_Z_S, MyMaterials::AirKiller(), sname.str().c_str()) ;
      sname.str("");

      sname << "Interface_Negative_Z_"
            << moduleName
            << "_PV";
      Interface_Negative_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                     0,
                                                                     -(AbsSizeZ+InterfaceSizeZ)/2.0),
                                                                     Interface_Negative_Z_LV, sname.str().c_str(),
                                                                     moduleLV, false, 0, checkOverlaps) ;
      sname.str("");
      //
      G4VSolid        *Gap_Abs_Interface_Negative_Z_S;
      G4LogicalVolume *Gap_Abs_Interface_Negative_Z_LV;
      G4PVPlacement   *Gap_Abs_Interface_Negative_Z_PV;

      sname << "Gap_Abs_Interface_Negative_Z_"
            << moduleName
            << "_S";
      Gap_Abs_Interface_Negative_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * gapSize ) ;
      sname.str("");
      sname << "Gap_Abs_Interface_Negative_Z_"
            << moduleName
            << "_LV";
      Gap_Abs_Interface_Negative_Z_LV = new G4LogicalVolume (Gap_Abs_Interface_Negative_Z_S, GapAbsToInterfaceMaterial, sname.str().c_str()) ;
      sname.str("");

      sname << "Gap_Abs_Interface_Negative_Z_"
            << moduleName
            << "_PV";
      Gap_Abs_Interface_Negative_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                     0,
                                                                     +(InterfaceSizeZ-gapSize)/2.0),
                                                                     Gap_Abs_Interface_Negative_Z_LV, sname.str().c_str(),
                                                                     Interface_Negative_Z_LV, false, 0, checkOverlaps) ;
      sname.str("");
      G4VSolid        *Gap_Interface_Readout_Negative_Z_S;
      G4LogicalVolume *Gap_Interface_Readout_Negative_Z_LV;
      G4PVPlacement   *Gap_Interface_Readout_Negative_Z_PV;

      sname << "Gap_Interface_Readout_Negative_Z_"
            << moduleName
            << "_S";
      Gap_Interface_Readout_Negative_Z_S = new G4Box (sname.str().c_str(), 0.5 * InterfaceSizeX, 0.5 * InterfaceSizeY, 0.5 * gapSize ) ;
      sname.str("");
      sname << "Gap_Interface_Readout_Negative_Z_"
            << moduleName
            << "_LV";
      Gap_Interface_Readout_Negative_Z_LV = new G4LogicalVolume (Gap_Interface_Readout_Negative_Z_S, GapAbsToInterfaceMaterial, sname.str().c_str()) ;
      sname.str("");

      sname << "Gap_Interface_Readout_Negative_Z_"
            << moduleName
            << "_PV";
      Gap_Interface_Readout_Negative_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                     0,
                                                                     -(InterfaceSizeZ-gapSize)/2.0),
                                                                     Gap_Interface_Readout_Negative_Z_LV, sname.str().c_str(),
                                                                     Interface_Negative_Z_LV, false, 0, checkOverlaps) ;
      sname.str("");
      // end of negative z
      //-----------------------------------//

      // vis att.
      G4VisAttributes* Interface_VisAtt = new G4VisAttributes(red);  // color
      Interface_VisAtt->SetVisibility(true);      // absorbers are always visible
      Interface_VisAtt->SetForceWireframe(true);
      Interface_Positive_Z_LV->SetVisAttributes(Interface_VisAtt);
      Interface_Negative_Z_LV->SetVisAttributes(Interface_VisAtt);

      G4VisAttributes* Gap_VisAtt = new G4VisAttributes(cyan);  // color
      Gap_VisAtt->SetVisibility(true);      // absorbers are always visible
      Gap_VisAtt->SetForceWireframe(true);
      Gap_Abs_Interface_Positive_Z_LV->SetVisAttributes(Gap_VisAtt);
      Gap_Interface_Readout_Positive_Z_LV->SetVisAttributes(Gap_VisAtt);
      Gap_Abs_Interface_Negative_Z_LV->SetVisAttributes(Gap_VisAtt);
      Gap_Interface_Readout_Negative_Z_LV->SetVisAttributes(Gap_VisAtt);

      // end of the interfaces
      //--------------------------------------- //



      //--------------------------------------- //
      // The readout volumes
      // prepare pointers
      G4VSolid        *Readout_Positive_Z_S;
      G4LogicalVolume *Readout_Positive_Z_LV;
      G4PVPlacement   *Readout_Positive_Z_PV;
      sname << "Readout_Positive_Z_"
            << moduleName
            << "_S";
      Readout_Positive_Z_S = new G4Box (sname.str().c_str(), 0.5 * ReadoutSizeX, 0.5 * ReadoutSizeY, 0.5 * ReadoutSizeZ ) ;
      sname.str("");
      sname << "Readout_Positive_Z_"
            << moduleName
            << "_LV";
      Readout_Positive_Z_LV = new G4LogicalVolume (Readout_Positive_Z_S, MyMaterials::AirKiller(), sname.str().c_str()) ;
      sname.str("");

      sname << "Readout_Positive_Z_"
            << moduleName
            << "_PV";
      Readout_Positive_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                     0,
                                                                     +(AbsSizeZ+2.0*InterfaceSizeZ+ReadoutSizeZ)/2.0),
                                                                     Readout_Positive_Z_LV, sname.str().c_str(),
                                                                     moduleLV, false, 0, checkOverlaps) ;
      sname.str("");
      G4VSolid        *Readout_Negative_Z_S;
      G4LogicalVolume *Readout_Negative_Z_LV;
      G4PVPlacement   *Readout_Negative_Z_PV;

      sname << "Readout_Negative_Z_"
            << moduleName
            << "_S";
      Readout_Negative_Z_S = new G4Box (sname.str().c_str(), 0.5 * ReadoutSizeX, 0.5 * ReadoutSizeY, 0.5 * ReadoutSizeZ ) ;
      sname.str("");
      sname << "Readout_Negative_Z_"
            << moduleName
            << "_LV";
      Readout_Negative_Z_LV = new G4LogicalVolume (Readout_Negative_Z_S, MyMaterials::AirKiller(), sname.str().c_str()) ;
      sname.str("");

      sname << "Readout_Negative_Z_"
            << moduleName
            << "_PV";
      Readout_Negative_Z_PV = new G4PVPlacement (0, G4ThreeVector (0,
                                                                     0,
                                                                     -(AbsSizeZ+2.0*InterfaceSizeZ+ReadoutSizeZ)/2.0),
                                                                     Readout_Negative_Z_LV, sname.str().c_str(),
                                                                     moduleLV, false, 0, checkOverlaps) ;
      sname.str("");
      G4VisAttributes* Readout_VisAtt = new G4VisAttributes(green);  // color
      Readout_VisAtt->SetVisibility(true);      // absorbers are always visible
      Readout_VisAtt->SetForceWireframe(true);
      Readout_Positive_Z_LV->SetVisAttributes(Readout_VisAtt);
      Readout_Negative_Z_LV->SetVisAttributes(Readout_VisAtt);
      // end of the readout
      //--------------------------------------- //



      //--------------------------------------- //
      // Specific optical readout part
      // this is intentionally hardcoded (but the meaningful job is done in different routines)
      // because there is no way to modularize every possible readout
      // the user will have to take care of writing another function to create the volumes of
      // light concentrators and photo-detectors, and modify the CreateTree and SteppingAction classes
      // to have a personalized readout of the light transport
      // light guides
      ConstructLightGuides(Interface_Negative_Z_LV,moduleName,0);
      // on the other side, we need to rotate the light guides...
      G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
      rotationMatrix->rotateY(180.*deg);
      ConstructLightGuides(Interface_Positive_Z_LV,moduleName,rotationMatrix);
      // and pmts
      ConstructPMTs(Readout_Positive_Z_LV,moduleName,0);
      ConstructPMTs(Readout_Negative_Z_LV,moduleName,0);



      //--------------------------------------- //
      // The absorber
      //--------------------------------------- //
      // the absorber is 1, but if there is a cell longitudinal separation
      // it means that it is logically divided in two
      // now, if this separation is of type 0 or 1 (just air, or aluminization)
      // nothing need to be done here
      // otherwise, if the separation is type 2, with a foil of esr, the absorber need to be effectively
      // separated in two
      // to do this, it's created larger by a factor equal to the esr thickness
      // then a air volume is created in the middle

      // prepare pointers
      G4VSolid        *Absorber_S;
      G4LogicalVolume *Absorber_LV;
      G4PVPlacement   *Absorber_PV;
      // create the shape
      sname << absorber.GetName() << "_S"
            << "_" << moduleName ;
      Absorber_S    = new G4Box (sname.str().c_str(),
                             0.5 * absorber.GetSizeX(),
                             0.5 * absorber.GetSizeY(),
                             0.5 * (absorber.GetSizeZ()));
      sname.str("");
      // create the logical volume
      sname << absorber.GetName() << "_LV"
            << "_" << moduleName;
      Absorber_LV   = new G4LogicalVolume (Absorber_S,
                                         absorber.GetMaterial(),
                                         sname.str().c_str()) ;
      sname.str("");
      // create the physical volume (i.e. place it in space)
      sname << absorber.GetName() << "_PV"
            << "_" << moduleName;
      Absorber_PV  =  new G4PVPlacement (0,
                                     G4ThreeVector (absorber.GetPositionX(),
                                                    absorber.GetPositionY(),
                                                    absorber.GetPositionZ()),
                                     Absorber_LV,
                                     sname.str().c_str(),
                                     moduleLV,
                                     false, 0, checkOverlaps);
      sname.str("");
      G4VisAttributes* Absorber_VisAtt = new G4VisAttributes(blue);  // color
      Absorber_VisAtt->SetVisibility(true);      // absorbers are always visible
      Absorber_VisAtt->SetForceWireframe(wireFrame);
      Absorber_LV->SetVisAttributes(Absorber_VisAtt);


      // the absorber does not have an index of refrection, but it will be touched by optical photons
      // if no surface is defined, it will kill them. So we define a G4LogicalSkinSurface
      sname << "SkinSurf_" << absorber.GetName()
            << "_" << moduleName;
      G4OpticalSurface* absSurface = new G4OpticalSurface(sname.str().c_str());
      sname.str("");
      absSurface->SetType(dielectric_metal);
      absSurface->SetFinish(ground);
      absSurface->SetModel(unified);
      absSurface->SetSigmaAlpha(0.1);
      absSurface->SetMaterialPropertiesTable(MyMaterials::ABS_SURF());
      sname << "LSS_"
            << absorber.GetName() << "_LSS"
            << "_" << moduleName;
      G4LogicalSkinSurface *abs_LSS = new G4LogicalSkinSurface(sname.str().c_str(),Absorber_LV,absSurface);
      sname.str("");



      //-------------------------------//
      // Reflector diving the absorber
      //-------------------------------//
      G4VSolid        *esr_S;
      G4LogicalVolume *esr_LV;
      G4PVPlacement   *esr_PV;
      if(esr_thickness > 0) // only if reflector thickness is > 0, remember that it is set to 0 if the separation is made by aluminization
      {


        sname << "esr_S_"
              << absorber.GetName()
              << "_" << moduleName;
        esr_S    = new G4Box (sname.str().c_str(),
                               0.5 * absorber.GetSizeX(),
                               0.5 * absorber.GetSizeY(),
                               0.5 * (esr_thickness));
        sname.str("");
        // create the logical volume
        sname << "esr_LV_" << absorber.GetName()
              << "_" << moduleName;
        esr_LV   = new G4LogicalVolume (esr_S,
                                           MyMaterials::Air(),
                                           sname.str().c_str()) ;
        sname.str("");
        // create the physical volume (i.e. place it in space)
        sname << "esr_PV_"
              << absorber.GetName()
              << "_" << moduleName;;
        esr_PV  =  new G4PVPlacement (0,
                                       G4ThreeVector (0,
                                                      0,
                                                      cell_separator_position),
                                       esr_LV,
                                       sname.str().c_str(),
                                       Absorber_LV,
                                       false, 0, checkOverlaps);
        sname.str("");
        G4VisAttributes* esr_VisAtt = new G4VisAttributes(green);  // color
        esr_VisAtt->SetVisibility(true);      // absorbers are always visible
        esr_VisAtt->SetForceWireframe(wireFrame);
        esr_LV->SetVisAttributes(esr_VisAtt);
      }
      //--------------------------------------- //






      // run on cells of the absorber
      for(int iCell = 0; iCell < absorber.GetNumberOfCells() ; iCell++)
      {
        // cells are a collection of holes in an absorber
        // first, get the cell specifications
        Cell cell = absorber.GetCell(iCell);



        // now for each cell run on rows
        // get the number of rows of this cell
        G4int NumberOfLayers = cell.GetNumberOfLayers();
        for(int iLay = 0; iLay < NumberOfLayers ; iLay++)
        {
          //get the row
          row_t layer = cell.GetLayer(iLay);

          // loop on elements of the row, aka the columns, aka the individual holes (elements)
          int NumberOfLayerElements = layer.elements;
          for(int iEl = 0; iEl < NumberOfLayerElements ; iEl++)
          {

            // now, each element is made by a hole filled by a crystal.
            // however, we need also to specify the air gaps
            // the thickness of air gaps is defined by the user (airLayer key in cfg file)
            // the hole transvers dimensions are already made larger than the crystal, by a
            // factor 2*airLayer, to fill the crystal plus the air gap. But the hole cannot
            // be longer than the absorber (and hopefully the user didn't define a crystal, or the 3 cry, shorter
            // than a hole).
            // Here, we also need to create physical volumes for the air gaps, because we will need
            // optical surfaces (and also to put optical grease/glue if needed).



            //--------------------------------------- //
            // the hole
            G4VSolid        *hole_S;
            G4LogicalVolume *hole_LV;
            G4PVPlacement   *hole_PV;
            // create the hole shape
            sname << "hole_" << iLay << "_" << iEl << "_S"
                  << cell.GetName() << "_"
                  << absorber.GetName() << "_"
                  << moduleName;
            hole_S    = new G4Box (sname.str().c_str(),
                                   0.5 * layer.hole_dimension_x[iEl],
                                   0.5 * layer.hole_dimension_y[iEl],
                                   0.5 * layer.hole_dimension_z[iEl]);
            sname.str("");
            // create the logical volume
            sname << "hole_" << iLay << "_" << iEl << "_LV"
                  << cell.GetName() << "_"
                  << absorber.GetName() << "_"
                  << moduleName;
            hole_LV   = new G4LogicalVolume (hole_S,
                                             MyMaterials::Air(),     // always air
                                             sname.str().c_str()) ;
            sname.str("");
            // create the physical volume (i.e. place the hole in the absorber)
            sname << "hole_" << iLay << "_" << iEl << "_PV"
                  << cell.GetName() << "_"
                  << absorber.GetName() << "_"
                  << moduleName;
            hole_PV  =  new G4PVPlacement (0,
                                           G4ThreeVector (cell.GetPositionX()+layer.hole_pos_x[iEl],
                                                          cell.GetPositionY()+layer.hole_pos_y[iEl],
                                                          cell.GetPositionZ()+layer.hole_pos_z[iEl]),
                                           hole_LV,
                                           sname.str().c_str(),
                                           Absorber_LV,
                                           false, 0, checkOverlaps);
            sname.str("");
            // set the Visualization properties of the hole
            G4VisAttributes* Hole_VisAtt = new G4VisAttributes(air);
            Hole_VisAtt->SetVisibility(visibility);
            Hole_VisAtt->SetForceWireframe(wireFrame);
            hole_LV->SetVisAttributes(Hole_VisAtt);
            // end of the hole
            //--------------------------------------- //


            //--------------------------------------- //
            // create the crystal
            G4VSolid        *crystal_S;
            G4LogicalVolume *crystal_LV;
            G4PVPlacement   *crystal_PV;
            // G4cout << "-----------------------> here "<< G4endl;
            // create the shape
            sname << "crystal_" << iLay << "_" << iEl << "_S"
                  << cell.GetName() << "_"
                  << absorber.GetName() << "_"
                  << moduleName;
            crystal_S    = new G4Box (sname.str().c_str(),
                                      0.5 * layer.crystal_dimension_x[iEl],
                                      0.5 * layer.crystal_dimension_y[iEl],
                                      0.5 * layer.crystal_dimension_z[iEl]);
            sname.str("");

            // create the logical volume
            sname << "crystal_" << iLay << "_" << iEl << "_LV"
                  << cell.GetName() << "_"
                  << absorber.GetName() << "_"
                  << moduleName;
            crystal_LV   = new G4LogicalVolume (crystal_S,
                                                cell.GetCrystalMaterial(),
                                                sname.str().c_str()) ;
            sname.str("");

            // G4cout << "MATERIAL = " << crystal_LV->GetMaterial() << G4endl;

            // create the physical volume (i.e. place it in space)
            sname << "crystal_" << iLay << "_" << iEl << "_PV"
                  << cell.GetName() << "_"
                  << absorber.GetName() << "_"
                  << moduleName;
            crystal_PV  =  new G4PVPlacement (0,
                                           G4ThreeVector (layer.crystal_pos_x[iEl],
                                                          layer.crystal_pos_y[iEl],
                                                          layer.crystal_pos_z[iEl]),

                                           crystal_LV,
                                           sname.str().c_str(),
                                           hole_LV,   //mother is hole
                                           false, 0, checkOverlaps);
            sname.str("");

            G4VisAttributes* Cry_VisAtt = new G4VisAttributes(cell.GetCrystalColor());
            Cry_VisAtt->SetVisibility(visibility);
            Cry_VisAtt->SetForceWireframe(wireFrame);
            crystal_LV->SetVisAttributes(Cry_VisAtt);
            // end of the crystal
            //--------------------------------------- //


            // opt surf is depolishing of exit faces
            if(crystal_exit_depolishing > 0)
            {
              //
              // surf to cry air gap positive z interface
              sname << "surf_CryToPositiveInterface1" << iLay << "_" << iEl << "_"
                    << cell.GetName() << "_"
                    << absorber.GetName() << "_"
                    << moduleName;
              G4OpticalSurface* surf_CryToPositiveInterface1 = new G4OpticalSurface(sname.str().c_str());
              surf_CryToPositiveInterface1->SetType(dielectric_dielectric);
              surf_CryToPositiveInterface1->SetFinish(ground);
              surf_CryToPositiveInterface1->SetModel(unified);
              surf_CryToPositiveInterface1->SetSigmaAlpha(crystal_exit_depolishing);
              surf_CryToPositiveInterface1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
              G4LogicalBorderSurface* l_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                            crystal_PV,
                                                                            Gap_Abs_Interface_Positive_Z_PV,
                                                                            surf_CryToPositiveInterface1);
              sname.str("");
              //
              sname << "surf_CryToPositiveInterface2" << iLay << "_" << iEl << "_"
                    << cell.GetName() << "_"
                    << absorber.GetName() << "_"
                    << moduleName;
              G4OpticalSurface* surf_CryToPositiveInterface2 = new G4OpticalSurface(sname.str().c_str());
              surf_CryToPositiveInterface2->SetType(dielectric_dielectric);
              surf_CryToPositiveInterface2->SetFinish(ground);
              surf_CryToPositiveInterface2->SetModel(unified);
              surf_CryToPositiveInterface2->SetSigmaAlpha(crystal_exit_depolishing);
              surf_CryToPositiveInterface2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
              G4LogicalBorderSurface* l_border2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                            Gap_Abs_Interface_Positive_Z_PV,
                                                                            crystal_PV,
                                                                            surf_CryToPositiveInterface2);
              sname.str("");



              // surf to cry air gap negative z interface
              sname << "surf_CryToNegativeInterface1" << iLay << "_" << iEl << "_"
                    << cell.GetName() << "_"
                    << absorber.GetName() << "_"
                    << moduleName;
              G4OpticalSurface* surf_CryToNegativeInterface1 = new G4OpticalSurface(sname.str().c_str());
              surf_CryToNegativeInterface1->SetType(dielectric_dielectric);
              surf_CryToNegativeInterface1->SetFinish(ground);
              surf_CryToNegativeInterface1->SetModel(unified);
              surf_CryToNegativeInterface1->SetSigmaAlpha(crystal_exit_depolishing);
              surf_CryToNegativeInterface1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
              G4LogicalBorderSurface* l_border3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                            crystal_PV,
                                                                            Gap_Abs_Interface_Negative_Z_PV,
                                                                            surf_CryToNegativeInterface1);
              sname.str("");
              //
              sname << "surf_CryToNegativeInterface2" << iLay << "_" << iEl << "_"
                    << cell.GetName() << "_"
                    << absorber.GetName() << "_"
                    << moduleName;
              G4OpticalSurface* surf_CryToNegativeInterface2 = new G4OpticalSurface(sname.str().c_str());
              surf_CryToNegativeInterface2->SetType(dielectric_dielectric);
              surf_CryToNegativeInterface2->SetFinish(ground);
              surf_CryToNegativeInterface2->SetModel(unified);
              surf_CryToNegativeInterface2->SetSigmaAlpha(crystal_exit_depolishing);
              surf_CryToNegativeInterface2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
              G4LogicalBorderSurface* l_border4 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                            Gap_Abs_Interface_Negative_Z_PV,
                                                                            crystal_PV,
                                                                            surf_CryToNegativeInterface2);
              sname.str("");
            }








            // add air gap on internal side, only if there is an internal side...
            if(cell.GetAirLayer() > 0)
            {
              G4double intGapPosition;
              if(cell.GetPositionZ() != 0)
              {
                if(cell.GetPositionZ() >  0) // then the exit must be in positive position
                {
                  // extGapPosition = + (cell.GetHoleDimensionZ()/2.0 - cell.GetAirLayer() / 2.0 );
                  intGapPosition = - (cell.GetHoleDimensionZ()/2.0 - cell.GetAirLayer() / 2.0 );
                }
                else  // then the exit must be in negative position
                {
                  // extGapPosition = - (cell.GetHoleDimensionZ()/2.0 - cell.GetAirLayer() / 2.0 );
                  intGapPosition = + (cell.GetHoleDimensionZ()/2.0 - cell.GetAirLayer() / 2.0 );
                }


                //--------------------------------------- //
                // the internal air gap
                G4VSolid        *air_gap_int_S;
                G4LogicalVolume *air_gap_int_LV;
                G4PVPlacement   *air_gap_int_PV;

                // create the shape
                // section as the hole, thickness as the airLayer
                sname << "air_gap_int_" << iLay << "_" << iEl << "_S"
                      << cell.GetName() << "_"
                      << absorber.GetName() << "_"
                      << moduleName;
                air_gap_int_S  = new G4Box (sname.str().c_str(),
                                             0.5*layer.hole_dimension_x[iEl],
                                             0.5*layer.hole_dimension_y[iEl],
                                             0.5 * cell.GetAirLayer());
                sname.str("");
                // create the logical volume
                sname << "air_gap_int_" << iLay << "_" << iEl << "_LV"
                      << cell.GetName() << "_"
                      << absorber.GetName() << "_"
                      << moduleName;
                air_gap_int_LV  = new G4LogicalVolume (air_gap_int_S,
                                                        cell.GetIntGapMaterial(),
                                                        sname.str().c_str()) ;
                sname.str("");
                // create the physical volume (i.e. place it in space)
                sname << "air_gap_int_" << iLay << "_" << iEl << "_PV"
                      << cell.GetName() << "_"
                      << absorber.GetName() << "_"
                      << moduleName;
                air_gap_int_PV = new G4PVPlacement  (0,
                                                      G4ThreeVector (0,0,
                                                      intGapPosition),
                                                      air_gap_int_LV,
                                                      sname.str().c_str(),
                                                      hole_LV,    //mother is hole
                                                      false, 0, checkOverlaps);
                sname.str("");

                G4VisAttributes* air_gap_int_VisAtt = new G4VisAttributes(cell.GetIntGapMaterialColor());
                air_gap_int_VisAtt->SetVisibility(visibility);
                air_gap_int_VisAtt->SetForceWireframe(wireFrame);
                air_gap_int_LV->SetVisAttributes(air_gap_int_VisAtt);
                // end of air gap "back"
                //--------------------------------------- //

                // optical surfaces
                // define the optical surf that separates cells longitudinally
                // sepation can be
                // 0) nothing = air
                // 1) aluminization
                // 2) esr


                if(cell_separation_type == 1)  // aluminization
                {
                // CASE 0)
                // if exit surface is polished ----> nothing to do
                // if exit surface depolished  ----> set a depolished surface between the crystal and the air layer just added
                  if(crystal_exit_depolishing > 0)
                  {
                    //
                    // surf to int air gap
                    sname << "intSurf1_" << iLay << "_" << iEl << "_"
                          << cell.GetName() << "_"
                          << absorber.GetName() << "_"
                          << moduleName;
                    G4OpticalSurface* int_opt_surf1 = new G4OpticalSurface(sname.str().c_str());
                    int_opt_surf1->SetType(dielectric_dielectric);
                    int_opt_surf1->SetFinish(ground);
                    int_opt_surf1->SetModel(unified);
                    int_opt_surf1->SetSigmaAlpha(crystal_exit_depolishing);
                    int_opt_surf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
                    G4LogicalBorderSurface* l_border3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                                  crystal_PV,
                                                                                  air_gap_int_PV,
                                                                                  int_opt_surf1);
                    sname.str("");
                    //
                    sname << "intSurf2_" << iLay << "_" << iEl << "_"
                          << cell.GetName() << "_"
                          << absorber.GetName() << "_"
                          << moduleName;
                    G4OpticalSurface* int_opt_surf2 = new G4OpticalSurface(sname.str().c_str());
                    int_opt_surf2->SetType(dielectric_dielectric);
                    int_opt_surf2->SetFinish(ground);
                    int_opt_surf2->SetModel(unified);
                    int_opt_surf2->SetSigmaAlpha(crystal_exit_depolishing);
                    int_opt_surf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
                    G4LogicalBorderSurface* l_border4 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                                  air_gap_int_PV,
                                                                                  crystal_PV,
                                                                                  int_opt_surf2);
                    sname.str("");
                  }
                }

                // CASE 1)
                // aluminization
                if(cell_separation_type == 1)  // aluminization
                {
                  // aluminization means polished front painted
                  //
                  sname << "extSurf1_" << iLay << "_" << iEl << "_"
                        << cell.GetName() << "_"
                        << absorber.GetName() << "_"
                        << moduleName;
                  G4OpticalSurface* opt_surf1 = new G4OpticalSurface(sname.str().c_str());
                  opt_surf1->SetType(dielectric_dielectric);
                  opt_surf1->SetFinish(polishedfrontpainted);
                  opt_surf1->SetModel(unified);
                  // opt_surf1->SetSigmaAlpha(crystal_exit_depolishing);
                  opt_surf1->SetMaterialPropertiesTable(MyMaterials::ESR());
                  G4LogicalBorderSurface* l_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                                 crystal_PV,
                                                                                 air_gap_int_PV,
                                                                                 opt_surf1);
                  sname.str("");
                  //
                }


                if(cell_separation_type == 2)  // aluminization
                {
                  // CASE 2)
                  // esr
                  // if this is the case, there MUST be a esr volume
                  // first, the depolishing, if it's there
                  if(crystal_exit_depolishing > 0)
                  {
                    //
                    // surf to int air gap
                    sname << "intSurf1_" << iLay << "_" << iEl << "_"
                          << cell.GetName() << "_"
                          << absorber.GetName() << "_"
                          << moduleName;
                    G4OpticalSurface* int_opt_surf1 = new G4OpticalSurface(sname.str().c_str());
                    int_opt_surf1->SetType(dielectric_dielectric);
                    int_opt_surf1->SetFinish(ground);
                    int_opt_surf1->SetModel(unified);
                    int_opt_surf1->SetSigmaAlpha(crystal_exit_depolishing);
                    int_opt_surf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
                    G4LogicalBorderSurface* l_border3 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                                  crystal_PV,
                                                                                  air_gap_int_PV,
                                                                                  int_opt_surf1);
                    sname.str("");
                    //
                    sname << "intSurf2_" << iLay << "_" << iEl << "_"
                          << cell.GetName() << "_"
                          << absorber.GetName() << "_"
                          << moduleName;
                    G4OpticalSurface* int_opt_surf2 = new G4OpticalSurface(sname.str().c_str());
                    int_opt_surf2->SetType(dielectric_dielectric);
                    int_opt_surf2->SetFinish(ground);
                    int_opt_surf2->SetModel(unified);
                    int_opt_surf2->SetSigmaAlpha(crystal_exit_depolishing);
                    int_opt_surf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
                    G4LogicalBorderSurface* l_border4 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                                  air_gap_int_PV,
                                                                                  crystal_PV,
                                                                                  int_opt_surf2);
                    sname.str("");
                  }

                  // then, the esr surface between air_gap_int_PV and esr volume (just in one direction)
                  // aluminization means polished front painted
                  //
                  sname << "esrSurf1_" << iLay << "_" << iEl << "_"
                        << cell.GetName() << "_"
                        << absorber.GetName() << "_"
                        << moduleName;
                  G4OpticalSurface* esrSurf1 = new G4OpticalSurface(sname.str().c_str());
                  esrSurf1->SetType(dielectric_metal);
                  esrSurf1->SetFinish(polished);
                  esrSurf1->SetModel(unified);
                  esrSurf1->SetMaterialPropertiesTable(MyMaterials::ESR());
                  G4LogicalBorderSurface* l_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                                 air_gap_int_PV,
                                                                                 esr_PV,
                                                                                 esrSurf1);
                  sname.str("");


                }

              }
            } // end if airLayer


            // optical surfaces
            // optical photons are going to meet, in the cell
            // 1) surface crystal - hole when they exit from the side  -> surface could be
            //    depolished, so need for surface dielectric_dielectric with depolishing,
            //    if user requires it - back and forth!
            // 2) surface  crystal - ext gap material                  -> surface could be
            //    depolished, so need for surface dielectric_dielectric with depolishing,
            //    if user requires it - back and forth!
            // 3) surface  crystal - int gap material                  -> surface could be
            //    depolished, so need for surface dielectric_dielectric with depolishing,
            //    if user requires it - back and forth!
            // 4) surface int - reflector   -> if there is a reflector between modules
            // 5) surface int               -> if there is aluminization


            //---------------------------------------- //
            // 1) surface crystal - hole (crystal lateral surf)
            if(crystal_lateral_depolishing > 0)
            {
              sname << "latSurf1_" << iLay << "_" << iEl << "_"
                    << cell.GetName() << "_"
                    << absorber.GetName() << "_"
                    << moduleName;
              G4OpticalSurface* opt_surf1 = new G4OpticalSurface(sname.str().c_str());
              opt_surf1->SetType(dielectric_dielectric);
              opt_surf1->SetFinish(ground);
              opt_surf1->SetModel(unified);
              opt_surf1->SetSigmaAlpha(crystal_lateral_depolishing);
              opt_surf1->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
              G4LogicalBorderSurface* l_border1 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                            crystal_PV,
                                                                            hole_PV,
                                                                            opt_surf1);
              sname.str("");
              //
              sname << "latSurf2_" << iLay << "_" << iEl << "_"
                    << cell.GetName() << "_"
                    << absorber.GetName() << "_"
                    << moduleName;
              G4OpticalSurface* opt_surf2 = new G4OpticalSurface(sname.str().c_str());
              opt_surf2->SetType(dielectric_dielectric);
              opt_surf2->SetFinish(ground);
              opt_surf2->SetModel(unified);
              opt_surf2->SetSigmaAlpha(crystal_lateral_depolishing);
              opt_surf2->SetMaterialPropertiesTable(MyMaterials::crystal_depo_SURF());
              G4LogicalBorderSurface* l_border2 = new G4LogicalBorderSurface(sname.str().c_str(),
                                                                            hole_PV,
                                                                            crystal_PV,
                                                                            opt_surf2);
              sname.str("");
              //
            }
            // end of 1)
            //---------------------------------------- //
          }// end of loop on columns
        }//end loop on rows
      }// end loop on cells
    }//end loop on modules
  }








  //PG call the magnetic field initialisation
  if (B_field_intensity > 0.1 * tesla) ConstructField () ;

  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl ;
  return worldPV ;
}



//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void DetectorConstruction::initializeMaterials ()
{
  // define materials


  WoMaterial = NULL ;
  if      ( world_material == 1 ) WoMaterial = MyMaterials::Air () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre world material specifier " << world_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Wo. material: "<< WoMaterial << G4endl ;

  GlueMaterial = NULL;
  if      ( glue_interface == 0 ) GlueMaterial = MyMaterials::Air () ;
  else if ( glue_interface == 1 ) GlueMaterial = MyMaterials::OpticalGrease();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid glue material specifier " << glue_interface << G4endl ;
    exit (-1) ;
  }
  G4cout << "Glue material: "<< GlueMaterial << G4endl ;

  ConeMaterial  = NULL;
  if      ( cone_material == 0 ) ConeMaterial = MyMaterials::Air () ;
  else if ( cone_material == 1 ) ConeMaterial = MyMaterials::PLEX () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid cone material specifier " << cone_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Cone material: "<< ConeMaterial << G4endl ;

  GapAbsToInterfaceMaterial  = NULL;
  if      ( gap_abs_interface_material == 0 ) GapAbsToInterfaceMaterial = MyMaterials::Air () ;
  else if ( gap_abs_interface_material == 1 ) GapAbsToInterfaceMaterial = MyMaterials::OpticalGrease () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid GapAbsToInterfaceMaterial material specifier " << gap_abs_interface_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Cone material: "<< GapAbsToInterfaceMaterial << G4endl ;

  GapInterfaceToReadoutMaterial  = NULL;
  if      ( gap_interface_readout_material == 0 ) GapInterfaceToReadoutMaterial = MyMaterials::Air () ;
  else if ( gap_interface_readout_material == 1 ) GapInterfaceToReadoutMaterial = MyMaterials::OpticalGrease () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid GapInterfaceToReadoutMaterial material specifier " << gap_interface_readout_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "GapInterfaceToReadoutMaterial material: "<< GapInterfaceToReadoutMaterial << G4endl ;

  // AbMaterial = NULL ;
  // if      ( abs_material == 1 ) AbMaterial = MyMaterials::Brass () ;
  // else if ( abs_material == 2 ) AbMaterial = MyMaterials::Tungsten () ;
  // else if ( abs_material == 3 ) AbMaterial = MyMaterials::Lead () ;
  // else if ( abs_material == 4 ) AbMaterial = MyMaterials::Iron () ;
  // else if ( abs_material == 5 ) AbMaterial = MyMaterials::Aluminium () ;
  // else if ( abs_material == 6 ) AbMaterial = MyMaterials::CopperTungstenAlloy(W_fraction) ;
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid absorber material specifier " << abs_material << G4endl ;
  //   exit (-1) ;
  // }
  // G4cout << "Ab. material: "<< AbMaterial << G4endl ;

  ///*************88    ********************8  ******************


  // AbMaterial2 = NULL ;
  // if      ( Second_abs_material == 1 ) AbMaterial2 = MyMaterials::Brass () ;
  // else if ( Second_abs_material == 2 ) AbMaterial2 = MyMaterials::Tungsten () ;
  // else if ( Second_abs_material == 3 ) AbMaterial2 = MyMaterials::Lead () ;
  // else if ( Second_abs_material == 4 ) AbMaterial2 = MyMaterials::Iron () ;
  // else if ( Second_abs_material == 5 ) AbMaterial2 = MyMaterials::Aluminium () ;
  // else if ( Second_abs_material == 6 ) AbMaterial2 = MyMaterials::CopperTungstenAlloy(W_fraction) ;
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid absorber material specifier " << Second_abs_material << G4endl ;
  //   exit (-1) ;
  // }
  // G4cout << "Ab. material2: "<< AbMaterial2 << G4endl ;
  //
  PLEXMaterial = NULL ;
  PLEXMaterial = MyMaterials::PLEX () ;
  //
  G4cout << "PLEX. material: "<< PLEXMaterial << G4endl ;
  //
  PVCMaterial = NULL ;
  PVCMaterial = MyMaterials::PVC () ;
  //
  G4cout << "PVC. material: "<< PVCMaterial << G4endl ;
  //
  // PMTMaterial = NULL ;
  // PMTMaterial = MyMaterials::CuAir () ;
  //
  // G4cout << "PMT. material: "<< PMTMaterial << G4endl ;
  //
  // WiresMaterial = NULL ;
  // WiresMaterial = MyMaterials::Cu () ;
  //
  // G4cout << "Wires. material: "<< WiresMaterial << G4endl ;
  //
  // PlaneMaterial = NULL ;
  // if      ( lp_mat == 1 ) PlaneMaterial = MyMaterials::Brass () ;
  // else if ( lp_mat == 2 ) PlaneMaterial = MyMaterials::Tungsten () ;
  // else if ( lp_mat == 3 ) PlaneMaterial = MyMaterials::Lead () ;
  // else if ( lp_mat == 4 ) PlaneMaterial = MyMaterials::Iron () ;
  // else if ( lp_mat == 5 ) PlaneMaterial = MyMaterials::Aluminium () ;
  // else if ( lp_mat == 6 ) PlaneMaterial = MyMaterials::CopperTungstenAlloy(W_fraction) ;
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid absorber material specifier " << lp_mat << G4endl ;
  //   exit (-1) ;
  // }
  // G4cout << "Plane material: "<< PlaneMaterial << G4endl ;
  //
  //
  // // ***********  ****************  ****************   ***********
  //
  //
  // ClMaterial = NULL ;
  // if      ( fibre_material == 1 ) ClMaterial = MyMaterials::Quartz () ;
  // else if ( fibre_material == 2 ) ClMaterial = MyMaterials::SiO2_Ce () ;
  // else if ( fibre_material == 3 ) ClMaterial = MyMaterials::DSB_Ce () ;
  // else if ( fibre_material == 4 ) ClMaterial = MyMaterials::LuAG_Ce () ;
  // else if ( fibre_material == 5 ) ClMaterial = MyMaterials::YAG_Ce () ;
  // else if ( fibre_material == 6 ) ClMaterial = MyMaterials::GAGG_Ce_Mg() ;
  // else if ( fibre_material == 7 ) ClMaterial = MyMaterials::Water() ;
  //
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << fibre_material << G4endl ;
  //   exit (-1) ;
  // }
  // G4cout << "Cl. material: "<< ClMaterial << G4endl ;


  //____________________________________________________________

  // ClSSMaterial = NULL ;
  // if      ( fibre_material1 == 1 ) ClSSMaterial = MyMaterials::Quartz () ;
  // else if ( fibre_material1 == 2 ) ClSSMaterial = MyMaterials::SiO2_Ce () ;
  // else if ( fibre_material1 == 3 ) ClSSMaterial = MyMaterials::DSB_Ce () ;
  // else if ( fibre_material1 == 4 ) ClSSMaterial = MyMaterials::LuAG_Ce () ;
  // else if ( fibre_material1 == 5 ) ClSSMaterial = MyMaterials::YAG_Ce () ;
  // else if ( fibre_material1 == 6 ) ClSSMaterial = MyMaterials::GAGG_Ce_Mg() ;
  // else if ( fibre_material1 == 7 ) ClSSMaterial = MyMaterials::Water() ;
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << fibre_material1 << G4endl ;
  //   exit (-1) ;
  // }
  // G4cout << "ClSS. material: "<< ClSSMaterial << G4endl ;
  //
  // ClSSSMaterial = NULL ;
  // if      ( fibre_material2 == 1 ) ClSSSMaterial = MyMaterials::Quartz () ;
  // else if ( fibre_material2 == 2 ) ClSSSMaterial = MyMaterials::SiO2_Ce () ;
  // else if ( fibre_material2 == 3 ) ClSSSMaterial = MyMaterials::DSB_Ce () ;
  // else if ( fibre_material2 == 4 ) ClSSSMaterial = MyMaterials::LuAG_Ce () ;
  // else if ( fibre_material2 == 5 ) ClSSSMaterial = MyMaterials::YAG_Ce () ;
  // else if ( fibre_material2 == 6 ) ClSSSMaterial = MyMaterials::GAGG_Ce_Mg() ;
  // else if ( fibre_material2 == 7 ) ClSSSMaterial = MyMaterials::Water() ;
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << fibre_material1 << G4endl ;
  //   exit (-1) ;
  // }
  // G4cout << "2nd_Sect Corners fibres Material:"<< ClSSMaterial << G4endl ;
  //
  //
  // // *************    ***************   *****************  * *******
  // Cl3Material = NULL ;
  // if      ( Second_fibre_material == 1 ) Cl3Material = MyMaterials::Quartz () ;
  // else if ( Second_fibre_material == 2 ) Cl3Material = MyMaterials::SiO2_Ce () ;
  // else if ( Second_fibre_material == 3 ) Cl3Material = MyMaterials::DSB_Ce () ;
  // else if ( Second_fibre_material == 4 ) Cl3Material = MyMaterials::LuAG_Ce () ;
  // else if ( Second_fibre_material == 5 ) Cl3Material = MyMaterials::YAG_Ce () ;
  // else if ( Second_fibre_material == 6 ) Cl3Material = MyMaterials::GAGG_Ce_Mg() ;
  // else if ( Second_fibre_material == 7 ) Cl3Material = MyMaterials::Water() ;
  //
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Second_fibre_material << G4endl ;
  //   exit (-1) ;
  // }
  // G4cout << "Cl. material: "<< Cl3Material << G4endl ;
  //
  //
  // //____________________________________________________________
  //
  // Cl4SSMaterial = NULL ;
  // if      ( Second_fibre_material1 == 1 ) Cl4SSMaterial = MyMaterials::Quartz () ;
  // else if ( Second_fibre_material1 == 2 ) Cl4SSMaterial = MyMaterials::SiO2_Ce () ;
  // else if ( Second_fibre_material1 == 3 ) Cl4SSMaterial = MyMaterials::DSB_Ce () ;
  // else if ( Second_fibre_material1 == 4 ) Cl4SSMaterial = MyMaterials::LuAG_Ce () ;
  // else if ( Second_fibre_material1 == 5 ) Cl4SSMaterial = MyMaterials::YAG_Ce () ;
  // else if ( Second_fibre_material1 == 6 ) Cl4SSMaterial = MyMaterials::GAGG_Ce_Mg() ;
  // else if ( Second_fibre_material1 == 7 ) Cl4SSMaterial = MyMaterials::Water() ;
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Second_fibre_material1 << G4endl ;
  //   exit (-1) ;
  // }
  // G4cout << "ClSS. material: "<< Cl4SSMaterial << G4endl ;
  //
  //
  // Cl43SSMaterial = NULL ;
  // if      ( Second_fibre_material2 == 1 ) Cl43SSMaterial = MyMaterials::Quartz () ;
  // else if ( Second_fibre_material2 == 2 ) Cl43SSMaterial = MyMaterials::SiO2_Ce () ;
  // else if ( Second_fibre_material2 == 3 ) Cl43SSMaterial = MyMaterials::DSB_Ce () ;
  // else if ( Second_fibre_material2 == 4 ) Cl43SSMaterial = MyMaterials::LuAG_Ce () ;
  // else if ( Second_fibre_material2 == 5 ) Cl43SSMaterial = MyMaterials::YAG_Ce () ;
  // else if ( Second_fibre_material2 == 6 ) Cl43SSMaterial = MyMaterials::GAGG_Ce_Mg() ;
  // else if ( Second_fibre_material2 == 7 ) Cl43SSMaterial = MyMaterials::Water() ;
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Second_fibre_material1 << G4endl ;
  //   exit (-1) ;
  // }
  // G4cout << "2nd_Sect Corners fibres Material: "<< Cl4SSMaterial << G4endl ;





  //  *****************   ********************   ***************


  //______________________________________________________________

  // GaMaterial = NULL;
  // if     ( gap_material == 1 ) GaMaterial = MyMaterials::Air();
  // else if( gap_material == 2 ) GaMaterial = MyMaterials::OpticalGrease();
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid gap material specifier " << gap_material << G4endl;
  //   exit(-1);
  // }
  // G4cout << "Gap material: " << gap_material << G4endl;
  //
  //
  // DeMaterial = NULL;
  // if( det_material == 1 ) DeMaterial = MyMaterials::Silicon();
  // else
  // {
  //   G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid detector material specifier " << det_material << G4endl;
  //   exit(-1);
  // }
  // G4cout << "Detector material: " << det_material << G4endl;



  // if( fibre_absLength >= 0 )
  // {
  //   const G4int nEntries_ABS = 2;
  //   G4double PhotonEnergy_ABS[nEntries_ABS] = { 1.*eV, 10.*eV };
  //   G4double Absorption[nEntries_ABS] = { fibre_absLength*mm, fibre_absLength*mm };
  //
  //   ClMaterial->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
  //   ClMaterial->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH",PhotonEnergy_ABS,Absorption,nEntries_ABS);
  //   //_______________________________
  //
  //   ClSSMaterial->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
  //   ClSSMaterial->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH",PhotonEnergy_ABS,Absorption,nEntries_ABS);
  //   //___________________________________________________________
  // }

}


//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----


void DetectorConstruction::ConstructField ()
{
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::begin <<<<<<" << G4endl ;
  static G4TransportationManager * trMgr = G4TransportationManager::GetTransportationManager () ;

  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager * globalFieldMgr = trMgr->GetFieldManager () ;

  if (!B_field_IsInitialized)
  {
    // magnetic field parallel to the beam direction (w/ tilt)
    G4ThreeVector fieldVector (
      0.0522 * B_field_intensity,
      0.0522 * B_field_intensity,
      0.9973 * B_field_intensity
    ) ;

    B_field = new G4UniformMagField (fieldVector) ;
    globalFieldMgr->SetDetectorField (B_field) ;
    globalFieldMgr->CreateChordFinder (B_field) ;
    globalFieldMgr->GetChordFinder ()->SetDeltaChord (0.005 * mm) ;
    B_field_IsInitialized = true ;
  }
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::end <<< " << G4endl ;
  return ;
}

void DetectorConstruction::ConstructLightGuides(G4LogicalVolume *Interface_LV, std::string moduleName,G4RotationMatrix *rot)
{
  // now the light guides
  // Light guide structure
  // Truncated of cone, cut on the sides
  // goes from 20x20 mm2 (on the fibers) to 5 mm radius on the PMT
  // first, do a simple truncated cone
  std::stringstream sname;
  G4double lguide_edge = 20 *mm;
  G4double PLEX_depth = InterfaceSizeZ - 2.0*gapSize;
  G4double PMT_radius = 5 *mm;
  G4VSolid* coneS = new G4Cons("aCone",0,PMT_radius,0, 0.5*sqrt(2)*lguide_edge ,0.5*PLEX_depth,0.0 * deg,360.0* deg);
  // then prepare a hollow box, section of the hole 20x20 mm2
  G4Box* innerBoxS = new G4Box ("innerBoxS", 0.5*lguide_edge, 0.5*lguide_edge, 0.5*PLEX_depth) ;
  G4Box* outerBoxS = new G4Box ("outerBoxS", 0.5*sqrt(2)*lguide_edge, 0.5*sqrt(2)*lguide_edge, 0.5*PLEX_depth) ;
  G4VSolid* subtract = new G4SubtractionSolid("Hollow-Box", outerBoxS, innerBoxS,0,G4ThreeVector(0.,0.,0.));
  // subtract hollow box from truncated cone
  G4VSolid* coneSolid = new G4SubtractionSolid("coneSolid", coneS, subtract,0,  G4ThreeVector(0.,0.,0.));
  G4LogicalVolume * LguideLV = new G4LogicalVolume (coneSolid,ConeMaterial,"LguideLV");

  int nLGx = 3;
  int nLGy = 3;
  G4double lg_pitch = 20 * mm;
  G4double lg_x[3] = {-lg_pitch,0,lg_pitch};
  G4double lg_y[3] = {-lg_pitch,0,lg_pitch};
  for(int iLG = 0; iLG < nLGx; iLG++)
  {
    for(int jLG = 0; jLG < nLGy; jLG++)
    {
      sname << "LightGuide_"
            << iLG << "_" << jLG << "_"
            << moduleName
            << "_PV";
      new G4PVPlacement (rot, G4ThreeVector (lg_x[iLG],lg_y[jLG],0), LguideLV, sname.str().c_str(), Interface_LV, false, 0, checkOverlaps);
      sname.str("");
      G4VisAttributes* LG_VisAtt = new G4VisAttributes(green);  // color
      LG_VisAtt->SetVisibility(true);
      LG_VisAtt->SetForceWireframe(wireFrame);
      LguideLV->SetVisAttributes(LG_VisAtt);
    }
  }
}


void DetectorConstruction::ConstructPMTs(G4LogicalVolume *Readout_LV, std::string moduleName,G4RotationMatrix *rot)
{
  std::stringstream sname;
  G4double PMT_radius = 5*mm;
  G4double PMT_length = ReadoutSizeZ;
  G4VSolid * TubeS = new G4Tubs ("TubeS", 0., PMT_radius, 0.5*PMT_length, 0.*deg, 360.*deg);
  G4LogicalVolume * PMTLV = new G4LogicalVolume (TubeS,PLEXMaterial, "PMTLV");

  int nLGx = 3;
  int nLGy = 3;
  G4double lg_pitch = 20 * mm;
  G4double lg_x[3] = {-lg_pitch,0,lg_pitch};
  G4double lg_y[3] = {-lg_pitch,0,lg_pitch};

  for(int iLG = 0; iLG < nLGx; iLG++)
  {
    for(int jLG = 0; jLG < nLGy; jLG++)
    {
      int pmt_number = iLG*nLGy + jLG;
      sname << "PMT_"
            << pmt_number
            << "_"
            << moduleName;
      new G4PVPlacement (rot, G4ThreeVector (lg_x[iLG],lg_y[jLG],0), PMTLV, sname.str().c_str(), Readout_LV, false, 0, checkOverlaps);
      sname.str("");
      G4VisAttributes* LG_VisAtt = new G4VisAttributes(grey);  // color
      LG_VisAtt->SetVisibility(true);
      LG_VisAtt->SetForceWireframe(wireFrame);
      PMTLV->SetVisAttributes(LG_VisAtt);
    }
  }

}
