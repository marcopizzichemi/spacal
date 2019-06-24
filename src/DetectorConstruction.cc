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
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------

  ConfigFile config (configFileName) ;

  config.readInto (checkOverlaps, "checkOverlaps") ;
  config.readInto (world_material, "world_material") ;
  // config.readInto (abs_material, "abs_material") ;
  config.readInto (W_fraction, "W_fraction") ;
  // config.readInto (hole_radius, "hole_radius") ;
  // config.readInto (module_xy, "module_xy") ;
  // config.readInto (module_yx, "module_yx") ;
  // config.readInto (fibres_x, "fibres_x");
  // config.readInto (fibres_y, "fibres_y");
  // config.readInto (fibres_x1, "fibres_x1");
  // config.readInto (fibres_y1, "fibres_y1");
  //config.readInto (absorber_x, "absorber_x");
  //  config.readInto (absorber_y, "absorber_y");
  // config.readInto (postshower, "postshower") ;
  // config.readInto (second, "second") ;
  //  config.readInto (secondAbs, "secondAbs") ;
  // config.readInto (Second_abs_material, "Second_abs_material") ;
  // config.readInto (Second_module_xy, "Second_module_xy");
  // config.readInto (Second_module_yx, "Second_module_yx") ;
  // config.readInto (Second_W_fraction, "Second_W_fraction") ;
  // config.readInto (Second_fibres_x, "Second_fibres_x");
  // config.readInto (Second_fibres_y, "Second_fibres_y");
  // config.readInto (Second_fibres_x1, "Second_fibres_x1");
  // config.readInto (Second_fibres_y1, "Second_fibres_y1");
  //
  // config.readInto (fibre_scheme, "fibre_scheme") ;
  // config.readInto (Second_fibre_scheme, "Second_fibre_scheme") ;
  // config.readInto (fibre_material, "fibre_material") ;
  // config.readInto (fibre_material1, "fibre_material1") ;
  // config.readInto (fibre_material2, "fibre_material2") ;
  // config.readInto (Second_fibre_material, "Second_fibre_material") ;
  // config.readInto (Second_fibre_material1, "Second_fibre_material1") ;
  // config.readInto (Second_fibre_material2, "Second_fibre_material2") ;
  //
  // config.readInto (fibre_cladRIndex, "fibre_cladRIndex") ;
  // config.readInto (fibre_isSquare, "fibre_isSquare") ;
  // config.readInto (fibre_radius, "fibre_radius") ;
  // config.readInto (fibre_length, "fibre_length") ;
  // config.readInto (Second_fibre_radius, "Second_fibre_radius") ;
  // config.readInto (Second_fibre_length, "Second_fibre_length") ;
  // config.readInto (fibre_distance, "fibre_distance") ;
  // config.readInto (Second_fibre_distance, "Second_fibre_distance") ;
  // config.readInto (fibre_absLength, "fibre_absLength") ;
  //
  // config.readInto (gap_l, "gap_l") ;
  // config.readInto (gap_material, "gap_material") ;
  //
  // config.readInto (det_l, "det_l") ;
  // config.readInto (det_material, "det_material") ;
  //
  // config.readInto (depth, "depth") ;
  //
  // config.readInto (lead_plane, "lead_plane") ;
  // config.readInto (lp_dist, "lp_dist") ;
  // config.readInto (lp_depth, "lp_depth") ;
  // config.readInto (lp_mat, "lp_mat") ;
  // config.readInto (lp_x, "lp_x") ;
  // config.readInto (lp_y, "lp_y") ;
  // config.readInto (preconstr, "preconstr") ;


  config.readInto (surface_lg,"surface_lg");
  config.readInto (glue_interface,"glue_interface");
  config.readInto (cone_material,"cone_material");

  B_field_intensity = config.read<double>("B_field_intensity") * tesla ;


  // margin = std::max( 0.25*fibre_distance, 2.*fibre_radius );
  // G4double staggering = 0.5*fibre_distance*((fibre_scheme+1)%2);
  // G4double staggeredStart = 0.5 + 0.5*(fibre_scheme%2);
  //
  // std::cout << "staggeredStart: " << staggeredStart << std::endl;
  // if( fibre_scheme == 1 || fibre_scheme == 2 ) // dice-4
  // {
  //   fibreDistanceAlongX = fibre_distance;
  //   fibreDistanceAlongY = fibre_distance;
  //   nFibresAlongX = floor( (fibres_x - 2.*margin) / fibreDistanceAlongX ) + 1 ;
  //   nFibresAlongY = floor( (fibres_y - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
  //   nFibresAlongX1 = floor( (fibres_x1 - 2.*margin) / fibreDistanceAlongX ) + 1 ;
  //   nFibresAlongY1 = floor( (fibres_y1 - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
  //   startX = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX - 1.0) ) ;
  //   startY = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY - staggeredStart)) ;
  //   startX1 = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX1 - 1.0) ) ;
  //   startY1 = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY1 - staggeredStart)) ;
  // }
  // else if( fibre_scheme == 3 || fibre_scheme == 4 ) // dice-5
  // {
  //   fibreDistanceAlongX = 0.8660 * fibre_distance;
  //   fibreDistanceAlongY = fibre_distance;
  //   nFibresAlongX = floor( (fibres_x - 2.*margin)/ fibreDistanceAlongX ) + 1 ;
  //   nFibresAlongY = floor( (fibres_y - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
  //   nFibresAlongX1 = floor( (fibres_x1 - 2.*margin) / fibreDistanceAlongX ) + 1 ;
  //   nFibresAlongY1 = floor( (fibres_y1 - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
  //   startX = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX - 1.0) ) ;
  //   startY = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY - staggeredStart)) ;
  //   startX1 = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX1 - 1.0) ) ;
  //   startY1 = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY1 - staggeredStart)) ;
  // }
  // else if( fibre_scheme == 5 || fibre_scheme == 6 ) // chessboard
  // {
  //   fibreDistanceAlongX = 0.5 * fibre_distance;
  //   fibreDistanceAlongY = fibre_distance;
  //   nFibresAlongX = floor( (fibres_x - 2.*margin)              / fibreDistanceAlongX ) + 1 ;
  //   nFibresAlongY = floor( (fibres_y - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
  //   nFibresAlongX1 = floor( (fibres_x1 - 2.*margin)              / fibreDistanceAlongX ) + 1 ;
  //   nFibresAlongY1 = floor( (fibres_y1 - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
  //   startX = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX - 1.0) ) ;
  //   startY = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY - staggeredStart)) ;
  //   startX1 = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX1 - 1.0) ) ;
  //   startY1 = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY1 - staggeredStart)) ;
  // }
  //
  // //________________________________________________________________________________
  //
  //
  // margin2 = std::max( 0.25*Second_fibre_distance, 2.*Second_fibre_radius );
  // G4double staggering2 = 0.5*Second_fibre_distance*((Second_fibre_scheme+1)%2);
  // G4double staggeredStart2 = 0.5 + 0.5*(Second_fibre_scheme%2);
  //
  // std::cout << "staggeredStart: " << staggeredStart2 << std::endl;
  // if( Second_fibre_scheme == 1 || Second_fibre_scheme == 2 ) // dice-4
  // {
  //   Second_fibreDistanceAlongX = Second_fibre_distance;
  //   Second_fibreDistanceAlongY = Second_fibre_distance;
  //   Second_nFibresAlongX = floor( (Second_fibres_x - 2.*margin2) / Second_fibreDistanceAlongX ) + 1 ;
  //   Second_nFibresAlongY = floor( (Second_fibres_y - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
  //   Second_nFibresAlongX1 = floor( (Second_fibres_x1 - 2.*margin2) / Second_fibreDistanceAlongX ) + 1 ;
  //   Second_nFibresAlongY1 = floor( (Second_fibres_y1 - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
  //   Second_startX = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX - 1.0) ) ;
  //   Second_startY = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY - staggeredStart2)) ;
  //   Second_startX1 = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX1 - 1.0) ) ;
  //   Second_startY1 = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY1 - staggeredStart)) ;
  // }
  // else if( Second_fibre_scheme == 3 || Second_fibre_scheme == 4 ) // dice-5
  // {
  //   Second_fibreDistanceAlongX = 0.8660 * Second_fibre_distance;
  //   Second_fibreDistanceAlongY = Second_fibre_distance;
  //   Second_nFibresAlongX = floor( (Second_fibres_x - 2.*margin2)/ Second_fibreDistanceAlongX ) + 1 ;
  //   Second_nFibresAlongY = floor( (Second_fibres_y - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
  //   Second_nFibresAlongX1 = floor( (Second_fibres_x1 - 2.*margin2) / Second_fibreDistanceAlongX ) + 1 ;
  //   Second_nFibresAlongY1 = floor( (Second_fibres_y1 - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
  //   Second_startX = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX - 1.0) ) ;
  //   Second_startY = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY - staggeredStart2)) ;
  //   Second_startX1 = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX1 - 1.0) ) ;
  //   Second_startY1 = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY1 - staggeredStart)) ;
  //
  // }
  // else if( Second_fibre_scheme == 5 || Second_fibre_scheme == 6 ) // chessboard
  // {
  //   Second_fibreDistanceAlongX = 0.5 * Second_fibre_distance;
  //   Second_fibreDistanceAlongY = Second_fibre_distance;
  //   Second_nFibresAlongX = floor( (Second_fibres_x - 2.*margin2)              / Second_fibreDistanceAlongX ) + 1 ;
  //   Second_nFibresAlongY = floor( (Second_fibres_y - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
  //   Second_nFibresAlongX1 = floor( (Second_fibres_x1 - 2.*margin2)              / Second_fibreDistanceAlongX ) + 1 ;
  //   Second_nFibresAlongY1 = floor( (Second_fibres_y1 - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
  //   Second_startX = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX - 1.0) ) ;
  //   Second_startY = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY - staggeredStart2)) ;
  //   Second_startX1 = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX1 - 1.0) ) ;
  //   Second_startY1 = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY1 - staggeredStart)) ;
  // }

  //---------------------------------------//
  // CALORIMETER                           //
  //---------------------------------------//

  config.readInto(modules_nx,"modules_nx");
  config.readInto(modules_ny,"modules_ny");


  //---------------------------------------//
  // ABSORBER                              //
  //---------------------------------------//
  // read from config file


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

  // Absorber absorber;
  absorber.SetName(AbsName);
  // absorber.SetID(i+1);
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

  //visibility
  config.readInto(visibility,"crystalsVisibility");
  config.readInto(wireFrame,"wireFrame");


  // create cells
  // for(unsigned int iAbs = 0 ; iAbs < AbsNames.size(); iAbs++)
  // {
    for(unsigned int iCell = 0 ; iCell < CellNames.size(); iCell++)
    {
      Cell cell;
      cell.SetID(iCell);
      cell.SetName(CellNames[iCell]);
      // cell.SetDimensions (CellSizeX[iCell],
      //                     CellSizeY[iCell],
      //                     CellSizeZ[iCell]);
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
      // cell.
      cell.MakeCellStruture();
      absorber.AddCell(cell);
    }

    G4cout << "Absorber " << absorber.GetName()
          <<  " cells = " << absorber.GetNumberOfCells() << G4endl;
  // }

  /////////_____________________________________________________________________________





  //

  gap_abs_interface = config.read<int>("gap_abs_interface",0);
  gap_interface_readout = config.read<int>("gap_interface_readout",0);

  InterfaceSizeX = AbsSizeX;
  InterfaceSizeY = AbsSizeY;
  InterfaceSizeZ = 5*cm; //TEMP
  ReadoutSizeX = AbsSizeX;
  ReadoutSizeY = AbsSizeY;
  ReadoutSizeZ = 10*cm; //TEMP



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

  //-----------------------------------------------------
  //------------- Define colors --------------
  //-----------------------------------------------------
  G4Colour  white   (1.00, 1.00, 1.00) ;  // white
  G4Colour  grey    (0.50, 0.50, 0.50) ;  // grey
  G4Colour  black   (0.00, 0.00, 0.00) ;  // black
  G4Colour  red     (1.00, 0.00, 0.00) ;  // red
  G4Colour  green   (0.00, 1.00, 0.00) ;  // green
  G4Colour  blue    (0.00, 0.00, 1.00) ;  // blue
  G4Colour  cyan    (0.00, 1.00, 1.00) ;  // cyan
  G4Colour  air     (0.00, 1.00, 1.00) ;  // cyan
  G4Colour  magenta (1.00, 0.00, 1.00) ;  // magenta
  G4Colour  yellow  (1.00, 1.00, 0.00) ;  // yellow
  G4Colour  brass   (0.80, 0.60, 0.40) ;  // brass
  G4Colour  brown   (0.70, 0.40, 0.10) ;  // brown
  G4Colour  orange  (1.00, 0.33, 0.00) ;  // orange



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









  // }
  // else
  // {
    // new G4PVPlacement (0, G4ThreeVector (0., 0., 0.), calorimeterLV, "calorimeterPV", worldLV, false, 0, checkOverlaps) ;
    // new G4PVPlacement (0, G4ThreeVector (0., 0., 1.*(module_z + Second_module_z)), calorimeterLV, "postshower1PV", worldLV, false, 0, checkOverlaps) ;
    // new G4PVPlacement (0, G4ThreeVector (0., 0., 2.*(module_z + Second_module_z)), calorimeterLV, "postshower2PV", worldLV, false, 0, checkOverlaps) ;
  // }

  // Plane in front of Calo
  // if ( lead_plane == 1)
  // {
  //
  //   G4VSolid * PLANE = new G4Box ("PLANE", 0.5*lp_x, 0.5*lp_y, 0.5*lp_depth);
  //   G4LogicalVolume * PLANELV = new G4LogicalVolume (PLANE, PlaneMaterial, "PLANELV");
  //   new G4PVPlacement (0, G4ThreeVector (0.,0.,-0.5*(lp_dist + lp_depth)),PLANELV,"PLANEPV", worldLV, false, 0 ,checkOverlaps);
  // }

  // // Light guide structure
  // // Truncated of cone, cut on the sides
  // // goes from 20x20 mm2 (on the fibers) to 5 mm radius on the PMT
  // // first, do a simple truncated cone
  // G4VSolid* coneS = new G4Cons("aCone",0,PMT_radius,0, 0.5*sqrt(2)*lguide_edge ,0.5*PLEX_depth,0.0 * deg,360.0* deg);
  // // then prepare a hollow box, section of the hole 20x20 mm2
  // G4Box* innerBoxS = new G4Box ("innerBoxS", 0.5*lguide_edge, 0.5*lguide_edge, 0.5*PLEX_depth) ;
  // G4Box* outerBoxS = new G4Box ("outerBoxS", 0.5*sqrt(2)*lguide_edge, 0.5*sqrt(2)*lguide_edge, 0.5*PLEX_depth) ;
  // G4VSolid* subtract = new G4SubtractionSolid("Hollow-Box", outerBoxS, innerBoxS,0,G4ThreeVector(0.,0.,0.));
  // // subtract hollow box from truncated cone
  // G4VSolid* coneSolid = new G4SubtractionSolid("coneSolid", coneS, subtract,0,  G4ThreeVector(0.,0.,0.));
  // G4LogicalVolume * LguideLV = new G4LogicalVolume (coneSolid,ConeMaterial,"LguideLV");
  //
  // //air light guide
  // G4VSolid * air_LG_front_S = new G4Box ("air_LG_front_S ",0.5*module_xy, 0.5*module_yx,0.5*PLEX_depth);
  // G4VSolid * air_LG_back_S = new G4Box ("air_LG_back_S ",0.5*module_xy, 0.5*module_yx,0.5*PLEX_depth);
  // G4LogicalVolume *air_LG_front_LV = new G4LogicalVolume (air_LG_front_S,WoMaterial, "air_LG_front_LV");
  // G4LogicalVolume *air_LG_back_LV = new G4LogicalVolume (air_LG_back_S,WoMaterial, "air_LG_back_LV");




  // G4LogicalVolume * coneLV = new G4LogicalVolume (coneSolid,PLEXMaterial,"coneLV");
  // G4PVPlacement *conePV = new G4PVPlacement (0,G4ThreeVector(0,0,-30.0 *cm), coneLV, "conePV", worldLV, false, 0, checkOverlaps) ;


  // G4VSolid * PVCS = new G4Box ("PVCS", 0.5*module_xy, 0.5*module_yx, 0.5*PVC_depth);
  // G4LogicalVolume * PVCLV = new G4LogicalVolume (PVCS,PVCMaterial, "PVCLV");
  // G4VSolid * TubeS = new G4Tubs ("TubeS", 0., PMT_radius, 0.5*PMT_length, 0.*deg, 360.*deg);
  //
  // G4VSolid * WiresS = new G4Tubs ("WiresS", 0., Wires_radius, Wires_dist, 0.*deg, 360.*deg);
  //
  // G4LogicalVolume * PMTLV = new G4LogicalVolume (TubeS,PLEXMaterial, "PMTLV");
  // G4LogicalVolume * WiresLV = new G4LogicalVolume (WiresS, WiresMaterial, "WiresLV");
  //
  // G4VSolid * PVC_pmt_frontS = new G4Box ("PVC_pmt_frontS",0.5*module_xy, 0.5*module_yx,0.5*PMT_length);
  // G4VSolid * PVC_pmt_backS  = new G4Box ("PVC_pmt_backS" ,0.5*module_xy, 0.5*module_yx,0.5*PMT_length);
  // G4LogicalVolume * PVC_pmt_frontLV = new G4LogicalVolume (PVC_pmt_frontS,PVCMaterial, "PVC_pmt_frontLV");
  // G4LogicalVolume * PVC_pmt_backLV  = new G4LogicalVolume (PVC_pmt_backS,PVCMaterial, "PVC_pmt_frontLV");
  //
  //
  //
  // G4PVPlacement *PVC_pmt_frontPV = new G4PVPlacement (0, G4ThreeVector (0, 0,-0.5*(PLEX_dist+PMT_length+PLEX_depth*2+airGap*4 + airGap_abs) ),PVC_pmt_frontLV, "PVC_pmt_frontPV", worldLV, false, 0, checkOverlaps);
  // G4PVPlacement *PVC_pmt_backPV = new G4PVPlacement (0, G4ThreeVector (0, 0,+0.5*(PLEX_dist+PMT_length+PLEX_depth*2+airGap*4 + airGap_abs) ),PVC_pmt_backLV, "PVC_pmt_backPV", worldLV, false, 0, checkOverlaps);
  //
  // G4PVPlacement *air_LG_front_PV = new G4PVPlacement(0, G4ThreeVector(0 ,0 , -0.5*(PLEX_dist + PLEX_depth + airGap*2 + airGap_abs)), air_LG_front_LV, "air_LG_front_PV", worldLV, false, 0, checkOverlaps);
  // G4PVPlacement  *air_LG_back_PV = new G4PVPlacement(0, G4ThreeVector(0 ,0 , +0.5*(PLEX_dist + PLEX_depth + airGap*2 + airGap_abs)), air_LG_back_LV , "air_LG_back_PV" , worldLV, false, 0, checkOverlaps);
  //
  // G4VSolid *glue_S = new G4Box ("glue_S", 0.5*module_xy, 0.5*module_yx, 0.5*airGap);
  // G4LogicalVolume * glue_LV = new G4LogicalVolume (glue_S,GlueMaterial, "glue_LV");
  //
  // G4PVPlacement  *gluePV_front = new G4PVPlacement(0, G4ThreeVector(0 ,0 , -0.5*(PLEX_dist + airGap + airGap_abs )), glue_LV, "gluePV_front", worldLV, false, 0, checkOverlaps);
  // G4PVPlacement  *gluePV_back  = new G4PVPlacement(0, G4ThreeVector(0 ,0 , +0.5*(PLEX_dist + airGap + airGap_abs )), glue_LV, "gluePV_back" , worldLV, false, 0, checkOverlaps);
  //
  // G4PVPlacement  *glue_pmt_PV_front = new G4PVPlacement(0, G4ThreeVector(0 ,0 , -0.5*(PLEX_dist + 3*airGap + airGap_abs +PLEX_depth*2)), glue_LV, "glue_pmt_PV_front", worldLV, false, 0, checkOverlaps);
  // G4PVPlacement  *glue_pmt_PV_back  = new G4PVPlacement(0, G4ThreeVector(0 ,0 , +0.5*(PLEX_dist + 3*airGap + airGap_abs +PLEX_depth*2)), glue_LV, "glue_pmt_PV_back" , worldLV, false, 0, checkOverlaps);





  // if (preconstr == 1)
  // {
  //   // The Pre-Construction
  //   // mods for real spacal test beam structure
  //   //----------- ---------------
  //
  //   // new G4PVPlacement (0, G4ThreeVector (0.,0.,-0.5*(PVC_dist+PVC_depth)),PVCLV,"PVCPV", worldLV, false, 0, checkOverlaps);
  //
  //   // new G4PVPlacement (0, G4ThreeVector (0.,0.,-(PLEX_dist+PVC_depth+60)),PVCLV,"PVCPV", worldLV, false, 0, checkOverlaps);
  //
  //   startAX = 0. + 0.5*(absorber_x);
  //   nCellsAlongX =  floor( module_xy / absorber_x );
  //   nCellsAlongY  = floor ( module_yx / absorber_y);
  //   startAY = 0.+ 0.5*(absorber_y + (module_yx-module_xy));
  //
  //   // The PMTS
  //   // ------------------------------------- --------  ------
  //
  //
  //   int NX = 1;
  //
  //
  //   float pmt_x[3] = {-pmts_pitch,0,pmts_pitch};
  //   float pmt_y[3] = {-pmts_pitch,0,pmts_pitch};
  //   for(int iPMT = 0 ; iPMT < 3; iPMT++)
  //   {
  //     for(int jPMT = 0 ; jPMT < 3; jPMT++)
  //     {
  //       float PMTx = pmt_x[iPMT];
  //       float PMTy = pmt_y[jPMT];
  //
  //       int iP = 3* iPMT + jPMT;
  //       std::string Pname;
  //
  //       //-----------//
  //       // FRONT     //
  //       //-----------//
  //
  //       Pname = Form("PMTSPV FRONT %d", iP);
  //       new G4PVPlacement (0, G4ThreeVector (PMTx, PMTy,0 ),PMTLV, Pname, PVC_pmt_frontLV, false, 0, checkOverlaps);
  //
  //       Pname = Form("LguidePV FRONT %d", iP);
  //       G4PVPlacement *LG_front_PV = new G4PVPlacement (0, G4ThreeVector (PMTx,PMTy,0), LguideLV, Pname, air_LG_front_LV, false, 0, checkOverlaps);
  //
  //       //optical surface between cone and air cone
  //       if(surface_lg)
  //       {
  //         std::stringstream Surfname_front;
  //         Surfname_front << "Surface_" << Pname << "_air_front";
  //         G4OpticalSurface* reflector_surf_front = new G4OpticalSurface(Surfname_front.str().c_str());
  //         reflector_surf_front->SetType(dielectric_metal);
  //         reflector_surf_front->SetFinish(polished);
  //         reflector_surf_front->SetModel(unified);
  //         reflector_surf_front->SetMaterialPropertiesTable(MyMaterials::ESR());
  //         G4LogicalBorderSurface* reflectorLB_front;
  //         Surfname_front << "LB";
  //         reflectorLB_front = new G4LogicalBorderSurface(Surfname_front.str().c_str(),
  //         LG_front_PV,
  //         air_LG_front_PV,
  //         reflector_surf_front);
  //         //optical surface between air cone and cone
  //         std::stringstream SurfnameInv_front;
  //         SurfnameInv_front << "SurfaceInv_" << Pname << "_air_front";
  //         G4OpticalSurface* reflectorInv_surf_front = new G4OpticalSurface(SurfnameInv_front.str().c_str());
  //         reflectorInv_surf_front->SetType(dielectric_metal);
  //         reflectorInv_surf_front->SetFinish(polished);
  //         reflectorInv_surf_front->SetModel(unified);
  //         reflectorInv_surf_front->SetMaterialPropertiesTable(MyMaterials::ESR());
  //         G4LogicalBorderSurface* reflectorLBInv_front;
  //         SurfnameInv_front << "LB";
  //         reflectorLBInv_front = new G4LogicalBorderSurface(SurfnameInv_front.str().c_str(),
  //         air_LG_front_PV,
  //         LG_front_PV,
  //         reflectorInv_surf_front);
  //       }
  //
  //       //-----------//
  //       // BACK      //
  //       //-----------//
  //
  //       G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
  //       rotationMatrix->rotateY(180.*deg);
  //
  //       Pname = Form("PMTSPV BACK %d", iP);
  //       new G4PVPlacement (0, G4ThreeVector (PMTx, PMTy,0 ),PMTLV, Pname, PVC_pmt_backLV, false, 0, checkOverlaps);
  //
  //       Pname = Form("LguidePV BACK %d", iP);
  //       G4PVPlacement *LG_back_PV = new G4PVPlacement (rotationMatrix, G4ThreeVector (PMTx,PMTy, 0) , LguideLV, Pname, air_LG_back_LV, false, 0, checkOverlaps);
  //
  //
  //       if(surface_lg)
  //       {
  //         //optical surface between cone and air cone
  //         std::stringstream Surfname_back;
  //         Surfname_back << "Surface_" << Pname << "_air_back";
  //         G4OpticalSurface* reflector_surf_back = new G4OpticalSurface(Surfname_back.str().c_str());
  //         reflector_surf_back->SetType(dielectric_metal);
  //         reflector_surf_back->SetFinish(polished);
  //         reflector_surf_back->SetModel(unified);
  //         reflector_surf_back->SetMaterialPropertiesTable(MyMaterials::ESR());
  //         G4LogicalBorderSurface* reflectorLB_back;
  //         Surfname_back << "LB";
  //         reflectorLB_back = new G4LogicalBorderSurface(Surfname_back.str().c_str(),
  //         LG_back_PV,
  //         air_LG_back_PV,
  //         reflector_surf_back);
  //         //optical surface between air cone and cone
  //         std::stringstream SurfnameInv_back;
  //         SurfnameInv_back << "SurfaceInv_" << Pname << "_air_back";
  //         G4OpticalSurface* reflectorInv_surf_back = new G4OpticalSurface(SurfnameInv_back.str().c_str());
  //         reflectorInv_surf_back->SetType(dielectric_metal);
  //         reflectorInv_surf_back->SetFinish(polished);
  //         reflectorInv_surf_back->SetModel(unified);
  //         reflectorInv_surf_back->SetMaterialPropertiesTable(MyMaterials::ESR());
  //         G4LogicalBorderSurface* reflectorLBInv_back;
  //         SurfnameInv_back << "LB";
  //         reflectorLBInv_back = new G4LogicalBorderSurface(SurfnameInv_back.str().c_str(),
  //         air_LG_back_PV,
  //         LG_back_PV,
  //         reflectorInv_surf_back);
  //
  //       }
  //     }
  //   }
  //
  //   // OLD PMT PLACEMENT CODE
  //   // for (float xP = -0.5 * module_xy + module_xy/6; NX < 4; xP += module_xy/3)
  //   // {
  //   //   int NY = 1;
  //   //   for (float yP = -0.5 * module_yx + module_yx/6; NY < 4; yP += module_yx/3)
  //   //   {
  //   //     float PMTx = xP;
  //   //     float PMTy = yP;
  //   //
  //   //     int iP = 3* NX + NY;
  //   //
  //   //     std::string Pname;
  //   //     Pname = Form("PMTSPV %d", iP);
  //   //
  //   // std::string Wname;
  //   // Wname = Form("WiresSPV %d", iP);
  //   //
  //   //     new G4PVPlacement (0, G4ThreeVector (PMTx, PMTy,-0.5*(PLEX_dist+PMT_length+PLEX_depth*2)),PMTLV, Pname, worldLV, false, 0, checkOverlaps);
  //   // new G4PVPlacement (0, G4ThreeVector (PMTx, PMTy,-0.5*(PVC_dist + PVC_depth*2 )),WiresLV, Wname, worldLV, false, 0, checkOverlaps);
  //   //
  //   //     ++NY;
  //   //   }
  //   //   ++NX;
  //   // }
  // }

  // The holes
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  // G4VSolid * holeS;
  // if( !fibre_isSquare ) holeS = new G4Tubs ("holeS", fibre_radius, fibre_radius+hole_radius, 0.5*fibre_length, 0.*deg, 360.*deg) ;
  // else
  // {
  //   G4VSolid * temp1 = new G4Box ("temp1", fibre_radius+hole_radius, fibre_radius+hole_radius, 0.5*fibre_length) ;
  //   G4VSolid * temp2 = new G4Box ("temp2", fibre_radius, fibre_radius, 1.6*fibre_length) ;
  //   holeS = new G4SubtractionSolid("holeS",temp1,temp2,0,G4ThreeVector(0.,0.,0.));
  // }
  // G4LogicalVolume * holeLV = new G4LogicalVolume (holeS, WoMaterial, "holeLV") ;
  // //HoleParameterisation* holeParam = new HoleParameterisation(module_xy,fibre_radius,hole_radius,fibre_distance,fibre_length,WoMaterial);
  // //new G4PVParameterised("holeP", holeLV, absorberLV, kUndefined, holeParam->GetNHoles(), holeParam);
  //
  //
  // // the fibres
  // // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  // G4VSolid * fibre1S;
  // if( !fibre_isSquare ) fibre1S = new G4Tubs ("fibre1S", 0., fibre_radius, 0.5*fibre_length, 0.*deg, 360.*deg) ;
  // else                  fibre1S = new G4Box ("fibre1S", fibre_radius, fibre_radius, 0.5*fibre_length) ;
  // G4LogicalVolume * fibre1SLV = new G4LogicalVolume (fibre1S, ClMaterial, "fibre1LV") ;
  // //FibreParameterisation* fibreParam = new FibreParameterisation(module_xy,fibre_radius,fibre_distance,fibre_length,ClMaterial);
  // //new G4PVParameterised("fibreP", fibreLV, absorberLV, kUndefined, fibreParam->GetNFibres(), fibreParam);
  //
  //
  // //___________________________________________________________
  //
  //
  // G4VSolid * fibre12S;
  // if( !fibre_isSquare ) fibre12S = new G4Tubs ("fibre12S", 0., fibre_radius, 0.5*fibre_length, 0.*deg, 360.*deg) ;
  // else                  fibre12S = new G4Box ("fibre12S", fibre_radius, fibre_radius, 0.5*fibre_length) ;
  // G4LogicalVolume * fibre12SLV = new G4LogicalVolume (fibre12S, ClSSMaterial, "fibre12SLV") ;
  //
  // G4VSolid * fibre13S;
  // if( !fibre_isSquare ) fibre13S = new G4Tubs ("fibre13S", 0., fibre_radius, 0.5*fibre_length, 0.*deg, 360.*deg) ;
  // else                  fibre13S = new G4Box ("fibre13S", fibre_radius, fibre_radius, 0.5*fibre_length) ;
  // G4LogicalVolume * fibre13SLV = new G4LogicalVolume (fibre13S, ClSSSMaterial , "fibre13SLV") ;
  //
  //
  // G4VSolid * fibre2S;
  // if( !fibre_isSquare ) fibre2S = new G4Tubs ("fibre2S", 0., fibre_radius, 0.5*Second_fibre_length, 0.*deg, 360.*deg) ;
  // else                  fibre2S = new G4Box ("fibre2S", Second_fibre_radius, Second_fibre_radius, 0.5*Second_fibre_length) ;
  // G4LogicalVolume * fibre2SLV = new G4LogicalVolume (fibre2S, Cl3Material, "fibre2LV") ;
  // //FibreParameterisation* fibreParam = new FibreParameterisation(module_xy,fibre_radius,fibre_distance,fibre_length,ClMaterial);
  // //new G4PVParameterised("fibreP", fibreLV, absorberLV, kUndefined, fibreParam->GetNFibres(), fibreParam);
  //
  //
  // //___________________________________________________________
  //
  //
  // G4VSolid * fibre22S;
  // if( !fibre_isSquare ) fibre22S = new G4Tubs ("fibre22S", 0., Second_fibre_radius, 0.5*Second_fibre_length, 0.*deg, 360.*deg) ;
  // else                  fibre22S = new G4Box ("fibre22S", Second_fibre_radius, Second_fibre_radius, 0.5*Second_fibre_length) ;
  // G4LogicalVolume * fibre22SLV = new G4LogicalVolume (fibre22S, Cl4SSMaterial, "fibre2SSLV") ;
  //
  //
  // G4VSolid * fibre23S;
  // if( !fibre_isSquare ) fibre23S = new G4Tubs ("fibre23S", 0., Second_fibre_radius, 0.5*Second_fibre_length, 0.*deg, 360.*deg) ;
  // else                  fibre23S = new G4Box ("fibre23S", Second_fibre_radius, Second_fibre_radius, 0.5*Second_fibre_length) ;
  // G4LogicalVolume * fibre23SLV = new G4LogicalVolume (fibre23S, Cl43SSMaterial , "fibre23SLV") ;
  //
  // //_________________________________________________________________
  //

  // G4cout << "Sections: " << nCellsAlongX << "   " << nCellsAlongY << G4endl;



  //thin layer of air between absorbers
  // G4VSolid * airLayer = new G4Box ("airLayer", 0.5 * Second_module_xy, 0.5 * Second_module_yx, 0.5 * airGap_abs ) ;
  // G4LogicalVolume * airLayerLV;
  // airLayerLV = new G4LogicalVolume (airLayer, WoMaterial, "airLayerLV") ;
  // G4PVPlacement* airLayerPV = new G4PVPlacement (0, G4ThreeVector (0., 0., 0.), airLayerLV, "airLayerPV", calorimeterLV, false, 0, checkOverlaps) ;








  // //------------------------------------------------------------------//
  // // OLD PLACEMENT STYLE (commented out on 13/06/2019)                //
  // //------------------------------------------------------------------//
  //
  // // The Absorbers
  // // First Absorber
  // G4VSolid * absorber1S = new G4Box ("absorber1S", 0.5 * module_xy, 0.5 * (module_yx), 0.5 * module_z) ;
  // G4LogicalVolume * absorber1LV;
  // absorber1LV = new G4LogicalVolume (absorber1S, AbMaterial, "absorber1LV") ;
  // G4PVPlacement *absorber1PV =  new G4PVPlacement (0, G4ThreeVector (0., 0., - 0.5* (module_z + airGap_abs)), absorber1LV, "absorber1PV", calorimeterLV, false, 0, checkOverlaps) ;
  //
  // // Second Absorber
  // G4VSolid * absorber2S = new G4Box ("absorber2S", 0.5 * Second_module_xy, 0.5 * Second_module_yx, 0.5 * Second_module_z) ;
  // G4LogicalVolume * absorber2LV;
  // absorber2LV = new G4LogicalVolume (absorber2S, AbMaterial2, "absorber2LV") ;
  // G4PVPlacement *absorber2PV = new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5* (module_z + airGap_abs)), absorber2LV, "absorber2PV", calorimeterLV, false, 0, checkOverlaps) ;
  //
  // // fibres matrix filling
  // // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  // int index;
  // int indexion;
  // int index1;
  // // loop on x direction
  // int countX = 0 ; // for the staggering
  //
  // for (float x = -0.5*module_xy+startX; countX <nFibresAlongX; x += fibreDistanceAlongX)
  // {
  //   // loop on y direction
  //   int countY = 0 ; // for the staggering
  //   for (float y = - 0.5 * module_yx + startY; countY < nFibresAlongY; y += fibreDistanceAlongY)
  //   {
  //     float x_c = x;
  //     float y_c = y;
  //     //____________________________________________
  //     //   float x_cs = -0.5 * module_xy;
  //     //    float y_cs = -0.5 * module_xy;
  //     //________________________________________________________________
  //     // staggering
  //     if( (fibre_scheme%2) == 0 )
  //     y_c += 0.5*fibreDistanceAlongY*(countX%2) ;
  //
  //     int index = countX * nFibresAlongY + countY ;
  //
  //     CreateTree::Instance() -> fibresPosition -> Fill(index,x_c,y_c);
  //
  //
  //     std::string name;
  //
  //
  //     // **************  ***************
  //     // #5
  //     G4PVPlacement *fiberPV;
  //     G4PVPlacement *holePV;
  //     if (x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1)
  //     {
  //       name = Form("holePV %d",index);
  //       holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber1LV, false, 0, checkOverlaps) ;
  //
  //       name = Form("fibre12SPV %d",index);
  //       fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre12SLV, name, absorber1LV, false, 0, checkOverlaps) ;
  //     }
  //     // #2,#4,#6,#8
  //     else if (x >= -1.5*fibres_x1 && x<= fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1 || x >= fibres_x1 && x<= 1.5*fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1 || x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= -1.5*fibres_y1 && y <= -0.5*fibres_y1 ||  x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= 0.5*fibres_y1 && y <= 1.5*fibres_y1 )
  //     {
  //       name = Form("holePV %d",index);
  //       holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber1LV, false, 0, checkOverlaps) ;
  //
  //       name = Form("fibre1SPV %d",index);
  //       fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre1SLV, name, absorber1LV, false, 0, checkOverlaps) ;
  //     }
  //
  //     else
  //     {
  //       name = Form("holePV %d",index);
  //       holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber1LV, false, 0, checkOverlaps) ;
  //
  //       name = Form("fibre13SPV %d",index);
  //       fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre13SLV, name, absorber1LV, false, 0, checkOverlaps) ;
  //     }
  //
  //     indexion = index;
  //     index1 = index+1;
  //     CreateTree::Instance() -> fibresPosition_1st_Section -> Fill(indexion,x_c,y_c);
  //     ++countY ;
  //
  //
  //
  //     //optical surface between fiber and thin air layer
  //     std::stringstream Surfname;
  //     Surfname << "Surface_" << name << "_airLayer";
  //     G4OpticalSurface* reflector_surf = new G4OpticalSurface(Surfname.str().c_str());
  //     reflector_surf->SetType(dielectric_metal);
  //     reflector_surf->SetFinish(polished);
  //     reflector_surf->SetModel(unified);
  //     reflector_surf->SetMaterialPropertiesTable(MyMaterials::ESR());
  //     G4LogicalBorderSurface* reflectorLB;
  //     Surfname << "LB";
  //     reflectorLB = new G4LogicalBorderSurface(Surfname.str().c_str(),
  //     fiberPV,
  //     airLayerPV,
  //     reflector_surf);
  //     //optical surf from air gap to abs
  //     std::stringstream absName;
  //     absName << "Surface_hole_" << name << "_abs";
  //     G4OpticalSurface* absSurface = new G4OpticalSurface(absName.str().c_str());
  //     absSurface->SetType(dielectric_metal);
  //     absSurface->SetFinish(ground);
  //     absSurface->SetModel(unified);
  //     absSurface->SetSigmaAlpha(0.1);
  //     absSurface->SetMaterialPropertiesTable(MyMaterials::ABS_SURF());
  //     G4LogicalBorderSurface* absLB;
  //     absName << "LB";
  //     absLB = new G4LogicalBorderSurface(absName.str().c_str(),                                                                                       holePV,
  //     absorber1PV,
  //     absSurface);
  //
  //
  //
  //
  //
  //   } // loop on y direction
  //
  //   ++countX ;
  // } // loop on x direction
  //
  //
  //
  //
  //
  //
  // //***********************8**    ******************   8******************
  // // loop on x direction
  // int countX3 = 0 ; // for the staggering
  // int indexion4;
  // int indexion3;
  // for (float x = -0.5*Second_module_xy+Second_startX; countX3 < Second_nFibresAlongX; x += Second_fibreDistanceAlongX)
  // {
  //   // loop on y direction
  //   int countY3 = 0 ; // for the staggering
  //   for (float y = - 0.5 * Second_module_yx + Second_startY; countY3 < Second_nFibresAlongY; y += Second_fibreDistanceAlongY)
  //   {
  //     float x_c = x;
  //     float y_c = y;
  //     //____________________________________________
  //     //   float x_cs = -0.5 * module_xy;
  //     //    float y_cs = -0.5 * module_xy;
  //     //________________________________________________________________
  //     // staggering
  //     if( (Second_fibre_scheme%2) == 0 )
  //     y_c += 0.5*Second_fibreDistanceAlongY*(countX3%2) ;
  //
  //     int index = index1 + countX3 * Second_nFibresAlongY + countY3 ;
  //     int index2 = countX3 * Second_nFibresAlongY + countY3;
  //     CreateTree::Instance() -> fibresPosition -> Fill(index,x_c,y_c);
  //     CreateTree::Instance() -> fibresPosition_2nd_Section -> Fill(index2,x_c,y_c);
  //
  //     std::string name;
  //     G4PVPlacement *fiberPV;
  //     G4PVPlacement *holePV;
  //
  //
  //     // **************  ***************
  //     if (x >= -0.5*Second_fibres_x1 && x<= 0.5*Second_fibres_x1 && y >= -0.5*Second_fibres_y1 && y <= 0.5*Second_fibres_y1)
  //     {
  //       name = Form("holePV %d",index);
  //       holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber2LV, false, 0, checkOverlaps) ;
  //
  //       name = Form("fibre22SPV %d",index);
  //       fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre22SLV, name, absorber2LV, false, 0, checkOverlaps) ;
  //
  //
  //       // ************ *********************
  //     }
  //
  //     // #11,#13,#15,#17
  //     else if (x >= -1.5*fibres_x1 && x<= fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1 || x >= fibres_x1 && x<= 1.5*fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1 || x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= -1.5*fibres_y1 && y <= -0.5*fibres_y1 ||  x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= 0.5*fibres_y1 && y <= 1.5*fibres_y1 )
  //     {
  //       name = Form("holePV %d",index);
  //       holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber2LV, false, 0, checkOverlaps) ;
  //
  //       name = Form("fibre2SPV %d",index);
  //       fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre2SLV, name, absorber2LV, false, 0, checkOverlaps) ;
  //     }
  //
  //     else
  //     {
  //       name = Form("holePV %d",index);
  //       holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber2LV, false, 0, checkOverlaps) ;
  //
  //       name = Form("fibre23SPV %d",index);
  //       fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre23SLV, name, absorber2LV, false, 0, checkOverlaps) ;
  //     }
  //
  //     indexion4 = index+1;
  //     indexion3 = index2+1;
  //     ++countY3 ;
  //
  //     //optical surface between fiber and thin air
  //     std::stringstream Surfname;
  //     Surfname << "Surface_" << name << "_airLayer";
  //     G4OpticalSurface* reflector_surf = new G4OpticalSurface(Surfname.str().c_str());
  //     reflector_surf->SetType(dielectric_metal);
  //     reflector_surf->SetFinish(polished);
  //     reflector_surf->SetModel(unified);
  //     reflector_surf->SetMaterialPropertiesTable(MyMaterials::ESR());
  //     G4LogicalBorderSurface* reflectorLB;
  //     Surfname << "LB";
  //     reflectorLB = new G4LogicalBorderSurface(Surfname.str().c_str(),
  //     fiberPV,
  //     airLayerPV,
  //     reflector_surf);
  //     //
  //     std::stringstream absName;
  //     absName << "Surface_hole_" << name << "_abs";
  //     G4OpticalSurface* absSurface = new G4OpticalSurface(absName.str().c_str());
  //     absSurface->SetType(dielectric_metal);
  //     absSurface->SetFinish(ground);
  //     absSurface->SetModel(unified);
  //     absSurface->SetSigmaAlpha(0.1);
  //     absSurface->SetMaterialPropertiesTable(MyMaterials::ABS_SURF());
  //     G4LogicalBorderSurface* absLB;
  //     absName << "LB";
  //     absLB = new G4LogicalBorderSurface(absName.str().c_str(),                                                                             holePV,
  //     absorber2PV,
  //     absSurface);
  //   } // loop on y direction
  //   ++countX3 ;
  // } // loop on x direction
  // //------------------------------------------------------------------//
  // // end of OLD PLACEMENT STYLE (commented out on 13/06/2019)         //
  // //------------------------------------------------------------------//












  //------------------------------------------------------------------//
  // NEW PLACEMENT STYLE (started on 13/06/2019)                      //
  //------------------------------------------------------------------//
  //
  // New placement style is intended to reflect more accurately the SPACAL-RD test beam module
  // Old placement had (slightly) wrong fiber positioning.
  //
  // First, a bit of explanation on directions:
  // 1) in this simulation, the spacal longitutidal axis is oriented like z
  // 2) so the transverse plan is x-y.
  // 3) The beam of particle has always been shot from negative z (at least according to
  //    gps.mac files found in the code)
  // 4) The simulation calls FRONT the absorber placed in the negative z part. Its name is absorber1 and
  //    it is visualized in grey. So the particles are always entering the module from the FRONT
  // 5) The rows of crystals, or better the rows of absorber that are stacked on top of one another
  //    to form the module, are on the z-y plane, and are piled up along x.
  //
  // So, watching the module from the FRONT part, the z axis is pointing away from the viewer, and
  // the other two axis are placed like
  //
  //
  //     x
  //     ^
  //     |
  //     |
  //     |
  //     |
  //     ------------> y
  //
  //
  // From this viewpoint, the module looked like a pile up of rows with alternate number of elements.
  // First row from the bottom is made of 33 holes (we call it type A), second row 32 holes (type B),
  // third type A and so on...
  // This is the first inaccuracy of previous geometry, where all rows where made of 32 holes.
  //
  // The A-B-A-B scheme is maintained until the very last two row on the top of the module, where the last row
  // is a A instead then a B.
  // This would be a minor difference with reality and it is probably not worth implementing
  //
  // The holes are staggered between consecutive rows, both in real module and in the previous geometry
  //
  // The module is made of 2 sectors (FRONT and BACK), each made of 9 cells, in a 3x3 array. The real
  // cell dimension is no the same for all cells. This is a consequence of the type A-B scheme. In fact,
  // both rows are crossing 3 cells, so they are divided in 3 sectors, each sector defined by a given
  // type of crystal. In particular:
  //
  // 1) Type A rows are made of (11 + 11 + 11) crystals
  // 2) Type B rows are made of (10 + 11 + 11) crystals
  //
  // In previous geometry, the 10 crystals sector was on the right (in this viewpoint defined above)
  // and not on the left.
  //
  // new placement is based on Cell class (see the Cell.hh and Cell.cc files)


  // bool wireFrame = true;


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
      smodule << "module_" << iMod << "_" << jMod;
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
      G4double gapSize = 0.01 *mm;
      // first prepare the logical volumes of the light guides, then they will be placed in space with a loop
      // Light guide structure
      // Truncated of cone, cut on the sides
      // goes from 20x20 mm2 (on the fibers) to 5 mm radius on the PMT
      // first, do a simple truncated cone
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
      Interface_Positive_Z_LV = new G4LogicalVolume (Interface_Positive_Z_S, WoMaterial, sname.str().c_str()) ;
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

      // now the light guides



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
      Interface_Negative_Z_LV = new G4LogicalVolume (Interface_Negative_Z_S, WoMaterial, sname.str().c_str()) ;
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
      //-----------------------------------//


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
      // The readout
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
      Readout_Positive_Z_LV = new G4LogicalVolume (Readout_Positive_Z_S, WoMaterial, sname.str().c_str()) ;
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
      Readout_Negative_Z_LV = new G4LogicalVolume (Readout_Negative_Z_S, WoMaterial, sname.str().c_str()) ;
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
  if      ( gap_abs_interface == 0 ) GapAbsToInterfaceMaterial = MyMaterials::Air () ;
  else if ( gap_abs_interface == 1 ) GapAbsToInterfaceMaterial = MyMaterials::OpticalGrease () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid GapAbsToInterfaceMaterial material specifier " << gap_abs_interface << G4endl ;
    exit (-1) ;
  }
  G4cout << "Cone material: "<< GapAbsToInterfaceMaterial << G4endl ;

  GapInterfaceToReadoutMaterial  = NULL;
  if      ( gap_interface_readout == 0 ) GapInterfaceToReadoutMaterial = MyMaterials::Air () ;
  else if ( gap_interface_readout == 1 ) GapInterfaceToReadoutMaterial = MyMaterials::OpticalGrease () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid GapInterfaceToReadoutMaterial material specifier " << gap_interface_readout << G4endl ;
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
  // PLEXMaterial = NULL ;
  // PLEXMaterial = MyMaterials::PLEX () ;
  //
  // G4cout << "PLEX. material: "<< PLEXMaterial << G4endl ;
  //
  // PVCMaterial = NULL ;
  // PVCMaterial = MyMaterials::PVC () ;
  //
  // G4cout << "PVC. material: "<< PVCMaterial << G4endl ;
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

void ConstructInterface(G4LogicalVolume *moduleLV,  std::string moduleName)
{
  
}
