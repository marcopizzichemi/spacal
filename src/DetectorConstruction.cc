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

using namespace CLHEP;



DetectorConstruction::DetectorConstruction (const string& configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------

  ConfigFile config (configFileName) ;

  config.readInto (checkOverlaps, "checkOverlaps") ;

  config.readInto (world_material, "world_material") ;

  config.readInto (abs_material, "abs_material") ;
  config.readInto (W_fraction, "W_fraction") ;
  config.readInto (hole_radius, "hole_radius") ;
  config.readInto (module_xy, "module_xy") ;
  config.readInto (module_yx, "module_yx") ;
  config.readInto (fibres_x, "fibres_x");
  config.readInto (fibres_y, "fibres_y");
  config.readInto (fibres_x1, "fibres_x1");
  config.readInto (fibres_y1, "fibres_y1");
  //config.readInto (absorber_x, "absorber_x");
  //  config.readInto (absorber_y, "absorber_y");
  config.readInto (postshower, "postshower") ;
  config.readInto (second, "second") ;
  //  config.readInto (secondAbs, "secondAbs") ;
  config.readInto (Second_abs_material, "Second_abs_material") ;
  config.readInto (Second_module_xy, "Second_module_xy");
  config.readInto (Second_module_yx, "Second_module_yx") ;
  config.readInto (Second_W_fraction, "Second_W_fraction") ;
  config.readInto (Second_fibres_x, "Second_fibres_x");
  config.readInto (Second_fibres_y, "Second_fibres_y");
  config.readInto (Second_fibres_x1, "Second_fibres_x1");
  config.readInto (Second_fibres_y1, "Second_fibres_y1");

  config.readInto (fibre_scheme, "fibre_scheme") ;
  config.readInto (Second_fibre_scheme, "Second_fibre_scheme") ;
  config.readInto (fibre_material, "fibre_material") ;
  config.readInto (fibre_material1, "fibre_material1") ;
  config.readInto (fibre_material2, "fibre_material2") ;
  config.readInto (Second_fibre_material, "Second_fibre_material") ;
  config.readInto (Second_fibre_material1, "Second_fibre_material1") ;
  config.readInto (Second_fibre_material2, "Second_fibre_material2") ;

  config.readInto (fibre_cladRIndex, "fibre_cladRIndex") ;
  config.readInto (fibre_isSquare, "fibre_isSquare") ;
  config.readInto (fibre_radius, "fibre_radius") ;
  config.readInto (fibre_length, "fibre_length") ;
  config.readInto (Second_fibre_radius, "Second_fibre_radius") ;
  config.readInto (Second_fibre_length, "Second_fibre_length") ;
  config.readInto (fibre_distance, "fibre_distance") ;
  config.readInto (Second_fibre_distance, "Second_fibre_distance") ;
  config.readInto (fibre_absLength, "fibre_absLength") ;

  config.readInto (gap_l, "gap_l") ;
  config.readInto (gap_material, "gap_material") ;

  config.readInto (det_l, "det_l") ;
  config.readInto (det_material, "det_material") ;

  config.readInto (depth, "depth") ;

  config.readInto (lead_plane, "lead_plane") ;
  config.readInto (lp_dist, "lp_dist") ;
  config.readInto (lp_depth, "lp_depth") ;
  config.readInto (lp_mat, "lp_mat") ;
  config.readInto (lp_x, "lp_x") ;
  config.readInto (lp_y, "lp_y") ;
  config.readInto (preconstr, "preconstr") ;


  config.readInto (surface_lg,"surface_lg");
  config.readInto (glue_interface,"glue_interface");
  config.readInto (cone_material,"cone_material");

  B_field_intensity = config.read<double>("B_field_intensity") * tesla ;


  margin = std::max( 0.25*fibre_distance, 2.*fibre_radius );
  G4double staggering = 0.5*fibre_distance*((fibre_scheme+1)%2);
  G4double staggeredStart = 0.5 + 0.5*(fibre_scheme%2);

  std::cout << "staggeredStart: " << staggeredStart << std::endl;
  if( fibre_scheme == 1 || fibre_scheme == 2 ) // dice-4
  {
    fibreDistanceAlongX = fibre_distance;
    fibreDistanceAlongY = fibre_distance;
    nFibresAlongX = floor( (fibres_x - 2.*margin) / fibreDistanceAlongX ) + 1 ;
    nFibresAlongY = floor( (fibres_y - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
    nFibresAlongX1 = floor( (fibres_x1 - 2.*margin) / fibreDistanceAlongX ) + 1 ;
    nFibresAlongY1 = floor( (fibres_y1 - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
    startX = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX - 1.0) ) ;
    startY = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY - staggeredStart)) ;
    startX1 = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX1 - 1.0) ) ;
    startY1 = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY1 - staggeredStart)) ;
  }
  else if( fibre_scheme == 3 || fibre_scheme == 4 ) // dice-5
  {
    fibreDistanceAlongX = 0.8660 * fibre_distance;
    fibreDistanceAlongY = fibre_distance;
    nFibresAlongX = floor( (fibres_x - 2.*margin)/ fibreDistanceAlongX ) + 1 ;
    nFibresAlongY = floor( (fibres_y - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
    nFibresAlongX1 = floor( (fibres_x1 - 2.*margin) / fibreDistanceAlongX ) + 1 ;
    nFibresAlongY1 = floor( (fibres_y1 - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
    startX = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX - 1.0) ) ;
    startY = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY - staggeredStart)) ;
    startX1 = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX1 - 1.0) ) ;
    startY1 = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY1 - staggeredStart)) ;
  }
  else if( fibre_scheme == 5 || fibre_scheme == 6 ) // chessboard
  {
    fibreDistanceAlongX = 0.5 * fibre_distance;
    fibreDistanceAlongY = fibre_distance;
    nFibresAlongX = floor( (fibres_x - 2.*margin)              / fibreDistanceAlongX ) + 1 ;
    nFibresAlongY = floor( (fibres_y - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
    nFibresAlongX1 = floor( (fibres_x1 - 2.*margin)              / fibreDistanceAlongX ) + 1 ;
    nFibresAlongY1 = floor( (fibres_y1 - 2.*margin - staggering) / fibreDistanceAlongY ) + 1 ;
    startX = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX - 1.0) ) ;
    startY = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY - staggeredStart)) ;
    startX1 = 0.5 * ( module_xy - fibreDistanceAlongX * (nFibresAlongX1 - 1.0) ) ;
    startY1 = 0.5 * ( module_yx - fibreDistanceAlongY * (nFibresAlongY1 - staggeredStart)) ;
  }

  //________________________________________________________________________________


  margin2 = std::max( 0.25*Second_fibre_distance, 2.*Second_fibre_radius );
  G4double staggering2 = 0.5*Second_fibre_distance*((Second_fibre_scheme+1)%2);
  G4double staggeredStart2 = 0.5 + 0.5*(Second_fibre_scheme%2);

  std::cout << "staggeredStart: " << staggeredStart2 << std::endl;
  if( Second_fibre_scheme == 1 || Second_fibre_scheme == 2 ) // dice-4
  {
    Second_fibreDistanceAlongX = Second_fibre_distance;
    Second_fibreDistanceAlongY = Second_fibre_distance;
    Second_nFibresAlongX = floor( (Second_fibres_x - 2.*margin2) / Second_fibreDistanceAlongX ) + 1 ;
    Second_nFibresAlongY = floor( (Second_fibres_y - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
    Second_nFibresAlongX1 = floor( (Second_fibres_x1 - 2.*margin2) / Second_fibreDistanceAlongX ) + 1 ;
    Second_nFibresAlongY1 = floor( (Second_fibres_y1 - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
    Second_startX = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX - 1.0) ) ;
    Second_startY = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY - staggeredStart2)) ;
    Second_startX1 = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX1 - 1.0) ) ;
    Second_startY1 = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY1 - staggeredStart)) ;
  }
  else if( Second_fibre_scheme == 3 || Second_fibre_scheme == 4 ) // dice-5
  {
    Second_fibreDistanceAlongX = 0.8660 * Second_fibre_distance;
    Second_fibreDistanceAlongY = Second_fibre_distance;
    Second_nFibresAlongX = floor( (Second_fibres_x - 2.*margin2)/ Second_fibreDistanceAlongX ) + 1 ;
    Second_nFibresAlongY = floor( (Second_fibres_y - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
    Second_nFibresAlongX1 = floor( (Second_fibres_x1 - 2.*margin2) / Second_fibreDistanceAlongX ) + 1 ;
    Second_nFibresAlongY1 = floor( (Second_fibres_y1 - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
    Second_startX = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX - 1.0) ) ;
    Second_startY = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY - staggeredStart2)) ;
    Second_startX1 = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX1 - 1.0) ) ;
    Second_startY1 = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY1 - staggeredStart)) ;

  }
  else if( Second_fibre_scheme == 5 || Second_fibre_scheme == 6 ) // chessboard
  {
    Second_fibreDistanceAlongX = 0.5 * Second_fibre_distance;
    Second_fibreDistanceAlongY = Second_fibre_distance;
    Second_nFibresAlongX = floor( (Second_fibres_x - 2.*margin2)              / Second_fibreDistanceAlongX ) + 1 ;
    Second_nFibresAlongY = floor( (Second_fibres_y - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
    Second_nFibresAlongX1 = floor( (Second_fibres_x1 - 2.*margin2)              / Second_fibreDistanceAlongX ) + 1 ;
    Second_nFibresAlongY1 = floor( (Second_fibres_y1 - 2.*margin2 - staggering2) / Second_fibreDistanceAlongY ) + 1 ;
    Second_startX = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX - 1.0) ) ;
    Second_startY = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY - staggeredStart2)) ;
    Second_startX1 = 0.5 * ( Second_module_xy - Second_fibreDistanceAlongX * (Second_nFibresAlongX1 - 1.0) ) ;
    Second_startY1 = 0.5 * ( Second_module_yx - Second_fibreDistanceAlongY * (Second_nFibresAlongY1 - staggeredStart)) ;
  }


  /////////_____________________________________________________________________________




  module_z = fibre_length;
  PLEX_depth = 30;
  PLEX_dist = 2*fibre_length;
  PVC_depth = 20;
  PVC_dist = 200 + 2*fibre_length;
  PMT_radius = 5;
  PMT_length = 60;
  Wires_radius = 1.8;
  Wires_dist = 50;


  Second_module_z = Second_fibre_length;

  expHall_x = module_xy * 4 ;
  expHall_y = module_yx * 4 ;
  expHall_z = (module_z + Second_module_z ) * 4 + lp_dist + lp_depth ;
  if( postshower ) expHall_z *= 3;

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
    G4VSolid * worldS = new G4Box ("worldS", 0.5 * expHall_x, 0.5 * expHall_y, 0.5 * expHall_z) ;
    G4LogicalVolume * worldLV = new G4LogicalVolume (worldS, WoMaterial, "worldLV", 0, 0, 0) ;
    G4VPhysicalVolume * worldPV = new G4PVPlacement (0, G4ThreeVector (), worldLV, "World", 0, false, 0, checkOverlaps) ;
    //******************************************___ I'M here now ****************

    // The calorimeter
    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    G4VSolid * calorimeterS = new G4Box ("calorimeterS", 0.5 * module_xy, 0.5 * module_yx, 0.5 * (module_z + Second_module_z)) ;
    G4LogicalVolume * calorimeterLV = new G4LogicalVolume (calorimeterS, WoMaterial, "calorimeterLV") ;
    if( !postshower )
    {
      new G4PVPlacement (0, G4ThreeVector (0., 0., 0.), calorimeterLV, "calorimeterPV", worldLV, false, 0, checkOverlaps) ;
    }
    else
    {
      new G4PVPlacement (0, G4ThreeVector (0., 0., 0.), calorimeterLV, "calorimeterPV", worldLV, false, 0, checkOverlaps) ;
      new G4PVPlacement (0, G4ThreeVector (0., 0., 1.*(module_z + Second_module_z)), calorimeterLV, "postshower1PV", worldLV, false, 0, checkOverlaps) ;
      new G4PVPlacement (0, G4ThreeVector (0., 0., 2.*(module_z + Second_module_z)), calorimeterLV, "postshower2PV", worldLV, false, 0, checkOverlaps) ;
    }

    // Plane in front of Calo
    if ( lead_plane == 1)
    {
      G4VSolid * PLANE = new G4Box ("PLANE", 0.5*lp_x, 0.5*lp_y, 0.5*lp_depth);
      G4LogicalVolume * PLANELV = new G4LogicalVolume (PLANE, PlaneMaterial, "PLANELV");
      new G4PVPlacement (0, G4ThreeVector (0.,0.,-0.5*(lp_dist + lp_depth)),PLANELV,"PLANEPV", worldLV, false, 0 ,checkOverlaps);
    }





    // G4VSolid * PLEXS = new G4Box ("PLEXS", 0.5*module_xy, 0.5*module_yx, 0.5*PLEX_depth);


    // FIXME hardcoded for now...
    // airgap
    G4double airGap = 0.1*mm; // air gap light guide to pmts and to fibers
    G4double airGap_abs = 0.1*mm; //air gap between absorbers
    // Trapezoid shape for light guides
    G4double lguide_dx1 = 10*mm;
    G4double lguide_dy1 = 10*mm;
    G4double lguide_dx2 = 20*mm;
    G4double lguide_dy2 = 20*mm;
    G4double lguide_dz  = PLEX_depth;
    G4double pmts_pitch = 20*mm;

    // G4Trd* LguideS = new G4Trd("LguideS",  // name
                                   // 0.5*lguide_dx1, 0.5*lguide_dx2,
                                   // 0.5*lguide_dy1, 0.5*lguide_dy2, 0.5*lguide_dz); //its size
    // G4LogicalVolume * LguideLV = new G4LogicalVolume (LguideS,PLEXMaterial,"LguideLV");


    // create light guide structure
    G4VSolid* coneS = new G4Cons("aCone",0,PMT_radius,0, 0.5*sqrt(2)*20.0*mm ,0.5*PLEX_depth,0.0 * deg,360.0* deg);
    G4Box* innerBoxS = new G4Box ("innerBoxS", 0.5*20.0*mm, 0.5*20.0*mm, 0.5*PLEX_depth) ;
    G4Box* outerBoxS = new G4Box ("outerBoxS", 0.5*sqrt(2)*20.0*mm, 0.5*sqrt(2)*20.0*mm, 0.5*PLEX_depth) ;
    G4VSolid* subtract = new G4SubtractionSolid("Hollow-Box", outerBoxS, innerBoxS,
    0,  G4ThreeVector(0.,0.,0.));
    G4VSolid* coneSolid
    = new G4SubtractionSolid("coneSolid", coneS, subtract,
    0,  G4ThreeVector(0.,0.,0.));

    // G4LogicalVolume * LguideLV = new G4LogicalVolume (coneSolid,WoMaterial,"LguideLV");
    G4LogicalVolume * LguideLV = new G4LogicalVolume (coneSolid,ConeMaterial,"LguideLV"); //plexiglass
    //air light guide
    G4VSolid * air_LG_front_S = new G4Box ("air_LG_front_S ",0.5*module_xy, 0.5*module_yx,0.5*PLEX_depth);
    G4VSolid * air_LG_back_S = new G4Box ("air_LG_back_S ",0.5*module_xy, 0.5*module_yx,0.5*PLEX_depth);
    G4LogicalVolume *air_LG_front_LV = new G4LogicalVolume (air_LG_front_S,WoMaterial, "air_LG_front_LV");
    G4LogicalVolume *air_LG_back_LV = new G4LogicalVolume (air_LG_back_S,WoMaterial, "air_LG_back_LV");




    // G4LogicalVolume * coneLV = new G4LogicalVolume (coneSolid,PLEXMaterial,"coneLV");
    // G4PVPlacement *conePV = new G4PVPlacement (0,G4ThreeVector(0,0,-30.0 *cm), coneLV, "conePV", worldLV, false, 0, checkOverlaps) ;


    G4VSolid * PVCS = new G4Box ("PVCS", 0.5*module_xy, 0.5*module_yx, 0.5*PVC_depth);
    G4LogicalVolume * PVCLV = new G4LogicalVolume (PVCS,PVCMaterial, "PVCLV");
    G4VSolid * TubeS = new G4Tubs ("TubeS", 0., PMT_radius, 0.5*PMT_length, 0.*deg, 360.*deg);

    G4VSolid * WiresS = new G4Tubs ("WiresS", 0., Wires_radius, Wires_dist, 0.*deg, 360.*deg);

    G4LogicalVolume * PMTLV = new G4LogicalVolume (TubeS,PLEXMaterial, "PMTLV");
    G4LogicalVolume * WiresLV = new G4LogicalVolume (WiresS, WiresMaterial, "WiresLV");

    G4VSolid * PVC_pmt_frontS = new G4Box ("PVC_pmt_frontS",0.5*module_xy, 0.5*module_yx,0.5*PMT_length);
    G4VSolid * PVC_pmt_backS  = new G4Box ("PVC_pmt_backS" ,0.5*module_xy, 0.5*module_yx,0.5*PMT_length);
    G4LogicalVolume * PVC_pmt_frontLV = new G4LogicalVolume (PVC_pmt_frontS,PVCMaterial, "PVC_pmt_frontLV");
    G4LogicalVolume * PVC_pmt_backLV  = new G4LogicalVolume (PVC_pmt_backS,PVCMaterial, "PVC_pmt_frontLV");



    G4PVPlacement *PVC_pmt_frontPV = new G4PVPlacement (0, G4ThreeVector (0, 0,-0.5*(PLEX_dist+PMT_length+PLEX_depth*2+airGap*4 + airGap_abs) ),PVC_pmt_frontLV, "PVC_pmt_frontPV", worldLV, false, 0, checkOverlaps);
    G4PVPlacement *PVC_pmt_backPV = new G4PVPlacement (0, G4ThreeVector (0, 0,+0.5*(PLEX_dist+PMT_length+PLEX_depth*2+airGap*4 + airGap_abs) ),PVC_pmt_backLV, "PVC_pmt_backPV", worldLV, false, 0, checkOverlaps);

    G4PVPlacement *air_LG_front_PV = new G4PVPlacement(0, G4ThreeVector(0 ,0 , -0.5*(PLEX_dist + PLEX_depth + airGap*2 + airGap_abs)), air_LG_front_LV, "air_LG_front_PV", worldLV, false, 0, checkOverlaps);
    G4PVPlacement  *air_LG_back_PV = new G4PVPlacement(0, G4ThreeVector(0 ,0 , +0.5*(PLEX_dist + PLEX_depth + airGap*2 + airGap_abs)), air_LG_back_LV , "air_LG_back_PV" , worldLV, false, 0, checkOverlaps);

    G4VSolid *glue_S = new G4Box ("glue_S", 0.5*module_xy, 0.5*module_yx, 0.5*airGap);
    G4LogicalVolume * glue_LV = new G4LogicalVolume (glue_S,GlueMaterial, "glue_LV");

    G4PVPlacement  *gluePV_front = new G4PVPlacement(0, G4ThreeVector(0 ,0 , -0.5*(PLEX_dist + airGap + airGap_abs )), glue_LV, "gluePV_front", worldLV, false, 0, checkOverlaps);
    G4PVPlacement  *gluePV_back  = new G4PVPlacement(0, G4ThreeVector(0 ,0 , +0.5*(PLEX_dist + airGap + airGap_abs )), glue_LV, "gluePV_back" , worldLV, false, 0, checkOverlaps);

    G4PVPlacement  *glue_pmt_PV_front = new G4PVPlacement(0, G4ThreeVector(0 ,0 , -0.5*(PLEX_dist + 3*airGap + airGap_abs +PLEX_depth*2)), glue_LV, "glue_pmt_PV_front", worldLV, false, 0, checkOverlaps);
    G4PVPlacement  *glue_pmt_PV_back  = new G4PVPlacement(0, G4ThreeVector(0 ,0 , +0.5*(PLEX_dist + 3*airGap + airGap_abs +PLEX_depth*2)), glue_LV, "glue_pmt_PV_back" , worldLV, false, 0, checkOverlaps);





    if (preconstr == 1)
    {
      // The Pre-Construction
      // mods for real spacal test beam structure
      //----------- ---------------





      // new G4PVPlacement (0, G4ThreeVector (0.,0.,-0.5*(PVC_dist+PVC_depth)),PVCLV,"PVCPV", worldLV, false, 0, checkOverlaps);

      // new G4PVPlacement (0, G4ThreeVector (0.,0.,-(PLEX_dist+PVC_depth+60)),PVCLV,"PVCPV", worldLV, false, 0, checkOverlaps);

      startAX = 0. + 0.5*(absorber_x);
      nCellsAlongX =  floor( module_xy / absorber_x );
      nCellsAlongY  = floor ( module_yx / absorber_y);
      startAY = 0.+ 0.5*(absorber_y + (module_yx-module_xy));

      // The PMTS
      // ------------------------------------- --------  ------


      int NX = 1;


      float pmt_x[3] = {-pmts_pitch,0,pmts_pitch};
      float pmt_y[3] = {-pmts_pitch,0,pmts_pitch};
      for(int iPMT = 0 ; iPMT < 3; iPMT++)
      {
        for(int jPMT = 0 ; jPMT < 3; jPMT++)
        {
          float PMTx = pmt_x[iPMT];
          float PMTy = pmt_y[jPMT];

          int iP = 3* iPMT + jPMT;
          std::string Pname;

          //-----------//
          // FRONT     //
          //-----------//

          Pname = Form("PMTSPV FRONT %d", iP);
          new G4PVPlacement (0, G4ThreeVector (PMTx, PMTy,0 ),PMTLV, Pname, PVC_pmt_frontLV, false, 0, checkOverlaps);

          Pname = Form("LguidePV FRONT %d", iP);
          G4PVPlacement *LG_front_PV = new G4PVPlacement (0, G4ThreeVector (PMTx,PMTy,0), LguideLV, Pname, air_LG_front_LV, false, 0, checkOverlaps);

          //optical surface between cone and air cone
          if(surface_lg)
          {
            std::stringstream Surfname_front;
            Surfname_front << "Surface_" << Pname << "_air_front";
            G4OpticalSurface* reflector_surf_front = new G4OpticalSurface(Surfname_front.str().c_str());
            reflector_surf_front->SetType(dielectric_metal);
            reflector_surf_front->SetFinish(polished);
            reflector_surf_front->SetModel(unified);
            reflector_surf_front->SetMaterialPropertiesTable(MyMaterials::ESR());
            G4LogicalBorderSurface* reflectorLB_front;
            Surfname_front << "LB";
            reflectorLB_front = new G4LogicalBorderSurface(Surfname_front.str().c_str(),
                                                          LG_front_PV,
                                                          air_LG_front_PV,
                                                         reflector_surf_front);
            //optical surface between air cone and cone
            std::stringstream SurfnameInv_front;
            SurfnameInv_front << "SurfaceInv_" << Pname << "_air_front";
            G4OpticalSurface* reflectorInv_surf_front = new G4OpticalSurface(SurfnameInv_front.str().c_str());
            reflectorInv_surf_front->SetType(dielectric_metal);
            reflectorInv_surf_front->SetFinish(polished);
            reflectorInv_surf_front->SetModel(unified);
            reflectorInv_surf_front->SetMaterialPropertiesTable(MyMaterials::ESR());
            G4LogicalBorderSurface* reflectorLBInv_front;
            SurfnameInv_front << "LB";
            reflectorLBInv_front = new G4LogicalBorderSurface(SurfnameInv_front.str().c_str(),
                                                              air_LG_front_PV,
                                                              LG_front_PV,
                                                             reflectorInv_surf_front);
          }

          //-----------//
          // BACK      //
          //-----------//

          G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
          rotationMatrix->rotateY(180.*deg);

          Pname = Form("PMTSPV BACK %d", iP);
          new G4PVPlacement (0, G4ThreeVector (PMTx, PMTy,0 ),PMTLV, Pname, PVC_pmt_backLV, false, 0, checkOverlaps);

          Pname = Form("LguidePV BACK %d", iP);
          G4PVPlacement *LG_back_PV = new G4PVPlacement (rotationMatrix, G4ThreeVector (PMTx,PMTy, 0) , LguideLV, Pname, air_LG_back_LV, false, 0, checkOverlaps);


          if(surface_lg)
          {
            //optical surface between cone and air cone
            std::stringstream Surfname_back;
            Surfname_back << "Surface_" << Pname << "_air_back";
            G4OpticalSurface* reflector_surf_back = new G4OpticalSurface(Surfname_back.str().c_str());
            reflector_surf_back->SetType(dielectric_metal);
            reflector_surf_back->SetFinish(polished);
            reflector_surf_back->SetModel(unified);
            reflector_surf_back->SetMaterialPropertiesTable(MyMaterials::ESR());
            G4LogicalBorderSurface* reflectorLB_back;
            Surfname_back << "LB";
            reflectorLB_back = new G4LogicalBorderSurface(Surfname_back.str().c_str(),
                                                          LG_back_PV,
                                                          air_LG_back_PV,
                                                         reflector_surf_back);
            //optical surface between air cone and cone
            std::stringstream SurfnameInv_back;
            SurfnameInv_back << "SurfaceInv_" << Pname << "_air_back";
            G4OpticalSurface* reflectorInv_surf_back = new G4OpticalSurface(SurfnameInv_back.str().c_str());
            reflectorInv_surf_back->SetType(dielectric_metal);
            reflectorInv_surf_back->SetFinish(polished);
            reflectorInv_surf_back->SetModel(unified);
            reflectorInv_surf_back->SetMaterialPropertiesTable(MyMaterials::ESR());
            G4LogicalBorderSurface* reflectorLBInv_back;
            SurfnameInv_back << "LB";
            reflectorLBInv_back = new G4LogicalBorderSurface(SurfnameInv_back.str().c_str(),
                                                              air_LG_back_PV,
                                                              LG_back_PV,
                                                             reflectorInv_surf_back);

          }


        }
      }

      // OLD PMT PLACEMENT CODE
      // for (float xP = -0.5 * module_xy + module_xy/6; NX < 4; xP += module_xy/3)
      // {
      //   int NY = 1;
      //   for (float yP = -0.5 * module_yx + module_yx/6; NY < 4; yP += module_yx/3)
      //   {
      //     float PMTx = xP;
      //     float PMTy = yP;
      //
      //     int iP = 3* NX + NY;
      //
      //     std::string Pname;
      //     Pname = Form("PMTSPV %d", iP);
      //
          // std::string Wname;
          // Wname = Form("WiresSPV %d", iP);
      //
      //     new G4PVPlacement (0, G4ThreeVector (PMTx, PMTy,-0.5*(PLEX_dist+PMT_length+PLEX_depth*2)),PMTLV, Pname, worldLV, false, 0, checkOverlaps);
          // new G4PVPlacement (0, G4ThreeVector (PMTx, PMTy,-0.5*(PVC_dist + PVC_depth*2 )),WiresLV, Wname, worldLV, false, 0, checkOverlaps);
      //
      //     ++NY;
      //   }
      //   ++NX;
      // }
    }

    // The holes
    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    G4VSolid * holeS;
    if( !fibre_isSquare ) holeS = new G4Tubs ("holeS", fibre_radius, fibre_radius+hole_radius, 0.5*fibre_length, 0.*deg, 360.*deg) ;
    else
    {
      G4VSolid * temp1 = new G4Box ("temp1", fibre_radius+hole_radius, fibre_radius+hole_radius, 0.5*fibre_length) ;
      G4VSolid * temp2 = new G4Box ("temp2", fibre_radius, fibre_radius, 1.6*fibre_length) ;
      holeS = new G4SubtractionSolid("holeS",temp1,temp2,0,G4ThreeVector(0.,0.,0.));
    }
    G4LogicalVolume * holeLV = new G4LogicalVolume (holeS, WoMaterial, "holeLV") ;
    //HoleParameterisation* holeParam = new HoleParameterisation(module_xy,fibre_radius,hole_radius,fibre_distance,fibre_length,WoMaterial);
    //new G4PVParameterised("holeP", holeLV, absorberLV, kUndefined, holeParam->GetNHoles(), holeParam);


    // the fibres
    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    G4VSolid * fibre1S;
    if( !fibre_isSquare ) fibre1S = new G4Tubs ("fibre1S", 0., fibre_radius, 0.5*fibre_length, 0.*deg, 360.*deg) ;
    else                  fibre1S = new G4Box ("fibre1S", fibre_radius, fibre_radius, 0.5*fibre_length) ;
    G4LogicalVolume * fibre1SLV = new G4LogicalVolume (fibre1S, ClMaterial, "fibre1LV") ;
    //FibreParameterisation* fibreParam = new FibreParameterisation(module_xy,fibre_radius,fibre_distance,fibre_length,ClMaterial);
    //new G4PVParameterised("fibreP", fibreLV, absorberLV, kUndefined, fibreParam->GetNFibres(), fibreParam);


    //___________________________________________________________


    G4VSolid * fibre12S;
    if( !fibre_isSquare ) fibre12S = new G4Tubs ("fibre12S", 0., fibre_radius, 0.5*fibre_length, 0.*deg, 360.*deg) ;
    else                  fibre12S = new G4Box ("fibre12S", fibre_radius, fibre_radius, 0.5*fibre_length) ;
    G4LogicalVolume * fibre12SLV = new G4LogicalVolume (fibre12S, ClSSMaterial, "fibre12SLV") ;

    G4VSolid * fibre13S;
    if( !fibre_isSquare ) fibre13S = new G4Tubs ("fibre13S", 0., fibre_radius, 0.5*fibre_length, 0.*deg, 360.*deg) ;
    else                  fibre13S = new G4Box ("fibre13S", fibre_radius, fibre_radius, 0.5*fibre_length) ;
    G4LogicalVolume * fibre13SLV = new G4LogicalVolume (fibre13S, ClSSSMaterial , "fibre13SLV") ;


    G4VSolid * fibre2S;
    if( !fibre_isSquare ) fibre2S = new G4Tubs ("fibre2S", 0., fibre_radius, 0.5*Second_fibre_length, 0.*deg, 360.*deg) ;
    else                  fibre2S = new G4Box ("fibre2S", Second_fibre_radius, Second_fibre_radius, 0.5*Second_fibre_length) ;
    G4LogicalVolume * fibre2SLV = new G4LogicalVolume (fibre2S, Cl3Material, "fibre2LV") ;
    //FibreParameterisation* fibreParam = new FibreParameterisation(module_xy,fibre_radius,fibre_distance,fibre_length,ClMaterial);
    //new G4PVParameterised("fibreP", fibreLV, absorberLV, kUndefined, fibreParam->GetNFibres(), fibreParam);


    //___________________________________________________________


    G4VSolid * fibre22S;
    if( !fibre_isSquare ) fibre22S = new G4Tubs ("fibre22S", 0., Second_fibre_radius, 0.5*Second_fibre_length, 0.*deg, 360.*deg) ;
    else                  fibre22S = new G4Box ("fibre22S", Second_fibre_radius, Second_fibre_radius, 0.5*Second_fibre_length) ;
    G4LogicalVolume * fibre22SLV = new G4LogicalVolume (fibre22S, Cl4SSMaterial, "fibre2SSLV") ;


    G4VSolid * fibre23S;
    if( !fibre_isSquare ) fibre23S = new G4Tubs ("fibre23S", 0., Second_fibre_radius, 0.5*Second_fibre_length, 0.*deg, 360.*deg) ;
    else                  fibre23S = new G4Box ("fibre23S", Second_fibre_radius, Second_fibre_radius, 0.5*Second_fibre_length) ;
    G4LogicalVolume * fibre23SLV = new G4LogicalVolume (fibre23S, Cl43SSMaterial , "fibre23SLV") ;

    //_________________________________________________________________





    G4cout << "Sections: " << nCellsAlongX << "   " << nCellsAlongY << G4endl;




    // The Absorber

    // G4double airGap_abs = 0.1*mm;

    G4VSolid * absorber1S = new G4Box ("absorber1S", 0.5 * module_xy, 0.5 * (module_yx), 0.5 * module_z) ;
    G4LogicalVolume * absorber1LV;
    absorber1LV = new G4LogicalVolume (absorber1S, AbMaterial, "absorber1LV") ;
    G4PVPlacement *absorber1PV =  new G4PVPlacement (0, G4ThreeVector (0., 0., - 0.5* (module_z + airGap_abs)), absorber1LV, "absorber1PV", calorimeterLV, false, 0, checkOverlaps) ;



    // Second Absorber
    G4VSolid * absorber2S = new G4Box ("absorber2S", 0.5 * Second_module_xy, 0.5 * Second_module_yx, 0.5 * Second_module_z) ;
    G4LogicalVolume * absorber2LV;
    absorber2LV = new G4LogicalVolume (absorber2S, AbMaterial2, "absorber2LV") ;
    G4PVPlacement *absorber2PV = new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5* (module_z + airGap_abs)), absorber2LV, "absorber2PV", calorimeterLV, false, 0, checkOverlaps) ;

    //thin layer of air between absorbers

    G4VSolid * airLayer = new G4Box ("airLayer", 0.5 * Second_module_xy, 0.5 * Second_module_yx, 0.5 * airGap_abs ) ;
    G4LogicalVolume * airLayerLV;
    airLayerLV = new G4LogicalVolume (airLayer, WoMaterial, "airLayerLV") ;
    G4PVPlacement* airLayerPV = new G4PVPlacement (0, G4ThreeVector (0., 0., 0.), airLayerLV, "airLayerPV", calorimeterLV, false, 0, checkOverlaps) ;


    // fibres matrix filling
    // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
    int index;
    int indexion;
    int index1;
    // loop on x direction
    int countX = 0 ; // for the staggering

    for (float x = -0.5*module_xy+startX; countX <nFibresAlongX; x += fibreDistanceAlongX)
    {
      // loop on y direction
      int countY = 0 ; // for the staggering
      for (float y = - 0.5 * module_yx + startY; countY < nFibresAlongY; y += fibreDistanceAlongY)
      {
        float x_c = x;
        float y_c = y;
        //____________________________________________
        //   float x_cs = -0.5 * module_xy;
        //    float y_cs = -0.5 * module_xy;
        //________________________________________________________________
        // staggering
        if( (fibre_scheme%2) == 0 )
        y_c += 0.5*fibreDistanceAlongY*(countX%2) ;

        int index = countX * nFibresAlongY + countY ;

        CreateTree::Instance() -> fibresPosition -> Fill(index,x_c,y_c);


        std::string name;


        // **************  ***************
        // #5
        G4PVPlacement *fiberPV;
        G4PVPlacement *holePV;
        if (x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1)
        {
          name = Form("holePV %d",index);
          holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber1LV, false, 0, checkOverlaps) ;

          name = Form("fibre12SPV %d",index);
          fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre12SLV, name, absorber1LV, false, 0, checkOverlaps) ;
        }
        // #2,#4,#6,#8
        else if (x >= -1.5*fibres_x1 && x<= fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1 || x >= fibres_x1 && x<= 1.5*fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1 || x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= -1.5*fibres_y1 && y <= -0.5*fibres_y1 ||  x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= 0.5*fibres_y1 && y <= 1.5*fibres_y1 )
        {
          name = Form("holePV %d",index);
          holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber1LV, false, 0, checkOverlaps) ;

          name = Form("fibre1SPV %d",index);
          fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre1SLV, name, absorber1LV, false, 0, checkOverlaps) ;
        }

        else
        {
          name = Form("holePV %d",index);
          holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber1LV, false, 0, checkOverlaps) ;

          name = Form("fibre13SPV %d",index);
          fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre13SLV, name, absorber1LV, false, 0, checkOverlaps) ;
      }

      indexion = index;
      index1 = index+1;
      CreateTree::Instance() -> fibresPosition_1st_Section -> Fill(indexion,x_c,y_c);
      ++countY ;



      //optical surface between fiber and thin air layer
      std::stringstream Surfname;
      Surfname << "Surface_" << name << "_airLayer";
      G4OpticalSurface* reflector_surf = new G4OpticalSurface(Surfname.str().c_str());
      reflector_surf->SetType(dielectric_metal);
      reflector_surf->SetFinish(polished);
      reflector_surf->SetModel(unified);
      reflector_surf->SetMaterialPropertiesTable(MyMaterials::ESR());
      G4LogicalBorderSurface* reflectorLB;
      Surfname << "LB";
      reflectorLB = new G4LogicalBorderSurface(Surfname.str().c_str(),
                                               fiberPV,
                                               airLayerPV,
                                               reflector_surf);
      //optical surf from air gap to abs
      std::stringstream absName;
      absName << "Surface_hole_" << name << "_abs";
      G4OpticalSurface* absSurface = new G4OpticalSurface(absName.str().c_str());
      absSurface->SetType(dielectric_metal);
      absSurface->SetFinish(ground);
      absSurface->SetModel(unified);
      absSurface->SetSigmaAlpha(0.1);
      absSurface->SetMaterialPropertiesTable(MyMaterials::ABS_SURF());
      G4LogicalBorderSurface* absLB;
      absName << "LB";
      absLB = new G4LogicalBorderSurface(absName.str().c_str(),                                                                                       holePV,
                             absorber1PV,
                             absSurface);





    } // loop on y direction

    ++countX ;
  } // loop on x direction






  //***********************8**    ******************   8******************
  // loop on x direction
  int countX3 = 0 ; // for the staggering
  int indexion4;
  int indexion3;
  for (float x = -0.5*Second_module_xy+Second_startX; countX3 < Second_nFibresAlongX; x += Second_fibreDistanceAlongX)
  {
    // loop on y direction
    int countY3 = 0 ; // for the staggering
    for (float y = - 0.5 * Second_module_yx + Second_startY; countY3 < Second_nFibresAlongY; y += Second_fibreDistanceAlongY)
    {
      float x_c = x;
      float y_c = y;
      //____________________________________________
      //   float x_cs = -0.5 * module_xy;
      //    float y_cs = -0.5 * module_xy;
      //________________________________________________________________
      // staggering
      if( (Second_fibre_scheme%2) == 0 )
      y_c += 0.5*Second_fibreDistanceAlongY*(countX3%2) ;

      int index = index1 + countX3 * Second_nFibresAlongY + countY3 ;
      int index2 = countX3 * Second_nFibresAlongY + countY3;
      CreateTree::Instance() -> fibresPosition -> Fill(index,x_c,y_c);
      CreateTree::Instance() -> fibresPosition_2nd_Section -> Fill(index2,x_c,y_c);

      std::string name;
      G4PVPlacement *fiberPV;
      G4PVPlacement *holePV;


      // **************  ***************
      if (x >= -0.5*Second_fibres_x1 && x<= 0.5*Second_fibres_x1 && y >= -0.5*Second_fibres_y1 && y <= 0.5*Second_fibres_y1)
      {
        name = Form("holePV %d",index);
      holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber2LV, false, 0, checkOverlaps) ;

      name = Form("fibre22SPV %d",index);
      fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre22SLV, name, absorber2LV, false, 0, checkOverlaps) ;


      // ************ *********************
    }

    // #11,#13,#15,#17
    else if (x >= -1.5*fibres_x1 && x<= fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1 || x >= fibres_x1 && x<= 1.5*fibres_x1 && y >= -0.5*fibres_y1 && y <= 0.5*fibres_y1 || x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= -1.5*fibres_y1 && y <= -0.5*fibres_y1 ||  x >= -0.5*fibres_x1 && x<= 0.5*fibres_x1 && y >= 0.5*fibres_y1 && y <= 1.5*fibres_y1 )
    {
      name = Form("holePV %d",index);
      holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber2LV, false, 0, checkOverlaps) ;

      name = Form("fibre2SPV %d",index);
      fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre2SLV, name, absorber2LV, false, 0, checkOverlaps) ;
    }

    else
    {
      name = Form("holePV %d",index);
      holePV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), holeLV, name, absorber2LV, false, 0, checkOverlaps) ;

      name = Form("fibre23SPV %d",index);
      fiberPV = new G4PVPlacement (0,G4ThreeVector(x_c,y_c,0.), fibre23SLV, name, absorber2LV, false, 0, checkOverlaps) ;
    }

    indexion4 = index+1;
    indexion3 = index2+1;
    ++countY3 ;

    //optical surface between fiber and thin air
    std::stringstream Surfname;
    Surfname << "Surface_" << name << "_airLayer";
    G4OpticalSurface* reflector_surf = new G4OpticalSurface(Surfname.str().c_str());
    reflector_surf->SetType(dielectric_metal);
    reflector_surf->SetFinish(polished);
    reflector_surf->SetModel(unified);
    reflector_surf->SetMaterialPropertiesTable(MyMaterials::ESR());
    G4LogicalBorderSurface* reflectorLB;
    Surfname << "LB";
    reflectorLB = new G4LogicalBorderSurface(Surfname.str().c_str(),
                                             fiberPV,
                                             airLayerPV,
                                             reflector_surf);
    //
    std::stringstream absName;
    absName << "Surface_hole_" << name << "_abs";
    G4OpticalSurface* absSurface = new G4OpticalSurface(absName.str().c_str());
    absSurface->SetType(dielectric_metal);
    absSurface->SetFinish(ground);
    absSurface->SetModel(unified);
    absSurface->SetSigmaAlpha(0.1);
    absSurface->SetMaterialPropertiesTable(MyMaterials::ABS_SURF());
    G4LogicalBorderSurface* absLB;
    absName << "LB";
    absLB = new G4LogicalBorderSurface(absName.str().c_str(),                                                                             holePV,
                           absorber2PV,
                           absSurface);
  } // loop on y direction


  ++countX3 ;
} // loop on x direction



// temporary
// build the cones


// G4LogicalVolume *subtractLV = new G4LogicalVolume (subtract,PLEXMaterial,"subtractLV");
// G4PVPlacement *subtractPV = new G4PVPlacement (0,G4ThreeVector(0,0,-40.0 *cm), subtractLV, "subtractPV", worldLV, false, 0, checkOverlaps) ;







//**********    ***************  ************ ************ ****************


// fibre gap for photon counting
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
G4VSolid * gapLayerS = NULL;
G4VSolid * gapS      = NULL;
G4LogicalVolume * gapLayerLV = NULL;
G4LogicalVolume * gapLV      = NULL;
if( postshower )
{
  gapLayerS = new G4Box ("gapLayerS", 0.5*module_xy, 0.5*module_yx, 0.5*depth) ;
  gapS      = new G4Box (     "gapS", 0.5*module_xy, 0.5*module_yx, 0.5*(gap_l-depth)) ;
  gapLayerLV = new G4LogicalVolume (gapLayerS, GaMaterial, "gapLayerLV") ;
  gapLV      = new G4LogicalVolume (gapS,      GaMaterial,      "gapLV") ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*module_z+0.5*depth+0.5*Second_module_z),          gapLayerLV, "gapLayerPV", worldLV, false, 0, checkOverlaps) ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*module_z+0.5*Second_module_z+depth+0.5*(gap_l-depth)), gapLV,      "gapPV", worldLV, false, 0, checkOverlaps) ;
}


// Si detector for photon counting
// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
G4VSolid * detLayerS = NULL;
G4VSolid * detS      = NULL;
G4LogicalVolume * detLayerLV = NULL;
G4LogicalVolume * detLV      = NULL;
if( postshower )
{
  detLayerS = new G4Box ("detLayerS", 0.5*module_xy, 0.5*module_yx, 0.5*depth) ;
  detS      = new G4Box (     "detS", 0.5*module_xy, 0.5*module_yx, 0.5*(det_l-depth)) ;
  detLayerLV = new G4LogicalVolume (detLayerS, DeMaterial, "detLayerLV") ;
  detLV      = new G4LogicalVolume (detS,      DeMaterial,      "detLV") ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*module_z + 0.5*Second_module_z +gap_l+0.5*depth),          detLayerLV, "detLayerPV", worldLV, false, 0, checkOverlaps) ;
  new G4PVPlacement (0, G4ThreeVector (0., 0., 0.5*module_z+0.5*Second_module_z+gap_l+depth+0.5*(det_l-depth)), detLV,      "detPV", worldLV, false, 0, checkOverlaps) ;
}


//-----------------------------------------------------
//------------- Visualization attributes --------------
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

// G4VisAttributes* Vissubtract = new G4VisAttributes (blue) ;
// Vissubtract->SetVisibility (true) ;
// Vissubtract->SetForceWireframe (true) ;
// coneLV->SetVisAttributes (Vissubtract) ;


G4VisAttributes* VisAttWorld = new G4VisAttributes (black) ;
VisAttWorld->SetVisibility (true) ;
VisAttWorld->SetForceWireframe (true) ;
worldLV->SetVisAttributes (VisAttWorld) ;

G4VisAttributes* VisAttCalorimeter = new G4VisAttributes (yellow) ;
VisAttCalorimeter->SetVisibility (true) ;
VisAttCalorimeter->SetForceWireframe (true) ;
calorimeterLV->SetVisAttributes (VisAttCalorimeter) ;


G4VisAttributes* VisAttGlue = new G4VisAttributes (green) ;
VisAttGlue->SetVisibility (true) ;
VisAttGlue->SetForceWireframe (true) ;
glue_LV->SetVisAttributes (VisAttGlue) ;


G4VisAttributes* VisAttPLEX = new G4VisAttributes (green) ;
VisAttPLEX->SetVisibility (true) ;
VisAttPLEX->SetForceWireframe (false) ;
LguideLV->SetVisAttributes (VisAttPLEX) ;

G4VisAttributes* VisAttPVC = new G4VisAttributes (brown) ;
VisAttPVC->SetVisibility (true) ;
VisAttPVC->SetForceWireframe (false) ;
PVCLV->SetVisAttributes (VisAttPVC) ;

G4VisAttributes* VisAttWires = new G4VisAttributes (blue) ;
VisAttWires->SetVisibility (true) ;
VisAttWires->SetForceWireframe (false) ;
WiresLV->SetVisAttributes (VisAttWires) ;

G4VisAttributes* VisAttAbsorber = new G4VisAttributes (grey) ;
VisAttAbsorber->SetVisibility (true) ;
VisAttAbsorber->SetForceWireframe (false) ;
absorber1LV->SetVisAttributes (VisAttAbsorber) ;

// ********* **************** *********

G4VisAttributes* VisAttAbsorber2 = new G4VisAttributes (blue) ;
VisAttAbsorber2->SetVisibility (true) ;
VisAttAbsorber2->SetForceWireframe (false) ;
absorber2LV->SetVisAttributes (VisAttAbsorber2) ;

G4VisAttributes* VisAttAirLayer = new G4VisAttributes (cyan) ;
VisAttAirLayer->SetVisibility (true) ;
VisAttAirLayer->SetForceWireframe (false) ;
airLayerLV->SetVisAttributes (VisAttAirLayer) ;
//
// G4VisAttributes* VisAttAirLayer = new G4VisAttributes (cyan) ;
// VisAttAirLayer->SetVisibility (true) ;
// VisAttAirLayer->SetForceWireframe (false) ;
// airLayerLV->SetVisAttributes (VisAttAirLayer) ;

G4VisAttributes* VisAttPVC_pmt = new G4VisAttributes (red) ;
VisAttPVC_pmt->SetVisibility (true) ;
VisAttPVC_pmt->SetForceWireframe (true) ;
PVC_pmt_backLV->SetVisAttributes (VisAttPVC_pmt) ;
PVC_pmt_frontLV->SetVisAttributes (VisAttPVC_pmt) ;

G4VisAttributes* VisAtt_air_LG = new G4VisAttributes (blue) ;
VisAtt_air_LG->SetVisibility (true) ;
VisAtt_air_LG->SetForceWireframe (true) ;
air_LG_back_LV->SetVisAttributes (VisAtt_air_LG) ;
air_LG_front_LV->SetVisAttributes (VisAtt_air_LG) ;



// ********* ************** * ***************   *********


bool wireFrame = true;
G4VisAttributes* VisAttHole = new G4VisAttributes(air);
VisAttHole->SetVisibility(true);
VisAttHole->SetForceWireframe(wireFrame);
holeLV->SetVisAttributes(VisAttHole);

G4VisAttributes* VisAttfibre1S = new G4VisAttributes (yellow) ;
VisAttfibre1S->SetVisibility (true) ;
VisAttfibre1S->SetForceWireframe (wireFrame) ;
fibre1SLV->SetVisAttributes (VisAttfibre1S) ;

//______________________________

G4VisAttributes* VisAttfibre12S = new G4VisAttributes (red) ;
VisAttfibre12S->SetVisibility (true) ;
VisAttfibre12S->SetForceWireframe (wireFrame) ;
fibre12SLV->SetVisAttributes (VisAttfibre12S) ;


G4VisAttributes* VisAttfibre13S = new G4VisAttributes (green) ;
VisAttfibre13S->SetVisibility (true) ;
VisAttfibre13S->SetForceWireframe (wireFrame) ;
fibre13SLV->SetVisAttributes (VisAttfibre13S) ;

//******** ***********  * ***************    ****************
G4VisAttributes* VisAttfibre2S = new G4VisAttributes (yellow) ;
VisAttfibre2S->SetVisibility (true) ;
VisAttfibre2S->SetForceWireframe (wireFrame) ;
fibre2SLV->SetVisAttributes (VisAttfibre2S) ;

//______________________________

G4VisAttributes* VisAttfibre22S = new G4VisAttributes (red) ;
VisAttfibre22S->SetVisibility (true) ;
VisAttfibre22S->SetForceWireframe (wireFrame) ;
fibre22SLV->SetVisAttributes (VisAttfibre22S) ;


G4VisAttributes* VisAttfibre23S = new G4VisAttributes (green) ;
VisAttfibre23S->SetVisibility (true) ;
VisAttfibre23S->SetForceWireframe (wireFrame) ;
fibre23SLV->SetVisAttributes (VisAttfibre23S) ;


// ********** *********    *****************   *******************
//________________________________________________


if( postshower )
{
  G4VisAttributes* VisAttGapLayer = new G4VisAttributes(red);
  VisAttGapLayer->SetVisibility(true);
  VisAttGapLayer->SetForceWireframe(false);
  gapLayerLV->SetVisAttributes(VisAttGapLayer);

  G4VisAttributes* VisAttGap = new G4VisAttributes(blue);
  VisAttGap->SetVisibility(true);
  VisAttGap->SetForceWireframe(false);
  gapLV->SetVisAttributes(VisAttGap);

  G4VisAttributes* VisAttDetLayer = new G4VisAttributes(red);
  VisAttDetLayer->SetVisibility(true);
  VisAttDetLayer->SetForceWireframe(false);
  detLayerLV->SetVisAttributes(VisAttDetLayer);

  G4VisAttributes* VisAttDet = new G4VisAttributes(gray);
  VisAttDet->SetVisibility(true);
  VisAttDet->SetForceWireframe(false);
  detLV->SetVisAttributes(VisAttDet);
}

//PG call the magnetic field initialisation
if (B_field_intensity > 0.1 * tesla) ConstructField () ;



//-----------------------------------------------
//------------- Fast photon timing --------------
//-----------------------------------------------

std::vector<std::pair<double,double> > rIndVecCore;
std::vector<std::pair<double,double> > rIndVecClad;
std::vector<std::pair<double,double> > rIndVecAir;
std::vector<std::pair<double,double> > rIndVecGap;

G4MaterialPropertyVector* mpVec;


mpVec = ClMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
{
  std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  rIndVecCore.push_back(dummy);
}

//_______________________________________________________


mpVec = ClSSMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
{
  std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  rIndVecCore.push_back(dummy);
}
//  ***************     ******************   *******************
mpVec = Cl3Material->GetMaterialPropertiesTable()->GetProperty("RINDEX");
for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
{
  std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  rIndVecCore.push_back(dummy);
}

//_______________________________________________________


mpVec = Cl4SSMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
{
  std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  rIndVecCore.push_back(dummy);
}




// ****************  ********************   **************************


//__________________________________________-=+++++++++_+++++++

mpVec = WoMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
{
  std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  std::pair<double,double> dummy2(mpVec->GetLowEdgeEnergy(it)/eV,fibre_cladRIndex);
  rIndVecAir.push_back(dummy);
  rIndVecClad.push_back(dummy2);
}

mpVec = GaMaterial->GetMaterialPropertiesTable()->GetProperty("RINDEX");
for(unsigned int it = 0; it < mpVec->GetVectorLength(); ++it)
{
  std::pair<double,double> dummy(mpVec->GetLowEdgeEnergy(it)/eV,(*mpVec)[it]);
  rIndVecGap.push_back(dummy);
}


fib = FiberInit(fibre_length,fibre_radius,CreateTree::Instance()->attLengths,rIndVecCore,rIndVecClad,rIndVecAir,rIndVecGap) ;



G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl ;
return worldPV ;}



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
  if      ( glue_interface == 0 ) GlueMaterial = WoMaterial;
  else if ( glue_interface == 1 ) GlueMaterial = MyMaterials::OpticalGrease();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid glue material specifier " << glue_interface << G4endl ;
    exit (-1) ;
  }
  G4cout << "Glue material: "<< GlueMaterial << G4endl ;

  ConeMaterial  = NULL;
  if      ( cone_material == 0 ) ConeMaterial = WoMaterial;
  else if ( cone_material == 1 ) ConeMaterial = MyMaterials::PLEX () ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid cone material specifier " << cone_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Cone material: "<< ConeMaterial << G4endl ;

  AbMaterial = NULL ;
  if      ( abs_material == 1 ) AbMaterial = MyMaterials::Brass () ;
  else if ( abs_material == 2 ) AbMaterial = MyMaterials::Tungsten () ;
  else if ( abs_material == 3 ) AbMaterial = MyMaterials::Lead () ;
  else if ( abs_material == 4 ) AbMaterial = MyMaterials::Iron () ;
  else if ( abs_material == 5 ) AbMaterial = MyMaterials::Aluminium () ;
  else if ( abs_material == 6 ) AbMaterial = MyMaterials::CopperTungstenAlloy(W_fraction) ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid absorber material specifier " << abs_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Ab. material: "<< AbMaterial << G4endl ;

  ///*************88    ********************8  ******************


  AbMaterial2 = NULL ;
  if      ( Second_abs_material == 1 ) AbMaterial2 = MyMaterials::Brass () ;
  else if ( Second_abs_material == 2 ) AbMaterial2 = MyMaterials::Tungsten () ;
  else if ( Second_abs_material == 3 ) AbMaterial2 = MyMaterials::Lead () ;
  else if ( Second_abs_material == 4 ) AbMaterial2 = MyMaterials::Iron () ;
  else if ( Second_abs_material == 5 ) AbMaterial2 = MyMaterials::Aluminium () ;
  else if ( Second_abs_material == 6 ) AbMaterial2 = MyMaterials::CopperTungstenAlloy(W_fraction) ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid absorber material specifier " << Second_abs_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Ab. material2: "<< AbMaterial2 << G4endl ;

  PLEXMaterial = NULL ;
  PLEXMaterial = MyMaterials::PLEX () ;

  G4cout << "PLEX. material: "<< PLEXMaterial << G4endl ;

  PVCMaterial = NULL ;
  PVCMaterial = MyMaterials::PVC () ;

  G4cout << "PVC. material: "<< PVCMaterial << G4endl ;

  PMTMaterial = NULL ;
  PMTMaterial = MyMaterials::CuAir () ;

  G4cout << "PMT. material: "<< PMTMaterial << G4endl ;

  WiresMaterial = NULL ;
  WiresMaterial = MyMaterials::Cu () ;

  G4cout << "Wires. material: "<< WiresMaterial << G4endl ;

  PlaneMaterial = NULL ;
  if      ( lp_mat == 1 ) PlaneMaterial = MyMaterials::Brass () ;
  else if ( lp_mat == 2 ) PlaneMaterial = MyMaterials::Tungsten () ;
  else if ( lp_mat == 3 ) PlaneMaterial = MyMaterials::Lead () ;
  else if ( lp_mat == 4 ) PlaneMaterial = MyMaterials::Iron () ;
  else if ( lp_mat == 5 ) PlaneMaterial = MyMaterials::Aluminium () ;
  else if ( lp_mat == 6 ) PlaneMaterial = MyMaterials::CopperTungstenAlloy(W_fraction) ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid absorber material specifier " << lp_mat << G4endl ;
    exit (-1) ;
  }
  G4cout << "Plane material: "<< PlaneMaterial << G4endl ;


  // ***********  ****************  ****************   ***********


  ClMaterial = NULL ;
  if      ( fibre_material == 1 ) ClMaterial = MyMaterials::Quartz () ;
  else if ( fibre_material == 2 ) ClMaterial = MyMaterials::SiO2_Ce () ;
  else if ( fibre_material == 3 ) ClMaterial = MyMaterials::DSB_Ce () ;
  else if ( fibre_material == 4 ) ClMaterial = MyMaterials::LuAG_Ce () ;
  else if ( fibre_material == 5 ) ClMaterial = MyMaterials::YAG_Ce () ;
  else if ( fibre_material == 6 ) ClMaterial = MyMaterials::GAGG_Ce_Mg() ;
  else if ( fibre_material == 7 ) ClMaterial = MyMaterials::Water() ;

  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << fibre_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Cl. material: "<< ClMaterial << G4endl ;


  //____________________________________________________________

  ClSSMaterial = NULL ;
  if      ( fibre_material1 == 1 ) ClSSMaterial = MyMaterials::Quartz () ;
  else if ( fibre_material1 == 2 ) ClSSMaterial = MyMaterials::SiO2_Ce () ;
  else if ( fibre_material1 == 3 ) ClSSMaterial = MyMaterials::DSB_Ce () ;
  else if ( fibre_material1 == 4 ) ClSSMaterial = MyMaterials::LuAG_Ce () ;
  else if ( fibre_material1 == 5 ) ClSSMaterial = MyMaterials::YAG_Ce () ;
  else if ( fibre_material1 == 6 ) ClSSMaterial = MyMaterials::GAGG_Ce_Mg() ;
  else if ( fibre_material1 == 7 ) ClSSMaterial = MyMaterials::Water() ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << fibre_material1 << G4endl ;
    exit (-1) ;
  }
  G4cout << "ClSS. material: "<< ClSSMaterial << G4endl ;

  ClSSSMaterial = NULL ;
  if      ( fibre_material2 == 1 ) ClSSSMaterial = MyMaterials::Quartz () ;
  else if ( fibre_material2 == 2 ) ClSSSMaterial = MyMaterials::SiO2_Ce () ;
  else if ( fibre_material2 == 3 ) ClSSSMaterial = MyMaterials::DSB_Ce () ;
  else if ( fibre_material2 == 4 ) ClSSSMaterial = MyMaterials::LuAG_Ce () ;
  else if ( fibre_material2 == 5 ) ClSSSMaterial = MyMaterials::YAG_Ce () ;
  else if ( fibre_material2 == 6 ) ClSSSMaterial = MyMaterials::GAGG_Ce_Mg() ;
  else if ( fibre_material2 == 7 ) ClSSSMaterial = MyMaterials::Water() ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << fibre_material1 << G4endl ;
    exit (-1) ;
  }
  G4cout << "2nd_Sect Corners fibres Material:"<< ClSSMaterial << G4endl ;


  // *************    ***************   *****************  * *******
  Cl3Material = NULL ;
  if      ( Second_fibre_material == 1 ) Cl3Material = MyMaterials::Quartz () ;
  else if ( Second_fibre_material == 2 ) Cl3Material = MyMaterials::SiO2_Ce () ;
  else if ( Second_fibre_material == 3 ) Cl3Material = MyMaterials::DSB_Ce () ;
  else if ( Second_fibre_material == 4 ) Cl3Material = MyMaterials::LuAG_Ce () ;
  else if ( Second_fibre_material == 5 ) Cl3Material = MyMaterials::YAG_Ce () ;
  else if ( Second_fibre_material == 6 ) Cl3Material = MyMaterials::GAGG_Ce_Mg() ;
  else if ( Second_fibre_material == 7 ) Cl3Material = MyMaterials::Water() ;

  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Second_fibre_material << G4endl ;
    exit (-1) ;
  }
  G4cout << "Cl. material: "<< Cl3Material << G4endl ;


  //____________________________________________________________

  Cl4SSMaterial = NULL ;
  if      ( Second_fibre_material1 == 1 ) Cl4SSMaterial = MyMaterials::Quartz () ;
  else if ( Second_fibre_material1 == 2 ) Cl4SSMaterial = MyMaterials::SiO2_Ce () ;
  else if ( Second_fibre_material1 == 3 ) Cl4SSMaterial = MyMaterials::DSB_Ce () ;
  else if ( Second_fibre_material1 == 4 ) Cl4SSMaterial = MyMaterials::LuAG_Ce () ;
  else if ( Second_fibre_material1 == 5 ) Cl4SSMaterial = MyMaterials::YAG_Ce () ;
  else if ( Second_fibre_material1 == 6 ) Cl4SSMaterial = MyMaterials::GAGG_Ce_Mg() ;
  else if ( Second_fibre_material1 == 7 ) Cl4SSMaterial = MyMaterials::Water() ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Second_fibre_material1 << G4endl ;
    exit (-1) ;
  }
  G4cout << "ClSS. material: "<< Cl4SSMaterial << G4endl ;


  Cl43SSMaterial = NULL ;
  if      ( Second_fibre_material2 == 1 ) Cl43SSMaterial = MyMaterials::Quartz () ;
  else if ( Second_fibre_material2 == 2 ) Cl43SSMaterial = MyMaterials::SiO2_Ce () ;
  else if ( Second_fibre_material2 == 3 ) Cl43SSMaterial = MyMaterials::DSB_Ce () ;
  else if ( Second_fibre_material2 == 4 ) Cl43SSMaterial = MyMaterials::LuAG_Ce () ;
  else if ( Second_fibre_material2 == 5 ) Cl43SSMaterial = MyMaterials::YAG_Ce () ;
  else if ( Second_fibre_material2 == 6 ) Cl43SSMaterial = MyMaterials::GAGG_Ce_Mg() ;
  else if ( Second_fibre_material2 == 7 ) Cl43SSMaterial = MyMaterials::Water() ;
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Second_fibre_material1 << G4endl ;
    exit (-1) ;
  }
  G4cout << "2nd_Sect Corners fibres Material: "<< Cl4SSMaterial << G4endl ;





  //  *****************   ********************   ***************


  //______________________________________________________________

  GaMaterial = NULL;
  if     ( gap_material == 1 ) GaMaterial = MyMaterials::Air();
  else if( gap_material == 2 ) GaMaterial = MyMaterials::OpticalGrease();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid gap material specifier " << gap_material << G4endl;
    exit(-1);
  }
  G4cout << "Gap material: " << gap_material << G4endl;


  DeMaterial = NULL;
  if( det_material == 1 ) DeMaterial = MyMaterials::Silicon();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid detector material specifier " << det_material << G4endl;
    exit(-1);
  }
  G4cout << "Detector material: " << det_material << G4endl;



  if( fibre_absLength >= 0 )
  {
    const G4int nEntries_ABS = 2;
    G4double PhotonEnergy_ABS[nEntries_ABS] = { 1.*eV, 10.*eV };
    G4double Absorption[nEntries_ABS] = { fibre_absLength*mm, fibre_absLength*mm };

    ClMaterial->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
    ClMaterial->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH",PhotonEnergy_ABS,Absorption,nEntries_ABS);
    //_______________________________

    ClSSMaterial->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
    ClSSMaterial->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH",PhotonEnergy_ABS,Absorption,nEntries_ABS);
    //___________________________________________________________
  }

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
