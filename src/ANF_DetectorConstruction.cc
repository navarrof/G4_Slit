#include "ANF_DetectorConstruction.hh"
#include "G4GenericTrap.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Polyhedra.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ANF_DetectorConstruction::ANF_DetectorConstruction(): G4VUserDetectorConstruction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ANF_DetectorConstruction::~ANF_DetectorConstruction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ANF_DetectorConstruction::Construct()
{  

    G4double SlitSeparation = 100*CLHEP::um;
    // ----------------- Define World ----------------- // 

    G4double zet = 1.0;
    G4double amass = 1.01*CLHEP::g/CLHEP::mole;
    G4double density = CLHEP::universe_mean_density;
    G4double pressure = 3.0e-18*CLHEP::pascal;
    G4double temperature = 2.73*CLHEP::kelvin;
    G4Material * WorldMaterial = new G4Material("Galactic", zet, amass, density, 
                                                kStateGas, temperature, pressure);

    G4double worldXsize = 1.0*CLHEP::m;
    G4double worldYsize = 1.0*CLHEP::m;
    G4double worldZsize = 1.0*CLHEP::m;

    G4Box* worldSolid = new G4Box( "solid-world",       // name
                                    0.5*worldXsize,     // box half x-size
                                    0.5*worldYsize,     // box half y-size
                                    0.5*worldZsize     // box half z-size
                                     );

    G4LogicalVolume* worldLogical = new G4LogicalVolume( worldSolid,
                                                          WorldMaterial,
                                                          "logic-world" );

    G4VPhysicalVolume* worldPhysical = new G4PVPlacement( nullptr, 
                                                        G4ThreeVector(0.0,0.0,0.0),
                                                        worldLogical,
                                                        "world",
                                                        nullptr,
                                                        false, 
                                                        0.0);

    // -------------------------- Copper Base ---------------------------- //
    G4double CopperX =  64*CLHEP::mm;
    G4double CopperY = 64*CLHEP::mm; 
    G4double CopperZ = 19.9*CLHEP::mm;

    G4double A_cu  =  63.85*CLHEP::g/CLHEP::mole; 
    G4double density_cu = 8.96*CLHEP::g/CLHEP::cm3;
    G4Element* elCu   =  new G4Element("Copper","Cu",29.,A_cu);
    G4Material* CopperMaterial = new G4Material("Base_Copper",density_cu,1);
    CopperMaterial->AddElement(elCu, 1.0);

    G4Box* CopperSolidBase = new G4Box( "Copper-cube",
                                    CopperX/2.0,              
                                    CopperY/2.0,             
                                    CopperZ/2.);

    G4LogicalVolume* CopperLogicalBase = new G4LogicalVolume(CopperSolidBase,
                                                        CopperMaterial,
                                                        "CopperBase-Logical");

    G4VPhysicalVolume* CopperPhysical1 = new G4PVPlacement(nullptr,
                                                            G4ThreeVector(-(SlitSeparation+CopperY)/2.0,0.0,CopperZ/2.0),
                                                            CopperLogicalBase,
                                                            "CopperBase-Physical1",
                                                            worldLogical,
                                                            false, 
                                                            0.0);

    G4VPhysicalVolume* CopperPhysical2 = new G4PVPlacement(nullptr,
                                                            G4ThreeVector((SlitSeparation+CopperY)/2.0,0.0,CopperZ/2.0),
                                                            CopperLogicalBase,
                                                            "CopperBase-Physical1",
                                                            worldLogical,
                                                            false, 
                                                            0.0);
    
    static const double     pi  = 3.14159265358979323846;
    G4double TriangleX = 25.56*CLHEP::mm;
    G4double TriangleY = 64*CLHEP::mm;
    G4double TriangleZ = 67.39-CopperZ;

    G4RotationMatrix* rotationMatrix = new G4RotationMatrix();
    rotationMatrix->rotateX(180.*deg); 

    std::vector<G4TwoVector> polygon(8);
    polygon[0].set(-TriangleX/2.0,-TriangleY/2.0);
    polygon[1].set(-TriangleX/2.0,TriangleY/2.0); 
    polygon[2].set(TriangleX/2.0,TriangleY/2.0); 
    polygon[3].set(TriangleX/2.0,-TriangleY/2.0); 
    polygon[4].set(0,-TriangleY/2.0); 
    polygon[5].set(0,TriangleY/2.0);
    polygon[6].set(0,TriangleY/2.0);   // It is ok!
    polygon[7].set(0,-TriangleY/2.0);

    G4double SepTriang = 5.28*CLHEP::mm;
    G4GenericTrap* CopperPoly = new G4GenericTrap("Copper Top", TriangleZ,polygon);

    G4LogicalVolume* CopperLogicalTop = new G4LogicalVolume(CopperPoly,
                                                        CopperMaterial,
                                                        "CopperTop-Logical");

    
    G4VPhysicalVolume* CopperPhysicalTop1 = new G4PVPlacement(rotationMatrix,
                                                            G4ThreeVector(SlitSeparation/2+3/2.*SepTriang+3/2.0*TriangleX,0.0,-TriangleZ),
                                                            CopperLogicalTop,
                                                            "CopperTop-Physical1",
                                                            worldLogical,
                                                            false, 
                                                            0.0);

    G4VPhysicalVolume* CopperPhysicalTop2 = new G4PVPlacement(rotationMatrix,
                                                            G4ThreeVector(1.0/2.0*SepTriang+SlitSeparation/2+1/2.0*TriangleX,0.0,-TriangleZ),
                                                            CopperLogicalTop,
                                                            "CopperTop-Physical2",
                                                            worldLogical,
                                                            false, 
                                                            0.0);

    G4VPhysicalVolume* CopperPhysicalTop3 = new G4PVPlacement(rotationMatrix,
                                                            G4ThreeVector(-1.0/2.0*SepTriang-SlitSeparation/2-1/2.0*TriangleX,0.0,-TriangleZ),
                                                            CopperLogicalTop,
                                                            "CopperTop-Physical3",
                                                            worldLogical,
                                                            false, 
                                                            0.0);

    G4VPhysicalVolume* CopperPhysicalTop4 = new G4PVPlacement(rotationMatrix,
                                                            G4ThreeVector(-SlitSeparation/2-3/2.0*SepTriang-3.0/2.0*TriangleX,0.0,-TriangleZ),
                                                            CopperLogicalTop,
                                                            "CopperTop-Physical4",
                                                            worldLogical,
                                                            false, 
                                                            0.0);

    // ----------------------------- CARBON TOP --------------------------- //
    G4double A_C  =  12.00*CLHEP::g/CLHEP::mole; 
    G4double density_C = 2.265*CLHEP::g/CLHEP::cm3;
    G4Element* elC   =  new G4Element("Carbon","C",6.,A_C);
    G4Material* GraphiteMaterial = new G4Material("Top_Graphite",density_C,1);
    GraphiteMaterial->AddElement(elC, 1.0);

    G4double TriangleGrX = 25.56*CLHEP::mm + 7*CLHEP::mm;
    G4double TriangleGrY = 64*CLHEP::mm + 7*CLHEP::mm;
    G4double TriangleGrZ = 67.39-CopperZ + 7*CLHEP::mm;

    G4RotationMatrix* rotationMatrix2 = new G4RotationMatrix();
    rotationMatrix2->rotateX(180.*deg); 

    std::vector<G4TwoVector> polygon2(8);
    polygon2[0].set(-TriangleGrX/2.0,-TriangleGrY/2.0);
    polygon2[1].set(-TriangleGrX/2.0,TriangleGrY/2.0); 
    polygon2[2].set(TriangleGrX/2.0,TriangleGrY/2.0); 
    polygon2[3].set(TriangleGrX/2.0,-TriangleGrY/2.0); 
    polygon2[4].set(0,-TriangleGrY/2.0); 
    polygon2[5].set(0,TriangleGrY/2.0);
    polygon2[6].set(0,TriangleGrY/2.0);   // It is ok!
    polygon2[7].set(0,-TriangleGrY/2.0);

    G4double SepTriang2 = 5.28*CLHEP::mm;
    G4GenericTrap* GraphPoly = new G4GenericTrap("Carbon Top", TriangleGrZ,polygon2);

    G4SubtractionSolid* GraphTop_Solid = new G4SubtractionSolid("solid-CarbonTop",GraphPoly,CopperPoly,0,G4ThreeVector(0.0,0.0,0.0));
    G4LogicalVolume* GraphLogic = new G4LogicalVolume(GraphTop_Solid,GraphiteMaterial,"Graphite-Logical");

    
    G4VPhysicalVolume* GraphitePhysical1 = new G4PVPlacement(rotationMatrix2,
                                                            G4ThreeVector(SlitSeparation/2+3/2.*SepTriang+3/2.0*TriangleX,0.0,-TriangleZ),
                                                            GraphLogic,
                                                            "Graphite-Physical1",
                                                            worldLogical,
                                                            false, 
                                                            0.0);

    G4VPhysicalVolume* GraphitePhysical2 = new G4PVPlacement(rotationMatrix2,
                                                            G4ThreeVector(1.0/2.0*SepTriang+SlitSeparation/2+1/2.0*TriangleX,0.0,-TriangleZ),
                                                            GraphLogic,
                                                            "Graphite-Physical2",
                                                            worldLogical,
                                                            false, 
                                                            0.0);

    G4VPhysicalVolume* GraphitePhysical3 = new G4PVPlacement(rotationMatrix2,
                                                            G4ThreeVector(-1.0/2.0*SepTriang-SlitSeparation/2-1/2.0*TriangleX,0.0,-TriangleZ),
                                                            GraphLogic,
                                                            "Graphite-Physical3",
                                                            worldLogical,
                                                            false, 
                                                            0.0);

    G4VPhysicalVolume* GraphitePhysical4 = new G4PVPlacement(rotationMatrix2,
                                                            G4ThreeVector(-SlitSeparation/2-3/2.0*SepTriang-3.0/2.0*TriangleX,0.0,-TriangleZ),
                                                            GraphLogic,
                                                            "Graphite-Physical4",
                                                            worldLogical,
                                                            false, 
                                                            0.0);
    
    
  return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
