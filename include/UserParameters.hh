  G4int ZdivisionStereoLayers;
  G4int NStripsIn;
  G4int NStripsMid;
  G4int NStripsOut;
  
  G4double DetectorLength;
  G4double DetectorInnerRadius;
  G4double LayerInWidth;
  G4double LayerMidWidth;
  G4double LayerOutWidth;
  G4double StereoAngleInner;
  G4double StereoAngleOuter;

  
  

  //--------- SOURCE PARAMETERS, AND CORRESPONDING PORJECTION
  G4int ProjectionMatrixType; /*0 for Cartesian, 1 for Polar */
  G4double SourceLengthScale; /* Radius or Length (in rectangular source)*/
  G4double ProjectionPixelI; /* I-partitions*/
  G4double ProjectionPixelJ; /* J-partitions*/
 
  G4double World_sizeX;
  G4double World_sizeY;
  G4double World_sizeZ;