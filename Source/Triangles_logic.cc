
for (int i = 0; i < iaxis - 1; i++)
{
    for (int j = 0; j < jaxis - 1; j++)
    {
        
        // Makes boxes between four relief data points (above homogeneous ground),
        // takes height to be the lowest of the four (RMins[i][j])
        if (RMins[i][j] > 0.0*m)
        {
            G4VSolid* RBox = new G4Box("RBox", Pnt_Sep/2, Pnt_Sep/2, RMins[i][j]/2);
            G4LogicalVolume* RBoxLogicVolume = new G4LogicalVolume(RBox, GroundMat, "RBox");
            new G4PVPlacement(0,               // no rotation
                    G4ThreeVector((1-iaxis) * Pnt_Sep /2 + Pnt_Sep * i,(1-jaxis) * Pnt_Sep /2 + Pnt_Sep * j, DetectorDepth/2 + RMins[i][j]/2),
                    RBoxLogicVolume,
                    "RBox",         // its name
                    logicWorld,  // its mother  volume
                    false,           // no boolean operations
                    0);
        }
        
        
        // Makes surface continuous by placing tetraheadra on top of each box
        auto make_tets = [&](int anchorx, int anchory, double signx, double signy)
        {
            // anchorx, anchory - indices of the anchor point (minima of the adjacent data points)
            // signx, signy - (anchorx, anchory > i,j) ==> +1.0, else -1.0
            //              - relative direction of the other data points
            
            
            // make a tetrahedron and place it at the anchor point
            auto make_a_tet = [&](G4ThreeVector p1, G4ThreeVector p2, G4ThreeVector p3)
            {
                // p[i] - points of the tetrahedron relative to anchor
                G4VSolid* RTet = new G4Tet("RTet", G4ThreeVector {0, 0, 0}, p1, p2, p3, 0);
                G4LogicalVolume* RTetLogicVolume = new G4LogicalVolume(RTet, GroundMat, "RTet");
                new G4PVPlacement(0,               // no rotation
                        G4ThreeVector((2 * anchorx - iaxis) * Pnt_Sep /2, (2 * anchory - jaxis) * Pnt_Sep /2, DetectorDepth/2 + RMins[i][j]),
                        RTetLogicVolume,
                        "RTet",         // its name
                        logicWorld,  // its mother  volume
                        false,           // no boolean operations
                        0);
            };
            
            int n = 0;
            
            // x adjascent
            if  (Relief[i+signx][j] != RMins[i][j])
            {
                G4ThreeVector P1 = G4ThreeVector(signx * Pnt_Sep, 0, 0);
                G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);
                G4ThreeVector P3 = G4ThreeVector(signx * Pnt_Sep, 0, Relief[i+signx][j]-RMins[i][j]);
                
                //make_a_tet(P1, P2, P3);
                
                n++;
            }
            // y adjascent
            if  (Relief[i][j+signy] != RMins[i][j])
            {
                G4ThreeVector P1 = G4ThreeVector(0,signy * Pnt_Sep, 0);
                G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);
                G4ThreeVector P3 = G4ThreeVector(0, signy * Pnt_Sep, Relief[i][j+signy]-RMins[i][j]);
                
                //make_a_tet(P1, P2, P3);
                
                n++;
            }
            // across
            if  (Relief[i+signy][j+signx] != RMins[i][j])
            {
                // across and x adjascent
                G4ThreeVector P1 = G4ThreeVector(signx * Pnt_Sep, 0, 0);
                G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, 0);
                G4ThreeVector P3 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);
                
                //make_a_tet(P1, P2, P3);
                
                // across and y adjascent
                G4ThreeVector P4 = G4ThreeVector(0,signy * Pnt_Sep, 0);
                G4ThreeVector P5 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, 0);
                G4ThreeVector P6 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);
                
                //make_a_tet(P4, P5, P6);
                
                n++;
            }
            // cap
            if (n > 0)
            {
                G4ThreeVector P1 = G4ThreeVector(0, signy * Pnt_Sep, Relief[i][j+signy]-RMins[i][j]);
                G4ThreeVector P2 = G4ThreeVector(signx * Pnt_Sep, 0, Relief[i+signx][j]-RMins[i][j]);
                G4ThreeVector P3 = G4ThreeVector(signx * Pnt_Sep, signy * Pnt_Sep, Relief[i+signx][j+signy]-RMins[i][j]);
                
                //make_a_tet(P1, P2, P3);
            }
        };
        
        
        // Finds lowest data point to anchor tetrahedra
        if (Relief[i][j] == RMins[i][j])
        {
            int Anchorx = i; // ----remove anchor varaibles
            int Anchory = j;
            int Signx = 1.0;
            int Signy = 1.0;
            
            make_tets(Anchorx, Anchory, Signx, Signy);
        }
        else if (Relief[i+1][j] == RMins[i][j])
        {
            int Anchorx = i + 1;
            int Anchory = j;
            int Signx = - 1.0;
            int Signy = 1.0;
            
            make_tets(Anchorx, Anchory, Signx, Signy);
        }
        else if (Relief[i][j+1] == RMins[i][j])
        {
            int Anchorx = i;
            int Anchory = j + 1;
            int Signx = 1.0;
            int Signy = - 1.0;
            
            make_tets(Anchorx, Anchory, Signx, Signy);
        }
        else
        {
            int Anchorx = i + 1;
            int Anchory = j + 1;
            int Signx = - 1.0;
            int Signy = - 1.0;
            
            make_tets(Anchorx, Anchory, Signx, Signy);
        }
    }
}

