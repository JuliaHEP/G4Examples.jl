#---RoundCube definition---------------------------------------------------------------------------
# A RoundCube is a cube with rounded edges and vertices. The cube is defined by the side length a
# and the radius of the rounded edges and vertices r. 
function createRoundCube(a, r)
    ca = a/4
    cube = G4Box("cube",ca,ca,ca)
    cyl  = G4Tubs("cyl",0, r, ca,0, 2π)
    orb = G4Orb("orb", r)

    edge = G4Box("edge",r, r, ca)
    vert = G4Box("vert",r, r, r)

    acyl1 = G4SubtractionSolid("ucyl", CxxPtr(edge), CxxPtr(cyl), G4Transform3D(G4RotationMatrix(), G4ThreeVector(-r,-r,0)))
    acyl2 = G4SubtractionSolid("ucyl", CxxPtr(edge), CxxPtr(cyl), G4Transform3D(G4RotationMatrix(), G4ThreeVector(-r,r,0)))
    acyl3 = G4SubtractionSolid("ucyl", CxxPtr(edge), CxxPtr(cyl), G4Transform3D(G4RotationMatrix(), G4ThreeVector(r,r,0)))
    aorb = G4SubtractionSolid("uorb", CxxPtr(vert), CxxPtr(orb), G4Transform3D(G4RotationMatrix(), G4ThreeVector(-r,-r, -r)))

    t1 = G4Transform3D(G4RotationMatrix(), G4ThreeVector(ca,ca,0))
    t2 = G4Transform3D(G4RotationMatrix(0,π/2,0), G4ThreeVector(ca,0,ca))
    t3 = G4Transform3D(G4RotationMatrix(0,π/2,π/2), G4ThreeVector(0,ca,ca))
    t4 = G4Transform3D(G4RotationMatrix(), G4ThreeVector(ca,ca,ca))

    s1 = G4SubtractionSolid("s1", CxxPtr(cube), CxxPtr(acyl1), t1)
    s2 = G4SubtractionSolid("s2", CxxPtr(s1), CxxPtr(acyl2), t2)
    s3 = G4SubtractionSolid("s3", CxxPtr(s2), CxxPtr(acyl3), t3)
    s4 = G4SubtractionSolid("s4", CxxPtr(s3), CxxPtr(aorb), t4)

    t5 = G4Transform3D(G4RotationMatrix(0,0,0), G4ThreeVector(ca,ca,ca))
    t6 = G4Transform3D(G4RotationMatrix(0,π/2,0), G4ThreeVector(ca,ca,-ca))
    t7 = G4Transform3D(G4RotationMatrix(0,-π/2,0), G4ThreeVector(ca,-ca,ca))
    t8 = G4Transform3D(G4RotationMatrix(0,π/2,π/2), G4ThreeVector(ca,-ca,-ca))

    u1 = G4UnionSolid("u1", CxxPtr(cube), CxxPtr(s4), t5)
    u2 = G4UnionSolid("u2", CxxPtr(u1), CxxPtr(s4), t6)
    u3 = G4UnionSolid("u3", CxxPtr(u2), CxxPtr(s4), t7)
    u4 = G4UnionSolid("u4", CxxPtr(u3), CxxPtr(s4), t8)

    t9 = G4Transform3D(G4RotationMatrix(0,0,π), G4ThreeVector())
    u5 = G4UnionSolid("u5", CxxPtr(u4), CxxPtr(u4), t9)
    
    # remove finalizers
    for a in (cube, cyl, orb, edge, vert, acyl1, acyl2, acyl3, aorb, s1, s2, s3, s4, u1, u2, u3, u4)
        a.cpp_object = C_NULL
    end
    return u5
end
