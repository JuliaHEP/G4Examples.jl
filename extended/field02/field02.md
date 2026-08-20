```@meta
EditURL = "field02.lit"
```

# Extended/field02 Example
Test for investigation of tracking in electric field and field dependent electromagnetic processes.

The geometry consists of an "Absorber" that is a solid made of a given material.
Three parameters define the absorber :
- the material of the absorber,
- the thickness of an absorber,
- the radius of the absorber (the input face is a circle).

The primary kinematic consists of a single particle which hits the
absorber perpendicular to the input face. The type of the particle
and its energy are set

!!! note "Note that"
    You can also download this example as a
    [Jupyter notebook](field02.ipynb) and a plain
    [Julia source file](field02.jl).

#### Table of contents
```@contents
Pages = ["field02.md"]
Depth = 2:3
```

## Loading the necessary Julia modules
Load the `Geant4`, `Geant4.PhysicalConstants` and `Geant4.SystemOfUnits` modules. In addition we will use the `Parameters` module to handle the parameters of the detector.
We will also use the `FHist` and `Plots` modules to handle the histograms and plots.

````julia
using Geant4
using Geant4.SystemOfUnits
using Geant4.SystemOfUnits: volt
using Parameters
````

## Define the Detector
The Field02 detector is a simple cylinder filled with a given material. The detector is defined by the `Calorimeter` structure.

````julia
@with_kw mutable struct Calorimeter <: G4JLDetector
    # main input parameters
    material::String = "Kr20CO2"
    AbsorberThickness::Float64= 4cm
    AbsorberRadius::Float64 = 10cm
    AbsorberZ::Float64 = 36cm
    WorldSizeR::Float64 = 20cm
    WorldSizeZ::Float64 = 80cm
end

function Field02Construct(det::Calorimeter)::CxxPtr{G4VPhysicalVolume}
    (; material, AbsorberThickness, AbsorberRadius, AbsorberZ, WorldSizeR, WorldSizeZ) = det

    ##---Materials----------------------------------------------------------------------------------
    nist = G4NistManager!Instance()
    m_air = FindOrBuildMaterial(nist, "G4_AIR")

    C  = FindOrBuildElement(nist,6)
    O  = FindOrBuildElement(nist,8)
    Xe = FindOrBuildMaterial(nist, "G4_Xe")
    Kr = FindOrBuildMaterial(nist, "G4_Kr")

    _C0₂ = G4Material("CO2", 1.842 * mg / cm3, 2)
    AddElement(_C0₂, C, natoms=1)
    AddElement(_C0₂, O, natoms=2)
    CO₂ = FindOrBuildMaterial(nist, "CO2")

    Kr20CO2 = G4Material("Kr20CO2", 3.601 * mg / cm3, 2)
    AddMaterial(Kr20CO2, Kr, 0.89)
    AddMaterial(Kr20CO2, CO₂, 0.11)

    ##---Absorber----------------------------------------------------------------------------------
    m_absorber = FindOrBuildMaterial(nist, material)

    ##---Volumes------------------------------------------------------------------------------------
    worldS  = G4Tubs("world", 0, WorldSizeR, WorldSizeZ/2, 0, 360deg)
    worldLV = G4LogicalVolume(worldS, m_air, "World")
    worldPV = G4PVPlacement(nothing, G4ThreeVector(), worldLV, "World", nothing, false, 0, false)

    absorberS  = G4Tubs("absorber", 0, AbsorberRadius, AbsorberThickness/2, 0, 360deg)
    absorberLV = G4LogicalVolume(absorberS, m_absorber, "absorber")
    G4PVPlacement(nothing, G4ThreeVector(0, 0, AbsorberZ), absorberLV, "absorber", worldLV, false, 0, false)
    ##---Visualization attributes-------------------------------------------------------------------
    boxVisAtt = G4VisAttributes(G4Colour(1.0, 1.0, 1.0, 0.0))
    absorberVisAtt = G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.1))
    SetVisAttributes(worldLV, boxVisAtt)
    SetVisAttributes(absorberLV, absorberVisAtt)

    return worldPV
end
Geant4.getConstructor(::Calorimeter)::Function = Field02Construct
````

## Instantiate the detector with the parameters

````julia
det = Calorimeter(material="Kr20CO2",
                  AbsorberThickness=4cm,
                  AbsorberRadius=10cm,
                  AbsorberZ=0cm,
                  WorldSizeR=20cm,
                  WorldSizeZ=80cm)
````

````
Main.var"##277".Calorimeter
  material: String "Kr20CO2"
  AbsorberThickness: Float64 40.0
  AbsorberRadius: Float64 100.0
  AbsorberZ: Float64 0.0
  WorldSizeR: Float64 200.0
  WorldSizeZ: Float64 800.0

````

## Create the Primary Particle Generator

````julia
particlegun = G4JLGunGenerator(particle = "e-",
                               energy = 50MeV,
                               direction = G4ThreeVector(0,0,1),
                               position  = G4ThreeVector(0,0,0))
````

````
Geant4.G4JLGunGenerator("ParticleGun", Geant4.G4JLParticleGunData(nothing, "e-", G4ThreeVector(0.0,0.0,1.0), G4ThreeVector(0.0,0.0,0.0), 50.0), Geant4.var"#init#19"(), Geant4.var"#gen#20"(), Geant4.G4JLGeneratorAction[])
````

## Create Magnetic Field

````julia
efield = G4UniformElectricField(G4ThreeVector(0, 1e8volt/cm, 0))
````

````
Geant4.G4UniformElectricFieldAllocated(Ptr{Nothing}(0x0000000b8df6bf40))
````

## Create the Application

````julia
app = G4JLApplication(; detector = det,                               # detector with parameters
                        generator = particlegun,                      # primary particle generator
                        field = efield,                               # uniform magnetic field
                        nthreads = 0,                                 # # of threads (0 = no MT)
                        physics_type = FTFP_BERT,                     # what physics list to instantiate
                      );


configure(app)
initialize(app)
````

````

**************************************************************
 Geant4 version Name: geant4-11-04-patch-01 [MT]   (13-March-2026)
                       Copyright : Geant4 Collaboration
                      References : NIM A 506 (2003), 250-303
                                 : IEEE-TNS 53 (2006), 270-278
                                 : NIM A 835 (2016), 186-225
                             WWW : http://geant4.org/
**************************************************************


````

## Run the Application for 1 event with verbose tracking

````julia
ui`/tracking/verbose 1`
beamOn(app,1)
````

````

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 1,   Parent ID = 0
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0        0        0        0        50        0        0         0    absorber initStep
    1   0.0991    -15.3     10.5       203  0.00661     19.7      19.7    absorber eIoni
    2    0.243    -45.7     15.1       507   0.0151     30.7      50.5    absorber eIoni
    3    0.302    -68.6       17       736   0.0105       23      73.5    absorber eIoni
    4    0.361     -100     18.8  1.05e+03   0.0136     31.5       105       World Transportation
    5    0.364     -101     18.8  1.06e+03 0.000157     1.32       106       World eIoni
    6    0.378     -110     19.2  1.15e+03  0.00135     8.46       115       World eIoni
    7    0.444     -157       21  1.62e+03   0.0079     47.5       162       World eIoni
    8    0.469     -178     21.6  1.83e+03  0.00408     20.6       183       World eIoni
    9    0.493     -200     22.1  2.05e+03   0.0091     22.2       205  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 8,   Parent ID = 1
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    0.469     -178     21.6  0.000998        0        0         0       World initStep
    1    0.471     -178     21.6     0.214 0.000998   0.0569    0.0569       World eIoni
    2    0.485     -200     21.5       222  0.00548     22.7      22.8  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 7,   Parent ID = 1
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    0.444     -157       21   0.00212        0        0         0       World initStep
    1     0.44     -157       21     0.762  0.00212    0.178     0.178       World eIoni
    2   -0.418     -182     22.2       246  0.00425     24.6      24.8       World eIoni
    3   -0.727     -200     24.1       427  0.00275     18.3        43  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 6,   Parent ID = 1
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    0.378     -110     19.2   0.00121        0        0         0       World initStep
    1    0.379     -110     19.2     0.296  0.00121   0.0755    0.0755       World eIoni
    2    0.385     -200     19.2       902   0.0177       96      96.1  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 5,   Parent ID = 1
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    0.364     -101     18.8   0.00116        0        0         0       World initStep
    1     0.36     -101     18.8     0.276  0.00116   0.0711    0.0711       World eIoni
    2    0.337     -200     18.8       987     0.02      107       107  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 4,   Parent ID = 1
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    0.302    -68.6       17    0.0044        0        0         0    absorber initStep
    1    0.309    -68.7       17     0.596   0.0044    0.375     0.375    absorber eIoni
    2    0.318    -68.8       17      2.35 0.000777    0.175      0.55    absorber eIoni
    3    0.342     -100     16.9       314  0.00658     31.5        32       World Transportation
    4     5.98     -178     22.9  1.09e+03   0.0142     78.5       111       World eIoni
    5      6.8     -200     23.8  1.31e+03  0.00472     21.9       132  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 9,   Parent ID = 4
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0     5.98     -178     22.9   0.00136        0        0         0       World initStep
    1     5.98     -178     22.9     0.356  0.00136    0.089     0.089       World eIoni
    2     5.97     -200     22.9       219  0.00137     22.1      22.2  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 3,   Parent ID = 1
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    0.243    -45.7     15.1   0.00204        0        0         0    absorber initStep
    1     0.24    -45.7     15.1     0.186  0.00204    0.115     0.115    absorber eIoni
    2   -0.513    -49.1     14.7      34.7 0.000964     3.63      3.75    absorber eIoni
    3    -2.66    -90.7     14.2       450   0.0185     41.7      45.4    absorber eIoni
    4    -2.73    -95.2     14.2       495  0.00306     4.53      49.9    absorber eIoni
    5     -2.8    -99.8     14.2       541  0.00207     4.62      54.5    absorber eIoni
    6     -2.8     -100     14.2       543        0    0.146      54.7       World Transportation
    7    -2.99     -116     14.1       705   0.0058     16.2      70.9       World eIoni
    8    -3.56     -200       14  1.54e+03   0.0145     83.8       155  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 13,   Parent ID = 3
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    -2.99     -116     14.1   0.00358        0        0         0       World initStep
    1    -2.99     -116     14.1    0.0836 3.11e-05   0.0092    0.0092       World eIoni
    2    0.187     -139     23.8       232    0.011       27      27.1       World eIoni
    3    0.828     -192     33.8       755   0.0117     53.3      80.4       World eIoni
    4    0.886     -200     34.6       838   0.0017     8.38      88.8  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 15,   Parent ID = 13
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    0.828     -192     33.8   0.00103        0        0         0       World initStep
    1     0.83     -192     33.8     0.224  0.00103   0.0593    0.0593       World eIoni
    2    0.842     -200     33.7      83.4  0.00223     8.39      8.44  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 14,   Parent ID = 13
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    -2.99     -116     14.1    0.0012        0        0         0       World initStep
    1       -3     -116     14.1     0.288   0.0012   0.0738    0.0738       World eIoni
    2       -3     -200     14.1       838   0.0164       89        89  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 12,   Parent ID = 3
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0     -2.8    -99.8     14.2   0.00146        0        0         0    absorber initStep
    1     -2.8    -99.8     14.2      0.12  0.00146   0.0722    0.0722    absorber eIoni
    2     -2.8     -100     14.2      1.47 0.000137    0.135     0.207       World Transportation
    3    -3.07     -148     11.7       479  0.00318     47.9      48.1       World eIoni
    4    -3.87     -200     9.02     1e+03  0.00868     52.3       100  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 16,   Parent ID = 12
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    -3.07     -148     11.7   0.00106        0        0         0       World initStep
    1    -3.08     -148     11.7     0.236  0.00106    0.062     0.062       World eIoni
    2     -3.1     -200     11.7       523   0.0102       55      55.1  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 11,   Parent ID = 3
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    -2.73    -95.2     14.2  0.000996        0        0         0    absorber initStep
    1    -2.73    -95.2     14.2    0.0699 0.000996   0.0431    0.0431    absorber eIoni
    2    -2.97    -98.4     15.7      32.1  0.00342      4.4      4.45    absorber msc
    3    -3.12     -100     16.7      47.6  0.00024     1.86      6.31       World Transportation
    4    -4.07     -200     23.6  1.05e+03   0.0201      101       107  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 10,   Parent ID = 3
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    -2.66    -90.7     14.2   0.00404        0        0         0    absorber initStep
    1    -2.66    -90.7     14.2      0.52  0.00404    0.326     0.326    absorber eIoni
    2    -2.65    -91.6     14.2       9.1 0.000546     0.86      1.19    absorber eIoni
    3    -2.77     -100     14.2      92.9  0.00241     8.38      9.57       World Transportation
    4    -2.88     -140     14.6       490  0.00734     39.8      49.3       World eIoni
    5    -2.92     -200     14.8  1.09e+03   0.0108     60.3       110  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 17,   Parent ID = 10
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    -2.88     -140     14.6   0.00775        0        0         0       World initStep
    1     -2.9     -140     14.6      7.42  0.00513     1.46      1.46       World eIoni
    2    -2.48     -141     15.1      11.3 2.74e-05    0.777      2.24       World eIoni
    3   -0.421     -200     17.5       603   0.0084     59.5      61.7  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 18,   Parent ID = 17
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    -2.48     -141     15.1    0.0014        0        0         0       World initStep
    1    -2.49     -141     15.1     0.374   0.0014   0.0931    0.0931       World eIoni
    2    -2.51     -200     15.1       592   0.0115     60.8      60.8  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 2,   Parent ID = 1
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0   0.0991    -15.3     10.5     0.144        0        0         0    absorber initStep
    1     -3.2    -32.2     2.08       133   0.0199     33.5      33.5    absorber eBrem
    2    -13.1    -47.1      -20       282    0.013     32.5      66.1       World Transportation
    3    -17.2    -85.2    -29.3       663  0.00727     39.6       106       World eIoni
    4    -21.9     -199    -39.8   1.8e+03   0.0197      114       220  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = e-,   Track ID = 20,   Parent ID = 2
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0    -17.2    -85.2    -29.3   0.00108        0        0         0       World initStep
    1    -17.2    -85.2    -29.3     0.243  0.00108   0.0635    0.0635       World eIoni
    2    -17.2     -199    -29.3  1.14e+03   0.0267      128       128  OutOfWorld Transportation

*********************************************************************************************************
* G4Track Information:   Particle = gamma,   Track ID = 19,   Parent ID = 2
*********************************************************************************************************

Step#    X(mm)    Y(mm)    Z(mm) KinE(MeV)  dE(MeV) StepLeng TrackLeng  NextVolume ProcName
    0     -3.2    -32.2     2.08      37.1        0        0         0    absorber initStep
    1    -12.8    -18.2      -20      37.1        0     27.9      27.9       World Transportation
    2     -130      152     -288      37.1        0      338       366  OutOfWorld Transportation

````

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*

