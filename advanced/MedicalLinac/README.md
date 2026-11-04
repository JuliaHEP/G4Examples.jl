# Example MedicalLinac

Port of the Geant4 advanced example
[medical_linac](https://github.com/Geant4/geant4/tree/master/examples/advanced/medical_linac)
(S. Pozzi, B. Caccia and C. Mancini Terracciano). It simulates the head of a GE Saturn 43 medical
linear accelerator in the 12 MV photon mode with a 10 × 10 cm² field, and computes the dose in a
water phantom. The geometry and the measured dose distributions of this accelerator were published
by the EURADOS Working Group 6 (EURADOS Report 2020-05).

## Detector

The geometry is defined in `MedicalLinacDetector.jl`. The field side, the source to skin distance
(SSD) and the phantom side are keyword parameters of the `MedicalLinacDetector` structure. The
dimensions of the head are fixed.

- **Target**: a W block with a hole on the beam side, and a thin Ti disk.
- **Primary collimator**: three tubes (carbon steel, W-Ni-Cu alloy and Pb) with conical holes.
- **Flattening filter**: a steel tube and three cones.
- **Ionization chamber**: six Al layers, each followed by a kapton layer.
- **Vacuum window**: an Al disk.
- **Jaws**: four stacks of steel, Pb and W-Ni-Cu blocks inside air envelopes. The envelopes are
  rotated around the source so that the field side is 10 cm at SSD + 10 cm.
- **Phantom**: a water cube of 25.5 cm side at SSD = 90 cm, behind a 2.5 mm water layer.

Each element of the head is in its own `G4Region` with the production cuts of the C++ example.
These are 1 mm in the target, the flattening filter and the window, 1 cm in the primary
collimator, 0.1 mm in the ionization chamber, 2 cm in the jaws and 0.1 mm in the phantom.

## Primary particles

Electrons start 5 mm before the source, inside the hole of the target:

- **Energy:** Gaussian, 12 MeV mean with a standard deviation of 0.127 MeV.
- **Start point:** uniform in radius within 1 mm.
- **Direction:** Gaussian sin θ with a standard deviation of 0.003.

## Physics and variance reduction

The physics list is the option 4 of the standard EM physics with the decay physics. It is built
from QBBC with its hadronic physics removed.

As in the C++ example, the following particles are killed:

- particles that go backwards
- particles that enter a jaw envelope
- particles that leave the target, the upper and middle collimator or the flattening filter
  outside the cone of the flattening filter

In addition, the directional bremsstrahlung splitting of Geant4 (`G4EmParameters`) is used. Each
bremsstrahlung photon is split 1000 times, and only the copies that go towards a sphere of 15 cm
radius around the centre of the phantom surface are kept. The others survive with a probability
of 1/1000 and a weight of 1. The results are not biased.

The splitting is applied in all the regions of the head, and also to Compton scattering, the
photoelectric effect, pair production and annihilation. If it is applied only to the
bremsstrahlung in the target, a weight 1 photon sometimes scatters towards the phantom. A single
such photon deposits more energy in a voxel than all the split photons together, and shows up
as a spike in the dose.

## Output

The dose is scored with a `G4JLScoringMesh` of 51 × 51 × 51 voxels of 5 mm that covers the
phantom, the same mesh as in `run.mac` of the C++ example. The script shows:

- the geometry of the head and of the whole setup
- the percentage depth dose on the beam axis
- the lateral profiles at the depth of the maximum, 10 cm and 20 cm
- the dose map in the central plane

## Differences with the C++ example

- **EM physics:** option 4 instead of option 3, because `G4EmStandardPhysics_option3` is not
  available in Geant4.jl. Option 4 is also an option of the C++ example
  (`/physics/addPhysics standard_opt4`).
- **Primary collimator solids:** each piece is a `G4Polycone` instead of a tube minus a hole
  (`G4UnionSolid` and `G4SubtractionSolid`). The solids are the same, and the polycones can also be
  drawn cut in half.
- **Directional bremsstrahlung splitting** is added. With the same CPU time, the statistical errors
  are about 4 times smaller.
- **Phantom voxels:** the phantom is not divided into replicas. The dose is scored by the scoring
  mesh, as in the C++ example, so the replicas are not needed.
- **Messengers:** the UI commands that change the geometry, the beam and the physics are replaced by
  the parameters of the detector and the generator.
- **Kept from the C++ example:**
  - The envelopes of the X and Y jaws overlap in their air parts. Their overlaps are not checked.
  - The lower part of the primary collimator is not used in the cone cut,
    because of a typo in its volume name there.

## Validation

The C++ example was built against the same Geant4 version (11.4.1) and run with option 4 and
200000 electrons. The total dose in the phantom agrees with this port within 5% without splitting
and within 3% with splitting, consistent with the statistical uncertainty.

## Running

```sh
julia --project=. advanced/MedicalLinac/MedicalLinac.jl
```

The script uses 8 threads. 100000 electrons take about one minute. Increase `nevents` for
smoother curves.
