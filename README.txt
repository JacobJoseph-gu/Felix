!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! Felix
!
! Richard Beanland, Keith Evans & Rudolf A Roemer
!
! (C) 2013-17, all rights reserved
!
! Version: :VERSION:
! Date:    :DATE:
! Time:    :TIME:
! Status:  :RLSTATUS:
! Build:   :BUILD:
! Author:  :AUTHOR:
! 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Felix is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  Felix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with Felix.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  Contributors: Alex Hubert, Jacob Richardson   
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

README.txt file for felixsim bloch wave method diffraction pattern simulation software.

Input File:

IWriteFlag: Controls the amount of information printed to the screen or
logfile during simulation, values 0-11 increase amount with 0 resulting in no
information being printed appart from start and stop messages.  More output
allows for easier debugging but will slow the execution.

IImageFLAG: Determines what images will be produced by the software :

0 Montages of Diffraction patterns  one per thickness in .bin file format
1 Individual reflections saved in .bin file format, Reflections
will be bundled into seperate folders for each thickness.
2 Amplitude and phase images will be saved, phase and
amplitude images will be saved individually for each reflection and labelled
-A or -P respectively

Any combination (in any order) of up to three of these options can be specified to suit the users requirements

IOutputFLAG: Determines the amount of calculated Variables which are saved
for later use by felixdraw.

0 Nothing is saved (Fastest)
1 Ug Matrix will be saved in binary format
2 Eigenspectra is saved in binary format
3 Wavefunctions are saved in binary format

Any combination (in any order) of up to three of these options can be specified to suit the users requirements

IBinorTextFLAG: Select binary or text output files, binary are smaller and faster to read/write however text files are easier to import into other programs for later use

0 Binary
1 Text

IScatteringFactorMethod: Determines which method by which to calculate potentials
0 Kirkands Method (103 Elements)
1 Peng (98 Elements)
2 Doyle and Turner (only works for 68 elements)

IZolzFLAG: Choose to limit the simulation to the zeroth order laue zone.  

0 No (Includes HOLZ, slower) 
1 Yes (ZOLZ only, faster)

ICentralBeamFlag: Exclude the [000] beam from the final images which improves
relative intensity

0 No central beam
1 Central beam included

IMaskFLAG: Chooses between a circular or square input beam

0 Circular
1 Square

IAbsorbFLAG: Choose to include absorption in simulation

0 no
1 Proportional Model
2 Einstein Model (Perturbative approach)
3 Einstein Model (Exact Approach)

IAnisotropicDebyeWallerFLAG: Choose to use anisotropic debye waller factors if
available

0 no
1 yes

IBeamConvergence: Not Yet Implemented

IPseudoCubicFLAG: Indicate whether given directions are expressed in PseudoCubic
or Orthorhombic notation

0 Orthorhomnic
1 PseudoCubic

IXDirectionFLAG : If set to 0 the program will ignore any specified X Direction and take the shortest G Vector as the X direction

0 Ignore X Direction
1 Use X Direction

IPixelCount: Pixel Radius of images, simulation scales as the square of this
  number but primary parallelisation is over pixels (more pixels, more cores
  can be used effectively) 64 is a good for images, 128+ better for
  quantitative analysis

IMinReflectionPool : Controls the size of the reflection pool accessible 
   during diagonalisation

IMinStrongBeams : This paramater sets the minimum number of beams overwhich the
  diagonisation is preformed, increasing it will result in greater accuracy
  but will slow simultion

IMinWeakBeams: Minimum number of weak beams with which to perturb the Strong beams

RBSBmax : Maximum weak beam purturbation strength before the beam is 
  considered strong 

RBSPmax : Maximum weak beam purturbation of a prior weak beam

RConvergenceTolerance (%) : Not yet implemented 
	       
RDebyeWallerConstant: If no Debye waller factor is found in the .cif
  file this value will be used for all atomic species missing the factor, note
  this is the B factor not U

RAbsorptionPer: Defines the percentage of absorption applied to the
  potentials when using the propotional model, this values is used for all atomic species

RConvergenceAngle: Defines the convergence angle of the beam in units of half
  the minimum gvector magnitude, at a value of 1 all beams will touch at
  tjtheir edges, values greater than 1 would (experimentally) cause beams to
  overlap, in the simulation this will not occur

IIncidentBeamDirectionX: X Component of the incident beam direction (Zone
  axis) expressed in the crystal reference frame in real space

IIncidentBeamDirectionY: Y Component of the incident beam direction (Zone
  axis) expressed in the crystal reference frame in real space

IIncidentBeamDirectionZ: Z Component of the incident beam direction (Zone
  axis) expressed in the crystal reference frame in real space

IXDirectionX: X component of the chosen X-axis expressed in the crystal
  reference frame in reciprocal space

IXDirectionY: Y component of the chosen X-axis expressed in the crystal
  reference frame in reciprocal space

IXDirectionZ: Z component of the chosen X-axis expressed in the crystal
  reference frame in reciprocal space

INormalDirectionX: X component of the plane normal to the surface of crystal
  in real space

INormalDirectionY: Y component of the plane normal to the surface of crystal
  in real space

INormalDirectionZ: Z component of the plane normal to the surface of crystal
  in real space

RAccelerationVoltage: Acceleration voltage of the microscope expressed in KV

RInitialThickness: Lower bound thickness to be applied (Angstroms)

RFinalThickness: Upper Bound Thickness to be Appliced (Angstroms)

RDeltaThickness: Step between thickness (Angstroms)

IReflectOut: The number of the reflections to be included in the final
  image(s)

The Remainder of the options apply to FelixRefine only

IImageOutputFLAG: Choose whether FelixRefine will output images at the conclusion of the refinement

0 No
1 Yes

IDevFLAG: Unused awaiting Removal

IRefineModeFLAG: Choose the refinement variable(s)

0 Refine Debye-waller Factor
1 Refine Structure Factors (UGs)
2 Refine Thickness

This options 
