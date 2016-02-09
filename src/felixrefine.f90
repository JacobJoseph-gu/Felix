!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! felixrefine
!
! Richard Beanland, Keith Evans, Rudolf A Roemer and Alexander Hubert
!
! (C) 2013/14, all right reserved
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
!  This file is part of felixrefine.
!
!  felixrefine is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  felixrefine is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with felixrefine.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! $Id: Felixrefine.f90,v 1.89 2014/04/28 12:26:19 phslaz Exp $
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PROGRAM Felixrefine
 
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  !--------------------------------------------------------------------
  ! local variable definitions
  IMPLICIT NONE

  INTEGER(IKIND) :: IHours,IMinutes,ISeconds,IErr,IMilliSeconds,IIterationFLAG,&
       ind,jnd,ICalls,IIterationCount
  REAL(RKIND) :: StartTime, CurrentTime, Duration, TotalDurationEstimate,&
       RFigureOfMerit,SimplexFunction,RHOLZAcceptanceAngle
  INTEGER(IKIND) :: IStartTime, ICurrentTime ,IRate
  REAL(RKIND),DIMENSION(:,:),ALLOCATABLE :: RSimplexVolume
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RSimplexFoM,RIndependentVariableValues
  REAL(RKIND) :: RBCASTREAL,RStandardDeviation,RMean
  CHARACTER*40 my_rank_string ,SPrintString

  !-------------------------------------------------------------------
  ! constants
  CALL Init_Numbers
  
  !-------------------------------------------------------------------
  ! set the error value to zero, will change upon error
  IErr=0

  !--------------------------------------------------------------------
  ! MPI initialization
  CALL MPI_Init(IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Init()"
     GOTO 9999
  END IF

  ! Get the rank of the current process
  CALL MPI_Comm_rank(MPI_COMM_WORLD,my_rank,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Comm_rank()"
     GOTO 9999
  END IF

  ! Get the size of the current communicator
  CALL MPI_Comm_size(MPI_COMM_WORLD,p,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error in MPI_Comm_size()"
     GOTO 9999
  END IF

  !--------------------------------------------------------------------
  ! protocal feature startup
  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"--------------------------------------------------------------"
     PRINT*,"Felixrefine: ", RStr
     PRINT*,"          ", DStr
     PRINT*,"          ", AStr
     PRINT*,"          on rank= ", my_rank, " of ", p, " in total."
     PRINT*,"--------------------------------------------------------------"
  END IF
  ISoftwareMode =2 ! felixrefinemode

  !--------------------------------------------------------------------
  ! timing startup
  !--------------------------------------------------------------------

  CALL SYSTEM_CLOCK(count_rate=IRate)
  CALL SYSTEM_CLOCK(IStarttime)

  !--------------------------------------------------------------------
  ! INPUT section 
  !--------------------------------------------------------------------
  CALL ReadInput (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error in ReadInput"
     GOTO 9999
  END IF  
  
  ALLOCATE(RImageExpi(2*IPixelCount,2*IPixelCount,INoOfLacbedPatterns),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error in allocation of RImageExpi"
     GOTO 9999
  END IF

  CALL ReadExperimentalImages(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error in ReadExperimentalImages"
     GOTO 9999
  END IF

  !--------------------------------------------------------------------
  ! Initial simulation and variable setup
  !--------------------------------------------------------------------
  CALL MicroscopySettings( IErr )
  IF( IErr.NE.0 ) THEN
    PRINT*,"felixrefine(",my_rank,") error in MicroscopySettings"
    GOTO 9999
  ENDIF  
  
  CALL CrystalLatticeVectorDetermination(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,") error in CrystalLatticeVectorDetermination"
     GOTO 9999
  ENDIF

  ALLOCATE(RrVecMat(ITotalAtoms,THREEDIM),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error allocating RrVecMat"
     GOTO 9999
  ENDIF
  
  CALL AllAtomPositions(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in AllAtomPositions"
     GOTO 9999
  ENDIF

!zz temp deallocation to get it to work
DEALLOCATE(RFullPartialOccupancy,SMNP,MNP,RFullAtomicFracCoordVec,SFullAtomicNameVec, &
RFullIsotropicDebyeWallerFactor,IFullAtomNumber,IFullAnisotropicDWFTensor,&
RDWF,ROcc,IAtoms,IAnisoDWFT,RrVecMat)

!zz from diffractionpatterninitialisation/reflectiondetermination
  ind = 0!here acts as a flag
  jnd = 0
  IhklMaxValue = 15!RB starting value for maximum hkl, increments if necessary
  RHOLZAcceptanceAngle=TWODEG2RADIAN!RB maximum acceptance angle for HOLZ, suspect way too low
  DO WHILE (ind.EQ.0)
     jnd = jnd+1     
     CALL NewHKLMake(IhklMaxValue,RZDirC,RHOLZAcceptanceAngle,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine(",my_rank,")error in NewHKLMake()"
        GOTO 9999
     END IF
     IF(SIZE(Rhkl,DIM=1).LT.IMinReflectionPool) THEN
        IhklMaxValue = IhklMaxValue*2
        DEALLOCATE(Rhkl,STAT=ierr)!zz not happy about conditional deallocation with no corresponding allocation
        IF( IErr.NE.0 ) THEN
           PRINT*,"felixrefine(",my_rank,")error deallocating Rhkl"
           GOTO 9999
        END IF
        CYCLE
     ELSE
        ind = 1
     END IF
  END DO
  
  CALL SortHKL(Rhkl,SIZE(Rhkl,1),IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine(",my_rank,")error in SortHKL"
     GOTO 9999
  END IF

!zz temp deallocation to get it to work
DEALLOCATE(Rhkl)
  
  !--------------------------------------------------------------------
  ! Setup Simplex Variables
  !--------------------------------------------------------------------
  IF(IRefineModeSelectionArray(2).EQ.1) THEN 
     CALL SetupAtomicVectorMovements(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine (", my_rank, ") error in SetupAtomicVectorMovements"
        GOTO 9999
     END IF
  END IF

  !--------------------------------------------------------------------
  !  DetermineNumberofRefinementVariablesPerType
  INoofElementsForEachRefinementType(1) = &
       IRefineModeSelectionArray(1)*(INoofUgs*2+1)!RB +1 is for absorption
  INoofElementsForEachRefinementType(2) = &
       IRefineModeSelectionArray(2)*IAllowedVectors
  INoofElementsForEachRefinementType(3) = &
       IRefineModeSelectionArray(3)*SIZE(IAtomicSitesToRefine)
  INoofElementsForEachRefinementType(4) = &
       IRefineModeSelectionArray(4)*SIZE(IAtomicSitesToRefine)
  INoofElementsForEachRefinementType(5) = &
       IRefineModeSelectionArray(5)*SIZE(IAtomicSitesToRefine)*6
  INoofElementsForEachRefinementType(6) = &
       IRefineModeSelectionArray(6)*3
  INoofElementsForEachRefinementType(7) = &
       IRefineModeSelectionArray(7)*3
  INoofElementsForEachRefinementType(8) = &
       IRefineModeSelectionArray(8)
  INoofElementsForEachRefinementType(9) = &
       IRefineModeSelectionArray(9)
  INoofElementsForEachRefinementType(10) = &
       IRefineModeSelectionArray(10)
  INoofElementsForEachRefinementType(11) = &
       IRefineModeSelectionArray(11)
  
  IIndependentVariables = SUM(INoofElementsForEachRefinementType)
  IF(my_rank.EQ.0) THEN
    IF ( IIndependentVariables.EQ.1 ) THEN 
      PRINT*,"Only one independent variable"
	ELSE
      WRITE(SPrintString,FMT='(I3,1X,A21))') IIndependentVariables,"independent variables"
      PRINT*,TRIM(ADJUSTL(SPrintString))
    END IF
  END IF

  !allocations--------------------------------------------------------------------
  ALLOCATE(IIterativeVariableUniqueIDs(IIndependentVariables,5),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error allocating IIterativeVariableUniqueIDs"
     GOTO 9999
  ENDIF
  ALLOCATE(RIndependentVariableValues(IIndependentVariables),STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error allocating RIndependentVariableValues"
     GOTO 9999
  END IF

  !--------------------------------------------------------------------
  !  Assign IDs  
  IIterativeVariableUniqueIDs = 0
  ICalls = 0

  DO ind = 1,IRefinementVariableTypes !Loop over all possible iterative variables
     IF(IRefineModeSelectionArray(ind).EQ.1) THEN
        DO jnd = 1,INoofElementsForEachRefinementType(ind)
           ICalls = ICalls + 1
           IIterativeVariableUniqueIDs(ICalls,1) = ICalls
           CALL AssignArrayLocationsToIterationVariables(ind,jnd,IIterativeVariableUniqueIDs,IErr)
        END DO
     END IF
  END DO 
  DO ind=1,IIndependentVariables!RB debug
    PRINT*, IIterativeVariableUniqueIDs(ind,:)
  END DO
    
  !--------------------------------------------------------------------
  ! Initialise Simplex
  !--------------------------------------------------------------------

  ALLOCATE(RSimplexVolume(IIndependentVariables+1,IIndependentVariables),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation()"
     GOTO 9999
  END IF

  ALLOCATE(RSimplexFoM(IIndependentVariables),&
       STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation()"
     GOTO 9999
  END IF
  
  IIterationCount = 0

  CALL SimplexInitialisation(RSimplexVolume,RSimplexFoM,RIndependentVariableValues,IIterationCount,RStandardDeviation,RMean,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in SimplexInitialisation()"
     GOTO 9999
  END IF
     
  !--------------------------------------------------------------------
  ! Apply Simplex Method
  !--------------------------------------------------------------------

  CALL NDimensionalDownhillSimplex(RSimplexVolume,RSimplexFoM,&
       IIndependentVariables+1,&
       IIndependentVariables,IIndependentVariables,&
       RExitCriteria,IIterationCount,RStandardDeviation,RMean,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in NDimensionalDownhillSimplex()"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! Deallocate Memory
  !--------------------------------------------------------------------
  DEALLOCATE(IIterativeVariableUniqueIDs,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error deallocating IIterativeVariableUniqueIDs"
     GOTO 9999
  ENDIF

  DEALLOCATE(RImageExpi,STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error deallocating RImageExpi"
     GOTO 9999
  ENDIF

  !--------------------------------------------------------------------
  ! finish off
  !--------------------------------------------------------------------
  
  WRITE(my_rank_string,*) my_rank
    
  CALL SYSTEM_CLOCK(ICurrentTime)
  Duration=REAL(ICurrentTime-IStartTime)/REAL(IRate)
  IHours = FLOOR(Duration/3600.0D0)
  IMinutes = FLOOR(MOD(Duration,3600.0D0)/60.0D0)
  ISeconds = MOD(Duration,3600.0D0)-IMinutes*60
  IMilliSeconds = INT((Duration-(IHours*3600+IMinutes*60+ISeconds))*1000,IKIND)

  PRINT*, "felixrefine( ", TRIM(ADJUSTL(my_rank_string)), " ) ", &
       RStr, ", used time=", IHours, "hrs ", &
       IMinutes,"mins ",ISeconds,"secs ", IMilliSeconds,"millisecs"

  !--------------------------------------------------------------------
  ! Shut down MPI
  !--------------------------------------------------------------------

9999 &
  CALL MPI_Finalize(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"Felixrefine(", my_rank, ") error ", IErr, " in MPI_Finalize()"
     STOP
  ENDIF
  
  ! clean shutdown
  STOP
  
END PROGRAM Felixrefine
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE AssignArrayLocationsToIterationVariables(IIterativeVariableType,IVariableNo,IArrayToFill,IErr)
!NB IArrayToFill here is equivalent to IIterativeVariableUniqueIDs outside this subroutine
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IIterativeVariableType,IVariableNo,IErr,IArrayIndex,&
       IAnisotropicDebyeWallerFactorElementNo
  INTEGER(IKIND),DIMENSION(IIndependentVariables,5),INTENT(OUT) :: IArrayToFill  

!!$  Calculate How Many of Each Variable Type There are
!  CALL DetermineNumberofRefinementVariablesPerType(INoofElementsForEachRefinementType,IErr)
  
!!$  Where am I in the Array Right Now?
  IArrayIndex = SUM(INoofElementsForEachRefinementType(:(IIterativeVariableType-1)))+IVariableNo

  SELECT CASE(IIterativeVariableType)

  CASE(1) ! Ugs
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = &
          NINT(REAL(INoofUgs,RKIND)*(REAL(IVariableNo/REAL(INoofUgs,RKIND),RKIND)-&
          CEILING(REAL(IVariableNo/REAL(INoofUgs,RKIND),RKIND)))+REAL(INoofUgs,RKIND))

  CASE(2) ! Coordinates (x,y,z)
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IVariableNo

  CASE(3) ! Atomic Site Occupancies
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(IVariableNo)

  CASE(4) ! Isotropic Debye Waller Factors 
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(IVariableNo)

  CASE(5) ! Anisotropic Debye Waller Factors (a11-a33)
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = IAtomicSitesToRefine(INT(CEILING(REAL(IVariableNo/6.0D0,RKIND))))
     IAnisotropicDebyeWallerFactorElementNo = &
          NINT(6.D0*(REAL(IVariableNo/6.0D0,RKIND)-CEILING(REAL(IVariableNo/6.0D0,RKIND)))+6.0D0)

     SELECT CASE(IAnisotropicDebyeWallerFactorElementNo)

        CASE(1)
           IArrayToFill(IArrayIndex,4:5) = [1,1]
        CASE(2)
           IArrayToFill(IArrayIndex,4:5) = [2,1]
        CASE(3)
           IArrayToFill(IArrayIndex,4:5) = [2,2]
        CASE(4)
           IArrayToFill(IArrayIndex,4:5) = [3,1]
        CASE(5)
           IArrayToFill(IArrayIndex,4:5) = [3,2]
        CASE(6)
           IArrayToFill(IArrayIndex,4:5) = [3,3]

        END SELECT

  CASE(6) ! Lattice Parameters
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = &
          NINT(3.D0*(REAL(IVariableNo/3.0D0,RKIND)-CEILING(REAL(IVariableNo/3.0D0,RKIND)))+3.0D0)
     
  CASE(7) ! Lattice Angles
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     IArrayToFill(IArrayIndex,3) = &
          NINT(3.D0*(REAL(IVariableNo/3.0D0,RKIND)-CEILING(REAL(IVariableNo/3.0D0,RKIND)))+3.0D0)

  CASE(8) 
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType

  CASE(9)  
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType

  CASE(10)
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType

  CASE(11)
     IArrayToFill(IArrayIndex,2) = IIterativeVariableType
     
  END SELECT
  
END SUBROUTINE AssignArrayLocationsToIterationVariables

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RefinementVariableSetup(RIndependentVariableValues,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind,IVariableType
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(OUT) :: RIndependentVariableValues
  
  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RefinementVariableSetup(",my_rank,")"
  END IF
  
!!$  Fill the Independent Value array with values

  DO ind = 1,IIndependentVariables
     IVariableType = IIterativeVariableUniqueIDs(ind,2)
     SELECT CASE (IVariableType)
     CASE(1)
	    !Structure factor refinement, define in SymmetryRelatedStructureFactorDetermination
     CASE(2)
        RIndependentVariableValues(ind) = &
             RAllowedVectorMagnitudes(IIterativeVariableUniqueIDs(ind,3))
     CASE(3)
        RIndependentVariableValues(ind) = &
             RAtomicSitePartialOccupancy(IIterativeVariableUniqueIDs(ind,3))
     CASE(4)
        RIndependentVariableValues(ind) = &
             RIsotropicDebyeWallerFactors(IIterativeVariableUniqueIDs(ind,3))
     CASE(5)
        RIndependentVariableValues(ind) = &
             RAnisotropicDebyeWallerFactorTensor(&
             IIterativeVariableUniqueIDs(ind,3),&
             IIterativeVariableUniqueIDs(ind,4),&
             IIterativeVariableUniqueIDs(ind,5))
     CASE(6)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RIndependentVariableValues(ind) = RLengthX
        CASE(2)
           RIndependentVariableValues(ind) = RLengthY
        CASE(3)
           RIndependentVariableValues(ind) = RLengthZ
        END SELECT
     CASE(7)
        SELECT CASE(IIterativeVariableUniqueIDs(ind,3))
        CASE(1)
           RIndependentVariableValues(ind) = RAlpha
        CASE(2)
           RIndependentVariableValues(ind) = RBeta
        CASE(3)
           RIndependentVariableValues(ind) = RGamma
        END SELECT
     CASE(8)
        RIndependentVariableValues(ind) = &
             RConvergenceAngle
     CASE(9)
        RIndependentVariableValues(ind) = &
             RAbsorptionPercentage
     CASE(10)
        RIndependentVariableValues(ind) = &
             RAcceleratingVoltage
     CASE(11)
        RIndependentVariableValues(ind) = &
             RRSoSScalingFactor
     END SELECT
  END DO

END SUBROUTINE RefinementVariableSetup
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE StructureFactorRefinementSetup(RIndependentVariableValues,IIterationCount,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(OUT) :: RIndependentVariableValues
  INTEGER(IKIND),INTENT(IN) :: IIterationCount
  CHARACTER*200 :: SPrintString

  IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"StructureFactorRefinementSetup(",my_rank,")"
  END IF

  IF(IRefineModeSelectionArray(1).EQ.1) THEN
     DO ind = 1,INoofUgs !RB ignore the first one as it is the internal potential
        RIndependentVariableValues((ind-1)*2+1) = &!yy
             REAL(CSymmetryStrengthKey(ind+1),RKIND)!yy ind+1 instead of ind
        RIndependentVariableValues((ind-1)*2+2) = &
             AIMAG(CSymmetryStrengthKey(ind+1))!yy ind+1 instead of ind
     END DO
  END IF
  RIndependentVariableValues(2*INoofUgs+1) = RAbsorptionPercentage!RB absorption always included in structure factor refinement as last variable

!yyDO ind = 1,2*INoofUgs+1!yy
!yy  WRITE(SPrintString,FMT='(A3,I2,A1,F7.4)') "RB ",ind,":",RIndependentVariableValues(ind)
!yy  PRINT*,TRIM(ADJUSTL(SPrintString))
!yy END DO


END SUBROUTINE StructureFactorRefinementSetup

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RankSymmetryRelatedStructureFactor(IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE 
  
  INTEGER(IKIND) :: IErr,ind
  INTEGER(IKIND),DIMENSION(2) :: ILoc

  IF((IWriteFLAG.GE.10.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
     PRINT*,"RankSymmetryRelatedStructureFactor(",my_rank,")"
  END IF
  
  ALLOCATE(ISymmetryRelations(nReflections,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RankSymmetryRelatedStructureFactor(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of ISymmetryRelations"
     RETURN
  ENDIF
  
  CALL SymmetryRelatedStructureFactorDetermination (IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"RankSymmetryRelatedStructureFactor(", my_rank, ") error ", IErr, &
          " in SymmetryRelatedStructureFactorDetermination"
     RETURN
  ENDIF
  
  DO ind = 1,(SIZE(ISymmetryStrengthKey))
     ILoc = MINLOC(ABS(ISymmetryRelations-ind))
     ISymmetryStrengthKey(ind) = ind
     CSymmetryStrengthKey(ind) = CUgMatNoAbs(ILoc(1),ILoc(2))
!RB     PRINT*,"CSymmetryStrengthKey",ind,ISymmetryStrengthKey(ind),CSymmetryStrengthKey(ind)
  END DO
  
  CALL ReSortUgs(ISymmetryStrengthKey,CSymmetryStrengthKey,SIZE(CSymmetryStrengthKey,DIM=1))

END SUBROUTINE RankSymmetryRelatedStructureFactor

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SimplexInitialisation(RSimplexVolume,RSimplexFoM,RIndependentVariableValues,&
     IIterationCount,RStandardDeviation,RMean,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind,jnd,IExitFLAG
  LOGICAL :: LInitialSimulationFLAG = .TRUE. ! Its value is meaningless :)
  REAL(RKIND),DIMENSION(IIndependentVariables+1,IIndependentVariables),INTENT(OUT) :: RSimplexVolume
  REAL(RKIND),DIMENSION(IIndependentVariables+1),INTENT(OUT) :: RSimplexFoM
  REAL(RKIND) :: SimplexFunction,RSimplexDummy
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(INOUT) :: RIndependentVariableValues
  INTEGER(IKIND),INTENT(INOUT) :: IIterationCount
  REAL(RKIND),INTENT(OUT) :: RStandardDeviation,RMean
  REAL(RKIND) :: RStandardError,RStandardTolerance
  CHARACTER*200 :: SPrintString

  IF(IWriteFLAG.GE.10.AND.my_rank.EQ.0) THEN
     PRINT*,"SimplexInitialisation(",my_rank,")"
  END IF
      
  CALL FelixFunction(LInitialSimulationFLAG,IErr)!RB first thing!!
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(", my_rank, ") error in PerformInitialSimulation()"
     RETURN
  ENDIF
       
  CALL InitialiseWeightingCoefficients(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(", my_rank, ") error in InitialiseWeightingCoefficients()"
     RETURN
  ENDIF

  IF(my_rank.EQ.0) THEN   
     IThicknessCount= (RFinalThickness- RInitialThickness)/RDeltaThickness + 1
     IIterationCount = 0; !Initial Simulation is iteration zero
     IExitFLAG = 0; ! Do not exit
     IPreviousPrintedIteration = -IPrint
     CALL CreateImagesAndWriteOutput(IIterationCount,IExitFLAG,IErr) 
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexFunction(", my_rank, ") error ", IErr, &
             " in CreateImagesAndWriteOutput"
        RETURN
     ENDIF
  ELSE
     DEALLOCATE(Rhkl,STAT=IErr)  
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexInitialisation (", my_rank, ") error in Deallocation()"
        RETURN
     ENDIF
  END IF
  
  CALL RefinementVariableSetup(RIndependentVariableValues,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation(", my_rank, ") error in RefinementVariableSetup()"
     RETURN
  ENDIF
  
  IF(IRefineModeSelectionArray(1).EQ.1) THEN
     
     CALL RankSymmetryRelatedStructureFactor(IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexInitialisation(", my_rank, ") error in RankSymmetryRelatedStructureFactor()"
        RETURN
     ENDIF
     
     CALL StructureFactorRefinementSetup(RIndependentVariableValues,IIterationCount,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexInitialisation(", my_rank, ") error in StructureFactorRefinementSetup()"
        RETURN
     ENDIF
     
  ENDIF
!RB     PRINT*,"Deallocating CUgMat,CUgMatNoAbs,CUgMatPrime in felixrefine" NB Also deallocated in felixfunction!!!
  DEALLOCATE(RgSumMat,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixsim(", my_rank, ") error ", IErr, &
          " in Deallocation of RgSumMat"
     RETURN
  ENDIF

  DEALLOCATE(CUgMatNoAbs,STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation (", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF

  DEALLOCATE(CUgMatPrime,STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation (", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF
 
  DEALLOCATE(CUgMat,STAT=IErr)  
  IF( IErr.NE.0 ) THEN
     PRINT*,"SimplexInitialisation (", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF

!!$ RandomSequence

  IF(IContinueFLAG.EQ.0) THEN
     IF(my_rank.EQ.0) THEN
        CALL CreateRandomisedSimplex(RSimplexVolume,&
             RIndependentVariableValues,IErr)

        CALL MPI_BCAST(RSimplexVolume,(IIndependentVariables+1)*(IIndependentVariables),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)
     ELSE
        CALL MPI_BCAST(RSimplexVolume,(IIndependentVariables+1)*(IIndependentVariables),MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,IErr)

     END IF

     IPreviousPrintedIteration = -IPrint ! Ensures print out on first iteration

     DO ind = 1,(IIndependentVariables+1)
        
        IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
           PRINT*,"--------------------------------"
           WRITE(SPrintString,FMT='(A8,I2,A4,I3)') "Simplex ",ind," of ",IIndependentVariables+1
           PRINT*,TRIM(ADJUSTL(SPrintString))
 !          PRINT*,"-------- Simplex",ind,"of",IIndependentVariables+1
           PRINT*,"--------------------------------"
        END IF

        RSimplexDummy = SimplexFunction(RSimplexVolume(ind,:),1,0,IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"SimplexInitialisation(", my_rank, ") error in SimplexFunction()"
           RETURN
        ENDIF
        
        RStandardTolerance = RStandardError(RStandardDeviation,RMean,RSimplexDummy,IErr)
        
        RSimplexFoM(ind) =  RSimplexDummy
        
        IF((IWriteFLAG.GE.0.AND.my_rank.EQ.0).OR.IWriteFLAG.GE.10) THEN
 !         PRINT*,"--------------------------------"
          WRITE(SPrintString,FMT='(A16,F7.5))') "Figure of merit ",RSimplexFoM(ind)
          PRINT*,TRIM(ADJUSTL(SPrintString))
!          PRINT*,"-------- Figure of Merit" ,RSimplexFoM(ind)        
 !         PRINT*,"---------------------------------------------------------"
        END IF
     END DO
     
  ELSE
     
     CALL RecoverSavedSimplex(RSimplexVolume,RSimplexFoM,RStandardDeviation,RMean,IIterationCount,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"SimplexInitialisation (", my_rank, ") error in RecoverSavedSimplex()"
        RETURN
     ENDIF
     
  END IF
       
END SUBROUTINE SimplexInitialisation

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CreateRandomisedSimplex(RSimplexVolume,RIndependentVariableValues,IErr)

USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara 
  USE IChannels
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) :: IErr,ind
  REAL(RKIND),DIMENSION(:),ALLOCATABLE :: RRandomSigns,RRandomNumbers
  REAL(RKIND),DIMENSION(IIndependentVariables+1,IIndependentVariables),INTENT(OUT) :: RSimplexVolume
  REAL(RKIND),DIMENSION(IIndependentVariables),INTENT(INOUT) :: RIndependentVariableValues
  
  IF(IRefineModeSelectionArray(2).EQ.1) THEN
     DO ind = 1,(IIndependentVariables+1)
        ALLOCATE(RRandomSigns(IAllowedVectors),RRandomNumbers(IAllowedVectors),&
             STAT=IErr)       
        
!!$           Randomise Atomic Displacements
        CALL RandomSequence(RRandomNumbers,IAllowedVectors,ind,IErr)
        CALL RandomSequence(RRandomSigns,IAllowedVectors,2*ind,IErr)
        WHERE (RRandomSigns.LT.HALF)
           RRandomSigns = ONE
        ELSEWHERE
           RRandomSigns = -ONE
        END WHERE
        
        RSimplexVolume(ind,:IAllowedVectors) = &
             RRandomNumbers*RRandomSigns*RSimplexLengthScale
        
        DEALLOCATE(RRandomSigns,RRandomNumbers)
        
        ALLOCATE(RRandomSigns(IIndependentVariables-IAllowedVectors),&
             RRandomNumbers(IIndependentVariables-IAllowedVectors),&
             STAT=IErr)
        
!!$           Randomise Everything else
        
        CALL RandomSequence(RRandomNumbers,&
             IIndependentVariables-IAllowedVectors,ind,IErr)
        CALL RandomSequence(RRandomSigns,&
             IIndependentVariables-IAllowedVectors,2*ind,IErr)
        WHERE (RRandomSigns.LT.HALF)
           RRandomSigns = ONE
        ELSEWHERE
           RRandomSigns = -ONE
        END WHERE
        
        RSimplexVolume(ind,(IAllowedVectors+1):) = &
             RIndependentVariableValues((IAllowedVectors+1):)*&
             (1+(RRandomNumbers*RRandomSigns*RSimplexLengthScale))
        
        DEALLOCATE(RRandomSigns,RRandomNumbers)
        
     END DO
     
  ELSE
     ALLOCATE(RRandomSigns(IIndependentVariables),&
          RRandomNumbers(IIndependentVariables),&
          STAT=IErr)
     
     DO ind = 1,(IIndependentVariables+1)
        
        CALL RandomSequence(RRandomNumbers,IIndependentVariables,ind,IErr)
        CALL RandomSequence(RRandomSigns,IIndependentVariables,2*ind,IErr)
        WHERE (RRandomSigns.LT.HALF)
           RRandomSigns = ONE
        ELSEWHERE
           RRandomSigns = -ONE
        END WHERE
        
        
        RSimplexVolume(ind,:) = &
             RIndependentVariableValues(:)*&
             (1+(RRandomNumbers*RRandomSigns*RSimplexLengthScale))
     END DO
     
        DEALLOCATE(RRandomSigns,RRandomNumbers)
     
  END IF
  

END SUBROUTINE CreateRandomisedSimplex

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE InitialiseAtomicVectorMagnitudes(IVariableID,RCorrectedMovement,IErr)
  
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Creates pseudo random movements of atoms using allowed vectors
!!$  % to initialise the simplex, proposed movements which exit the unit
!!$  $ cell are corrected to bring the atom back in on the opposite side
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara
  
  USE IChannels
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,IVariableID
  REAL(RKIND) :: &
       RNegativeMovement,RPositiveMovement,RCorrectedMovement,RANDOMNUMBER
  RNegativeMovement = RSimplexLengthScale*(-1.0_RKIND)
  RPositiveMovement = RSimplexLengthScale

  IF(RANDOMNUMBER(IVariableID,IErr).LT.0.5_RKIND) THEN
     CALL OutofUnitCellCheck(IVariableID,RNegativeMovement,RCorrectedMovement,IErr)
  ELSE
     CALL OutofUnitCellCheck(IVariableID,RPositiveMovement,RCorrectedMovement,IErr)
  END IF

END SUBROUTINE InitialiseAtomicVectorMagnitudes

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RandomSequence(RRandomSequence,IRandomSequenceLength,ISeedModifier,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Sets up a pseudo random sequence and selects a number
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,Ivalues(1:8), k,IRandomSequenceLength,ISeedModifier
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       seed
  REAL(RKIND),DIMENSION(IRandomSequenceLength) :: &
       RRandomSequence
  
  CALL DATE_AND_TIME(VALUES=Ivalues)

  IValues = IValues*ISeedModifier
!!$  CALL SYSTEM_CLOCK(
  IF (IRandomFLAG.EQ.0) THEN
     CALL RANDOM_SEED(size=k)
     allocate(seed(1:k))
     seed(:) = Ivalues(8)
     CALL RANDOM_SEED(put=seed)
  ELSE
     CALL RANDOM_SEED(size=k)
     ALLOCATE(seed(1:k))
     seed(:) = IFixedSeed*ISeedModifier
     CALL RANDOM_SEED(put=seed)
  END IF
   
  DEALLOCATE(seed)

  CALL RANDOM_NUMBER(RRandomSequence)
  
!!$  RANDOMSEQUENCE = RRandomNumberSequence(IRequestedNumber)
  
END SUBROUTINE  RANDOMSEQUENCE

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

REAL(RKIND) FUNCTION RANDOMNUMBER(IRequestedNumber,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Sets up a pseudo random sequence and selects a number
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,values(1:8), k,IRequestedNumber
  INTEGER(IKIND), DIMENSION(:), ALLOCATABLE :: &
       seed
  REAL(RKIND),DIMENSION(IRequestedNumber) :: &
       RRandomNumberSequence
  
  CALL DATE_AND_TIME(values=values)
  
  IF (IRandomFLAG.EQ.0) THEN
     CALL RANDOM_SEED(size=k)
     allocate(seed(1:k))
     seed(:) = values(8)
     CALL RANDOM_SEED(put=seed)
  ELSE
     CALL RANDOM_SEED(size=k)
     ALLOCATE(seed(1:k))
     seed(:) = IFixedSeed*IRequestedNumber
     CALL RANDOM_SEED(put=seed)
  END IF
   
  CALL RANDOM_NUMBER(RRandomNumberSequence)
  
  RANDOMNUMBER = RRandomNumberSequence(IRequestedNumber)
  
END FUNCTION RANDOMNUMBER

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE OutofUnitCellCheck(IVariableID,RProposedMovement,RCorrectedMovement,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Checks that vector movement applied by the simplex initialisation
!!$  % does not move an atom out fo the unit cell, and if it does
!!$  % the atom is moved back into the unit cell on the opposite side
!!$  % as if the atom had moved from one unit cell into the neighbouring one
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: ind,IErr,IVariableID,IAtomID,IVectorID
  REAL(RKIND),DIMENSION(THREEDIM) :: RProposedAtomicCoordinate,RDummyMovement
  REAL(RKIND),INTENT(IN) :: RProposedMovement
  REAL(RKIND),INTENT(OUT) :: RCorrectedMovement
  
  IVectorID = IIterativeVariableUniqueIDs(IVariableID,3)
  
  IAtomID = IAllowedVectorIDs(IVectorID)
 
  RProposedAtomicCoordinate(:) = RAtomSiteFracCoordVec(IAtomID,:) + &
       RProposedMovement*RAllowedVectors(IVectorID,:)

  RDummyMovement = RProposedMovement

  IF(ANY(RProposedAtomicCoordinate.GT.ONE).OR.ANY(RProposedAtomicCoordinate.LT.ZERO)) THEN
     DO ind = 1,THREEDIM
        IF (RProposedAtomicCoordinate(ind).GT.ONE) THEN
           RDummyMovement(ind) = (ONE-RAtomSiteFracCoordVec(IAtomID,ind))/RAllowedVectors(IVectorID,ind)
        ELSEIF(RProposedAtomicCoordinate(ind).LT.ZERO) THEN
           RDummyMovement(ind) = (-RAtomSiteFracCoordVec(IAtomID,ind))/RAllowedVectors(IVectorID,ind)
        ELSE
           RDummyMovement(ind) = RProposedMovement
        END IF
     END DO
  END IF

  IF(RProposedMovement.LT.ZERO) THEN
     RCorrectedMovement = MAXVAL(RDummyMovement)
  ELSE
     RCorrectedMovement = MINVAL(RDummyMovement)
  END IF

END SUBROUTINE OutofUnitCellCheck

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE ApplyNewStructureFactors(IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % Subroutine to place iteratively determined Structure factors
!!$  % to Ug Matrix
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ind
  COMPLEX(CKIND),DIMENSION(nReflections,nReflections) :: CUgMatDummy

!!$  Dummy Matrix to contain new iterative values
  
   CUgMatDummy = CZERO

!!$  Populate Ug Matrix with new iterative elements, maintain Hermiticity
  DO ind = 1,INoofUgs
     WHERE(ISymmetryRelations.EQ.ISymmetryStrengthKey(ind))
        CUgMatDummy = CSymmetryStrengthKey(ind)
     END WHERE
     WHERE(ISymmetryRelations.EQ.-ISymmetryStrengthKey(ind))!RB 
        CUgMatDummy = CONJG(CSymmetryStrengthKey(ind))
     END WHERE
  END DO

  WHERE(ABS(CUgMatDummy).GT.TINY)
     CUgMatNoAbs = CUgMatDummy!RB
  END WHERE

!!$  CUgMatNoAbs now contains the new values from the iterative process
! N.B. CUgMatNoAbs is without absorption and therefore Hermitian
  
END SUBROUTINE ApplyNewStructureFactors

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE CreateIdentityMatrix(IIdentityMatrix,ISize,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % This Subroutine creates an identity matrix of size
!!$  % ISize * ISize
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER(IKIND) :: &
       IErr,ISize,ind
  INTEGER(IKIND),DIMENSION(ISize,ISize) :: &
       IIdentityMatrix

  IIdentityMatrix = 0

  DO ind = 1,ISize
     IIdentityMatrix(ind,ind) = 1
  END DO

END SUBROUTINE CreateIdentityMatrix

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE RecoverSavedSimplex(RSimplexVolume,RSimplexFoM,RStandardDeviation,RMean,IIterationCount,IErr)

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!$  % This Subroutine reads the fr-simplex.txt file from a previous
!!$  % refinement run, and recreates the simplex volume and tolerances
!!$  % allowing for the continuation of a previous refinement
!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI

  IMPLICIT NONE

  INTEGER(IKIND) :: &
       IErr,ind,IIterationCount
  REAL(RKIND),DIMENSION(IIndependentVariables+1,IIndependentVariables) :: &
       RSimplexVolume
  REAL(RKIND),DIMENSION(IIndependentVariables+1) :: &
       RSimplexFoM
  REAL(RKIND) :: &
       RStandardDeviation,RMean
  CHARACTER*200 :: &
       CSizeofData,SFormatString,filename

  WRITE(filename,*) "fr-Simplex.txt"

  OPEN(UNIT=IChOutSimplex,STATUS='UNKNOWN',&
        FILE=TRIM(ADJUSTL(filename)))
  
  WRITE(CSizeofData,*) IIndependentVariables+1
  WRITE(SFormatString,*) "("//TRIM(ADJUSTL(CSizeofData))//"(1F6.3,1X),A1)"

  DO ind = 1,(IIndependentVariables+1)
     READ(IChOutSimplex,FMT=SFormatString) RSimplexVolume(ind,:),RSimplexFoM(ind)
  END DO
    
  READ(IChOutSimplex,FMT="(2(1F6.3,1X),I5.1,I5.1,A1)") RStandardDeviation,RMean,IStandardDeviationCalls,IIterationCount

  CLOSE(IChOutSimplex)

END SUBROUTINE RecoverSavedSimplex

!!$  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE SetupAtomicVectorMovements(IErr)

  USE MyNumbers
  
  USE CConst; USE IConst; USE RConst
  USE IPara; USE RPara; USE SPara; USE CPara
  USE BlochPara

  USE IChannels

  USE MPI
  USE MyMPI
  
  IMPLICIT NONE

  INTEGER(IKIND) :: IErr,ICount,jnd,ind,ISpaceGrp
  INTEGER(IKIND),DIMENSION(:),ALLOCATABLE :: IVectors
  
  CALL ConvertSpaceGroupToNumber(ISpaceGrp,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in ConvertSpaceGroupToNumber"
     RETURN
  ENDIF

  ALLOCATE(IVectors(SIZE(SWyckoffSymbols)),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IVectors"
     RETURN
  ENDIF
  
!XX PRINT*, "Wyckoff Symbols: ",SWyckoffSymbols!XX
  DO ind = 1,SIZE(SWyckoffSymbols)!NB SIZE(SWyckoffSymbols)=IAtomicSitesToRefine
     CALL CountAllowedMovements(ISpaceGrp,SWyckoffSymbols(ind),IVectors(ind),IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine (", my_rank, ") error in CountAllowedMovements "
        RETURN
     ENDIF    
  END DO
 !XX PRINT*, IVectors," IVectors"!XX
  
  IAllowedVectors = SUM(IVectors)
  
  ALLOCATE(IAllowedVectorIDs(IAllowedVectors),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation() of IAllowedVectorIDs"
     RETURN
  ENDIF
  
  ICount = 0
  
  DO ind = 1,SIZE(SWyckoffSymbols)
     DO jnd = 1,IVectors(ind)
        ICount = ICount + 1
        IAllowedVectorIDs(ICount) = IAtomicSitesToRefine(ind)
     END DO
  END DO
  
  ALLOCATE(RAllowedVectors(IAllowedVectors,THREEDIM),STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation() of RAllowedVectors"
     RETURN
  ENDIF
  
  ALLOCATE(RAllowedVectorMagnitudes(IAllowedVectors),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in Allocation() of RAllowedVectorMagnitudes"
     RETURN
  ENDIF
  
  RAllowedVectorMagnitudes = ZERO
  
  DO ind = 1,SIZE(SWyckoffSymbols)
     CALL DetermineAllowedMovements(ISpaceGrp,SWyckoffSymbols(ind),&
          RAllowedVectors(SUM(IVectors(:(ind-1)))+1:SUM(IVectors(:(ind))),:),&
          IVectors(ind),IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"felixrefine (", my_rank, ") error in DetermineAllowedMovements"
        RETURN
     ENDIF
     
  END DO
  
  !--------------------------------------------------------------------
  ! Save Atomic Coordinates  
  !--------------------------------------------------------------------
  
  ALLOCATE(RInitialAtomSiteFracCoordVec(&
       SIZE(RAtomSiteFracCoordVec,DIM=1),&
       SIZE(RAtomSiteFracCoordVec,DIM=2)),&
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"felixrefine (", my_rank, ") error in ALLOCATE() of RInitialAtomSiteFracCoordVec "
     RETURN
  ENDIF
  
  RInitialAtomSiteFracCoordVec = RAtomSiteFracCoordVec
  
END SUBROUTINE SetupAtomicVectorMovements
