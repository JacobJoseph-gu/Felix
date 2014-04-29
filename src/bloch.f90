!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
! FelixSim
!
! Richard Beanland, Keith Evans and Rudolf A Roemer
!
! (C) 2013/14, all rights reserved
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
!  This file is part of FelixSim.
!
!  FelixSim is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!  
!  FelixSim is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!  
!  You should have received a copy of the GNU General Public License
!  along with FelixSim.  If not, see <http://www.gnu.org/licenses/>.
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SUBROUTINE BlochCoefficientCalculation(ind,jnd,IErr)

  USE MyNumbers
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara
  USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  IMPLICIT NONE
  
  INTEGER(IKIND) ind,jnd,hnd,knd,&
       ierr,IThickness, &
       IThicknessIndex, ILowerLimit, &
       IUpperLimit
       
  REAL(RKIND) Rx0,Ry0, RThickness

  CHARACTER*40 surname
  
  COMPLEX(CKIND), DIMENSION(:,:), ALLOCATABLE :: &
       CGeneralSolutionMatrix, CGeneralEigenVectors
  COMPLEX(CKIND),DIMENSION(:),ALLOCATABLE :: &
       CGeneralEigenValues

  Rx0=(ind-IPixelCount-0.5D0)*RDeltaK ! x-position in the disk
  
  Ry0=(jnd-IPixelCount-0.5D0)*RDeltaK ! y-position in the disk
    
  ! we are inside the mask
  IPixelComputed= IPixelComputed + 1

  !--------------------------------------------------------------------
  ! protocol progress
  !--------------------------------------------------------------------
  
  IF(IWriteFLAG.GE.10) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, "): working on pixel (", ind, ",", jnd,") of (", &
          2*IPixelCount, ",", 2*IPixelCount, ") in total."
  ENDIF
       
  !--------------------------------------------------------------------
  ! calculate deviation parameter Sg for the tilted Ewald spheres
  !--------------------------------------------------------------------
  
  ! TiltedK used to be called Kprime2
  ! the vector of the incoming tilted beam

  CALL CalculateKVectors(Rx0,Ry0,IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " In Calculation of KVectors"
     RETURN
  ENDIF

  ! Compute the deviation parameter for ALL reflections
  ! within RBSMaxGVecAmp

  CALL DeviationParameterCalculation(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in Calculation of Deviation Parameter"
     RETURN
  ENDIF

  ! select only those beams where the Ewald sphere is close to the
  ! reciprocal lattice, i.e. within RBSMaxDeviationPara

  CALL DetermineStrongAndWeakBeams(IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in Determination of Strong and Weak beams"
     RETURN
  ENDIF
 
  ! select the highest reflection that corresponds to a strong beam
  nBeams= IStrongBeamIndex

  !--------------------------------------------------------------------
  ! ALLOCATE memory for eigen problem
  !--------------------------------------------------------------------
  
  !Eigen Problem Solving
  ALLOCATE( &
       CBeamProjectionMatrix(nBeams,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CBeamProjectionMatrix"
     RETURN
  ENDIF
  ALLOCATE( &
       CDummyBeamMatrix(nBeams,nReflections), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CDummyBeamMatrix"
     RETURN
  ENDIF
  ALLOCATE( &
       CUgMatEffective(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CUgMatEffective"
     RETURN
  ENDIF
  
  !Allocate General Solution Specific Arrays
  
  IF(IZolzFLAG.EQ.0) THEN
     
     ALLOCATE( &
          CGeneralSolutionMatrix(2*nBeams,2*nBeams), & 
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables General Solution Matrix"
        PRINT*,"Failure Occured at Thickness,ChunkPixel,nBeams = ",IPixelCountTotal,nBeams
        RETURN
     ENDIF
     ALLOCATE( &
          CGeneralEigenVectors(2*nBeams,2*nBeams), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables CEigenVectors"
        RETURN
     ENDIF
     ALLOCATE(&
          CGeneralEigenValues(2*nBeams), &
          STAT=IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
             " in ALLOCATE() of DYNAMIC variables CEigenVectors"
        RETURN
     ENDIF
  END IF
  
  
  ALLOCATE( & 
       CEigenVectors(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CEigenVectors"
     RETURN
  ENDIF

  ALLOCATE( &
       CEigenValues(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CEigenValues"
     RETURN
  ENDIF
  ALLOCATE( &
       CInvertedEigenVectors(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CInvertedEigenVectors"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CAlphaWeightingCoefficients(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CAlphaWeightingCoefficients"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CEigenValueDependentTerms(nBeams,nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CEigenValueDependentTerms"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CWaveFunctions(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CWaveFunctions"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       RWaveIntensity(nBeams), &
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables RWaveIntensity"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     
     RETURN
  ENDIF
  ALLOCATE( &
       CPsi0(nBeams), & 
       STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error ", IErr, &
          " in ALLOCATE() of DYNAMIC variables CPsi0"
     PRINT*,"Failure Occured at Thickness,Chunk,Pixel,nBeams = ",IPixelCountTotal,nBeams
     RETURN
  ENDIF

  !--------------------------------------------------------------------
  ! construct the effective UgMat (for strong beams only at the moment)
  !--------------------------------------------------------------------
  
  IF(IWriteFLAG.GE.10) THEN 
     PRINT*,"BlochCoefficientCalculation(", my_rank, &
          ") using n(Strong)Beams= ", nBeams, " beams", &
          " with nWeakBeams=", IWeakBeamIndex
  ENDIF
  
  !--------------------------------------------------------------------
  ! back to eigen problem solution
  !--------------------------------------------------------------------
  
  ! compute the effective Ug matrix by selecting only those beams
  ! for which IStrongBeamList has an entry
  
  CBeamProjectionMatrix= CZERO
  DO knd=1,nBeams
     CBeamProjectionMatrix(knd,IStrongBeamList(knd))=CONE
  ENDDO
  
  CUgMatEffective= &
       MATMUL( &
       CBeamProjectionMatrix, &
       MATMUL(CUgMat,TRANSPOSE(CBeamProjectionMatrix)) &
       )
  
  IF (IZolzFLAG.EQ.0) THEN
  
     
     ! General Solution from page 457 of 
     ! "Diffraction of Electrons by Perfect Crystals"
     ! by A.J.F.Metherell
     
     !FILL 0 sub-Maxtrix 
     
     CGeneralSolutionMatrix = CZERO

     CGeneralSolutionMatrix(1:nBeams,(nBeams+1):(nBeams*2)) = CUgMatEffective

     DO hnd = 1,nBeams

        ! FILL I sub-Matrix
        
        CGeneralSolutionMatrix(hnd+nBeams,hnd) = CONE

        ! Fill B sub-Matrix
        CGeneralSolutionMatrix(hnd+nBeams,hnd+nBeams) = &
            -2*RgVecMatT(IStrongBeamList(hnd),3) !2*gz Terms
        
        ! Calculate Beta Values for D sub-Matrix
        
        CGeneralSolutionMatrix(hnd,hnd+nBeams) = &
             (RBigK**2 - & !K^2
             ( &
             (RTiltedK(1))**2 + & !kx^2
             (RTiltedK(2))**2 + & !ky^2
             2*(RTiltedK(1))*RgVecMatT(IStrongBeamList(hnd),1) + & !2*kx*gx
             2*(RTiltedK(2))*RgVecMatT(IStrongBeamList(hnd),2) + & !2*ky*gy
             RgVecMatT(IStrongBeamList(hnd),1)**2 + &  !gx^2
             RgVecMatT(IStrongBeamList(hnd),2)**2 + &  !gx^2
             RgVecMatT(IStrongBeamList(hnd),3)**2 & !gx^2
             ))
     END DO
  ELSE
     
     CUgMatEffective = CUgMatEffective/(TWO*RBigK)
     
      ! set the diagonal parts of the matrix to be equal to 
     ! strong beam deviation parameters (*2 BigK) 
     DO hnd=1,nBeams
        CUgMatEffective(hnd,hnd) = RDevPara(IStrongBeamList(hnd))
     ENDDO
     
  END IF

  !PRINT*,"SIZE of CUgMatEffective = ",SIZE(CUgMatEffective,DIM=1),SIZE(CUgMatEffective,DIM=2)
   
  !--------------------------------------------------------------------
  ! diagonalize the UgMatEffective
  !--------------------------------------------------------------------
   
  IF (IZolzFLAG.EQ.0) THEN
     CALL EigenSpectrum(2*nBeams, &
          CGeneralSolutionMatrix, &
          CGeneralEigenValues(:), CGeneralEigenVectors(:,:), &
          IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in EigenSpectrum()"
        RETURN
     ENDIF

     CEigenValues = CGeneralEigenValues((nBeams+1):(2*nBeams))
     CEigenVectors = CGeneralEigenVectors((nBeams+1):(nBeams*2),1:nBeams)

     DEALLOCATE(&
          CGeneralEigenVectors, &
          CGeneralEigenValues, &
          CGeneralSolutionMatrix)
  ELSE
     CALL EigenSpectrum(nBeams, &
          CUgMatEffective, &
          CEigenValues(:), CEigenVectors(:,:), &
          IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in EigenSpectrum()"
        RETURN
     ENDIF
  END IF
 
  DO IThicknessIndex=1,IThicknessCount,1
     
     RThickness = RInitialThickness + (IThicknessIndex-1)*RDeltaThickness 
     IThickness = RInitialThickness + (IThicknessIndex-1)*RDeltaThickness 
     
     CALL CreateWaveFunctions(rthickness,IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in CreateWavefunction()"
        RETURN
     ENDIF

     IF (IOutputFLAG.GE.3) THEN
        IF(IPixelComputed.EQ.1) THEN
           ! wave functions
           CALL OpenData_MPI(IChOutWF_MPI, "WF", surname, IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in OpenDataMPI()"
              RETURN
           ENDIF
        ELSE
           
           ! wave functions
           CALL OpenDataForAppend_MPI(IChOutWF_MPI, "WF", surname, IErr)
           IF( IErr.NE.0 ) THEN
              PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in OpenDataForAppend_ MPI()"
              RETURN
           ENDIF
        END IF
        CALL WriteDataC_MPI(IChOutWF_MPI, ind,jnd, &
             CFullWaveFunctions(:), &
             nReflections, 1, IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in WriteDataC_MPI of IChOutWF()"
           RETURN
        ENDIF
        CALL MPI_FILE_CLOSE(IChOutWF_MPI, IErr)
        IF( IErr.NE.0 ) THEN
           PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in MPI_FILE_CLOSE of IChOutWF()"
           RETURN
        ENDIF
        
     END IF
     
     !Collection Wave Intensities from all thickness for later writing

     IF(IImageFLAG.LE.1) THEN
        RIndividualReflections(ind,jnd,1:IReflectOut,IThicknessIndex) = &
             RFullWaveIntensity(1:IReflectOut)
     ELSE
        CAmplitudeandPhase(ind,jnd,1:IReflectOut,IThicknessIndex) = &
             CFullWavefunctions(1:IReflectOut)
     END IF

  END DO

  
  
  !--------------------------------------------------------------------
  ! OUTPUT EIGENsystem data for given pixel
  !--------------------------------------------------------------------
  
  IMAXCBuffer = 2*13*SIZE(CEigenVectors)+2*13*SIZE(CEigenValues)+3*SIZE(IStrongBeamList)+3*6*ADD_OUT_INFO
  
  IF(IOutputFLAG.GE.1) THEN
     CALL WriteEigenSystem_MPI(IChOutES_MPI, ind,jnd,nReflections,nBeams, &
          CEigenValues,CEigenVectors, IStrongBeamList,nBeams,nBeams, 1, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in WriteEigenSystem_ MPI()"
        RETURN
     ENDIF
  ENDIF
  
  IMAXCBuffer = 2*14*SIZE(CUgMatEffective)+7*6*ADD_OUT_INFO
  
  ! UgMatEffective
  IF(IOutputFLAG.GE.2) THEN
     CALL WriteDataC_MPI(IChOutUM_MPI, ind,jnd, &
          CUgMatEffective(:,:), nBeams*nBeams, 1, IErr)
     IF( IErr.NE.0 ) THEN
        PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in WriteDataC_ MPI() of IChOutUM"
        RETURN
     ENDIF
  ENDIF
  
  
  !--------------------------------------------------------------------
  ! DEALLOCATE eigen problem memory
  !--------------------------------------------------------------------
  
  DEALLOCATE( &
       CUgMatEffective,CPsi0,&
       CInvertedEigenVectors, CAlphaWeightingCoefficients, &
       CEigenValues,CEigenVectors,CEigenValueDependentTerms, &
       CBeamProjectionMatrix, CDummyBeamMatrix,CWavefunctions, &
       RWaveIntensity,STAT=IErr)
  IF( IErr.NE.0 ) THEN
     PRINT*,"BlochCoefficientCalculation(", my_rank, ") error in Deallocation()"
     RETURN
  ENDIF
  
END SUBROUTINE BlochCoefficientCalculation

SUBROUTINE CreateWavefunctions(rthickness,IErr)

 USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  INTEGER ind,jnd,knd,hnd,ierr, ifullind, iuniind,gnd,ichnk
  REAL(RKIND) rthickness 
   
  !--------------------------------------------------------------------
  ! calculate wavefunctions
  !--------------------------------------------------------------------
  
  CPsi0= CZERO
  IF(nBeams .GE. 0) CPsi0(1) = CONE
  
  ! Invert the EigenVector matrix
  CInvertedEigenVectors= CONJG(TRANSPOSE(CEigenVectors(:,:)))
  
  !From EQ 6.32 in Kirkland Advance Computing in EM
  CAlphaWeightingCoefficients = MATMUL(CInvertedEigenVectors(1:nBeams,1:nBeams),CPsi0) 
  
  CEigenValueDependentTerms= ZERO
  
  DO hnd=1,nBeams !IReflectOut 
     
     ! This needs to be a diagonal matrix
     CEigenValueDependentTerms(hnd,hnd) = &
          EXP(CIMAGONE*RThickness*CEigenValues(hnd)) 
     
  ENDDO
  
  ! EQ 6.35 in Kirkland Advance Computing in EM
  ! C-1*C*alpha 
  
  CWaveFunctions(:) = &
       MATMUL( &
       MATMUL(CEigenVectors(1:nBeams,1:nBeams),CEigenValueDependentTerms), & 
       CAlphaWeightingCoefficients(:) &
       )
  
  DO hnd=1,nBeams
     RWaveIntensity(hnd)= &
          CONJG(CWaveFunctions(hnd)) * CWaveFunctions(hnd)
  ENDDO
  
  
  !PRINT*,"This Far"
  
  !--------------------------------------------------------------------
  ! rePADDing of wave function and intensities with zero's 
  !--------------------------------------------------------------------
  
  CFullWaveFunctions=CZERO
  RFullWaveIntensity=ZERO
  
  DO knd=1,nBeams
     CFullWaveFunctions(IStrongBeamList(knd))=CWaveFunctions(knd)
     RFullWaveIntensity(IStrongBeamList(knd))=RWaveIntensity(knd)
  ENDDO
  
END SUBROUTINE CreateWavefunctions

SUBROUTINE CalculateKVectors(Rx0,Ry0,IErr)
  
  USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  IMPLICIT NONE
  
  REAL(RKIND) Rx0,Ry0
  INTEGER(IKIND) :: IErr
  
  RTiltedK(1)= Rx0
  RTiltedK(2)= Ry0
  RTiltedK(3)= SQRT(RBigK**2 - Rx0**2 - Ry0**2)
  
END SUBROUTINE CalculateKVectors

SUBROUTINE DeviationParameterCalculation(IErr)

USE MyNumbers
  
  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI

  INTEGER(IKIND) knd,IErr
  
  DO knd=1,nReflections
     ! DevPara used to be called Sg in the book
     
     RDevPara(knd)= &
          -( RBigK + DOT(RgVecMatT(knd,:),RTiltedK(:)) /RBigK) + &
          SQRT( &
          ( RBigK**2 + DOT(RgVecMatT(knd,:),RTiltedK(:)) )**2 /RBigK**2 - &
          (RgVecMag(knd)**2 + &
          2.0D0* DOT(RgVecMatT(knd,:),RTiltedK(:))) &
          )
  END DO

END SUBROUTINE DeviationParameterCalculation

SUBROUTINE DetermineStrongAndWeakBeams(IErr)

  USE CConst; USE IConst
  USE IPara; USE RPara; USE CPara; USE SPara
  USE IChannels
  USE BlochPara
  
  USE MPI
  USE MyMPI
  
  INTEGER(IKIND) ind,knd,IErr,IMinimum,IMaximum,ICheck,jnd
  REAL(RKIND) RDummySg(nReflections)

  !Determine RBSMaxDeviationPara

  RDummySg = ABS(RDevPara)

  DO ind=1,IMinStrongBeams
     IMinimum = MINLOC(RDummySg,1)
     IF(ind.EQ.IMinStrongBeams) THEN
        RBSMaxDeviationPara = ABS(RDummySg(IMinimum))
     ELSE
        RDummySg(IMinimum) = 1000000 !Large number
     END IF
  END DO
  
  IStrongBeamIndex=0
  IWeakBeamIndex=0
  DO knd=1,nReflections
     IF( ABS(RDevPara(knd)) .LE. RBSMaxDeviationPara ) THEN
        IStrongBeamIndex= IStrongBeamIndex +1
        IStrongBeamList(IStrongBeamIndex)= knd
     ENDIF
  ENDDO
  
  RDummySg = ABS(RMeanInnerCrystalPotential/RDevPara)
  
  jnd = 0

  !Determine RBSBethePara

  DO ind=1,nReflections
     ICheck = 0
     IMaximum = MAXLOC(RDummySg,1)

     DO knd = 1,IStrongBeamIndex
        IF(IMaximum.EQ.IStrongBeamList(knd)) THEN
           ICheck = 1
           EXIT
        END IF
     END DO

     IF(ICheck.EQ.0) THEN
        jnd = jnd+1
     END IF

     IF(jnd.EQ.IMinWeakBeams) THEN
        RBSBethePara = (RDummySg(IMaximum))
     ELSE
        RDummySg(IMaximum) = 0.D0 !Large number
     END IF
  END DO

  IWeakBeamIndex=0
  DO knd=1,nReflections
     IF( (ABS(RDevPara(knd)) .GT. RBSMaxDeviationPara).AND. &
          (ABS(RMeanInnerCrystalPotential/RDevPara(knd)) .GE. RBSBethePara) ) THEN
        IWeakBeamIndex= IWeakBeamIndex +1
        IWeakBeamList(IWeakBeamIndex)= knd
     ENDIF
  ENDDO

END SUBROUTINE DetermineStrongAndWeakBeams