SUBROUTINE MultilayerT_PBC_fct_constL(fileNum,DFT_path,write_path, matlsInds,Ls,Q0,skD,numPer, dT, kCond)
!Calculate deviational temperature and heat flux 
!through super lattice, assuming blackbody boundary conditions with dT1 and
!dT2 as temperatures at either end, and specularity spec at each interface

! periodic system of 4 layers;  
! Ls(1) = thickness of physical layer;  Ls(2) = thickness of layer represented at interface
! in this piece of code, our structure is L1/int/L1/int/;  the interface can be anything, but we'll make it be a QWP of length L2

USE Constants
USE omp_lib
IMPLICIT NONE

INTERFACE
    FUNCTION GetMatlParams(DFT_path) RESULT(allMatls)
		USE Constants
		IMPLICIT NONE
		
		! REAL, INTENT(IN) :: roughness
		CHARACTER(longStr), INTENT(IN) :: DFT_path
        TYPE(material), DIMENSION(numMatls) :: allMatls
    END FUNCTION
	
	SUBROUTINE GetRTMats(matlParams, Ls, write_path, &
					R12_d,T12_d,R21_d,T21_d,g1mCoeff_d,g2pCoeff_d, &
					R12_s,T12_s,R21_s,T21_s,g1mCoeff_s,g2pCoeff_s, &
					numG12d,numG21d,numR12d,numR21d,numT12d,numT21d, &
					numG12s,numG21s,numR12s,numR21s,numT12s,numT21s, fBoth,numFBoth,weightsMatch, spec1, spec2)
	! these matrices are represented by a list of indices and their values (row, col, val)
		USE Constants
		IMPLICIT NONE

		TYPE(material), DIMENSION(2), INTENT(IN) :: matlParams
		REAL, DIMENSION(2), INTENT(IN) :: Ls
		CHARACTER(longStr), INTENT(IN) :: write_path
		TYPE(sIndVal), INTENT(OUT) :: R12_d(:), T12_d(:), R21_d(:), T21_d(:), g1mCoeff_d(:), g2pCoeff_d(:), &
							R12_s(:), T12_s(:), R21_s(:), T21_s(:), g1mCoeff_s(:), g2pCoeff_s(:)
		INTEGER, INTENT(OUT) :: numG12d,numG21d,numR12d,numR21d,numT12d,numT21d,&
								numG12s,numG21s,numR12s,numR21s,numT12s,numT21s, numFBoth,fBoth(:,:)
		REAL *8, INTENT(OUT) :: spec1(:), spec2(:), weightsMatch(:)
	END SUBROUTINE
	
	SUBROUTINE GetInterfaceBM_PBC(matlParams,Ls,Qs,skD, spec1, spec2, dT12, &
							G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
							MatABCD, MatEFGH, VecJK, Amat1, f11Coeffs, f12Coeffs, f13Coeffs, Amat2, f21Coeffs, f22Coeffs, f23Coeffs)
		USE Constants
		USE omp_lib
		IMPLICIT NONE
		
		TYPE(material), DIMENSION(2), INTENT(IN) :: matlParams
					! material parameters from DFT for 2 adjacent blocks			
		REAL, DIMENSION(2), INTENT(IN) :: Qs, Ls			
					! heat source magnitude from laser at surface of slab, length of slabs
					! for the two slabs of interest
		REAL, INTENT(IN) :: skD, spec1(:,:), spec2(:,:), dT12
					! specularity for each mode
		TYPE(sIndVal), DIMENSION(:), INTENT(IN) :: R12_d, T12_d, T21_d, R21_d, R12_s, T12_s, T21_s, R21_s, G12_d, G12_s, G21_d, G21_s
		REAL*8, INTENT(OUT) :: MatABCD(:,:), MatEFGH(:,:), VecJK(:,:), &
								Amat1(:,:), f11Coeffs(:,:), f12Coeffs(:,:), f13Coeffs(:,:), &
								Amat2(:,:), f21Coeffs(:,:), f22Coeffs(:,:), f23Coeffs(:,:)
	END SUBROUTINE
	
	SUBROUTINE CheckBCs_Pbc(fileNum,write_path, matlParams, dT21, Qs, skD, Ls, totL, numPer, &
						 fX1Coeffs, fX2Coeffs, B1wvecs_s, B1wvecs_d, P2wvecs_s, P2wvecs_d, &
						 G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, spec1, spec2, &
						 fBoth, numFBoth, weightsMatch, &
						 dTCoeffs_s, dTCoeffs_d, kCond)

		USE Constants
		IMPLICIT NONE
								 
		CHARACTER(1), INTENT(IN) :: fileNum
		CHARACTER(longStr), INTENT(IN) :: write_path
		TYPE(material), DIMENSION(2), INTENT(IN) :: matlParams
		REAL, INTENT(IN) :: Ls(:), Qs(:), skD, totL, dT21
		REAL*8, DIMENSION(:,:,:), INTENT(IN) :: fX1Coeffs, fX2Coeffs
		REAL*8, DIMENSION(:,:,:), INTENT(IN) :: B1wvecs_s, B1wvecs_d, P2wvecs_s, P2wvecs_d
		REAL*8, DIMENSION(:,:), INTENT(IN) :: dTCoeffs_d, dTCoeffs_s, spec1, spec2
		TYPE(sIndVal), DIMENSION(:), INTENT(IN) :: R12_d, T12_d, T21_d, R21_d, R12_s, T12_s, T21_s, R21_s, G12_d, G12_s, G21_d, G21_s
		INTEGER, INTENT(IN) :: fBoth(:,:), numFBoth, numPer
		REAL*8, DIMENSION(:), INTENT(IN) :: weightsMatch
		REAL*8, DIMENSION(:), INTENT(OUT) :: kCond
	END SUBROUTINE
END INTERFACE



INTEGER, INTENT(IN) :: matlsInds(:), numPer
CHARACTER(1), INTENT(IN) :: fileNum
CHARACTER(longStr), INTENT(IN) :: DFT_path, write_path
REAL*8, INTENT(IN) :: Q0, skD, Ls(:)
REAL*8, DIMENSION(:), INTENT(OUT) :: kCond


TYPE(material), DIMENSION(numMatls) :: matlParams

INTEGER :: nL, maxNumRows, m, m1,m2, matN, dummy1, dummy2
REAL, DIMENSION(:), ALLOCATABLE :: Qs
INTEGER, DIMENSION(:), ALLOCATABLE :: numW, numWd
REAL*8, DIMENSION(:,:), ALLOCATABLE :: spec1, spec2, weights1, weights2

INTEGER :: map11Max,map22Max, numG12s,numR12s,numT12s,numG21s,numR21s,numT21s, numG12d,numR12d,numT12d,numG21d,numR21d,numT21d, numFBoth
CHARACTER(shortStr) :: TempString
CHARACTER(longStr) :: RT_fileRoot
LOGICAL :: filesExist
TYPE(sIndVal), ALLOCATABLE, DIMENSION(:) :: G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
								G12_st, R12_st, T12_st, G21_st, R21_st, T21_st, G12_dt, R12_dt, T12_dt, G21_dt, R21_dt, T21_dt
INTEGER, ALLOCATABLE, DIMENSION(:,:) :: fBoth
REAL*8, ALLOCATABLE, DIMENSION(:) :: weightsMatch

INTEGER :: w0,w1,w2,w3, w0d,w1d,w2d,w3d, wMax,wMaxd, numRows,numRows1,numRows0
REAL*8, DIMENSION(:,:), ALLOCATABLE :: P1wvec, BEndwvec, matErr
										! all are column vectors
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: B1P2s, B1wvecs, P2wvecs, tempB1_d, tempP2_d, &
										B1wvecs_s, B1wvecs_d, P2wvecs_s, P2wvecs_d

REAL*8, DIMENSION(:,:), ALLOCATABLE :: MatABCD_t1, MatEFGH_t1, VecJK_t1, MatABCD_t0, MatEFGH_t0, VecJK_t0, &
											temp1, temp2
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: AmatX, fX1Coeffs, fX2Coeffs, fX3Coeffs, fX1Coeffs_d, fX2Coeffs_d


REAL*8, DIMENSION(:,:), ALLOCATABLE :: dTCoeffs, dTCoeffs_s, dTCoeffs_d
REAL*8, DIMENSION(:), ALLOCATABLE :: uniqueFreqs1, uniqueFreqs2
INTEGER, DIMENSION(:), ALLOCATABLE :: numEachFreq1, numEachFreq2
REAL*8, ALLOCATABLE :: bcErr(:,:), intGs(:,:,:), intQs(:,:,:)
REAL, ALLOCATABLE :: cosMat(:,:), xs(:)
REAL*8, DIMENSION(:,:), ALLOCATABLE :: dTtemp
CHARACTER(10) :: timeT

REAL*8, DIMENSION(:,:), INTENT(OUT) :: dT


PRINT *, 'inside multilayer T fct PBC'
PRINT *, 'T1', dT1, 'T2', dT2
PRINT *, 'Ls shape', SHAPE(Ls)
PRINT *, 'Ls', Ls
PRINT *, 'matlInds', matlsInds
PRINT *, 'rms roughness', roughness
PRINT *, 'fileNum', fileNum

CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'start time: ', timeT
nL = SIZE(matlsInds,1)*numPer  ! numLayers
PRINT *, 'nL', nL
PRINT *, 'const L'

ALLOCATE(Qs(1:nL), numW(0:nL+1), numWd(0:nL+1))
! Qs(1) = Q0;
Qs = 0;

matlParams = GetMatlParams(DFT_path);
PRINT *, 'got matl params'

DO m = 1,nL
	m1 = 2-MODULO(m,2);
	numW(m) = matlParams(matlsInds(m1))%wTot;
	numWd(m) = matlParams(matlsInds(m1))%numUniqueFreq;
END DO
numW(0) = numW(2); numWd(0) = numWd(2);
numW(nL+1) = numW(1); numWd(nL+1) = numWd(1);

wMax = MAXVAL(numW); wMaxd = MAXVAL(numWd);



!!!!!!!!!!!!!!! R and T matrices !!!!!!!!!!!!!!!!!!!
! Get R/T matrices for material 1/material 2 interface (and then we know matrices for matl2/matl1 interface)
w1 = numW(1); w1d = numWd(1);  w2 = numW(2); w2d = numWd(2); 
PRINT *, 'int bm:  w1', w1, 'w2', w2, 'w1d', w1d, 'w2d', w2d

ALLOCATE(numEachFreq1(w1d),numEachFreq2(w2d),uniqueFreqs1(w1d),uniqueFreqs2(w2d))
		 
numEachFreq1 = matlParams(matlsInds(1))%numEachFreq(1:w1d);
numEachFreq2 = matlParams(matlsInds(2))%numEachFreq(1:w2d);
uniqueFreqs1 = matlParams(matlsInds(1))%uniqueFreqs(1:w1d);
uniqueFreqs2 = matlParams(matlsInds(2))%uniqueFreqs(1:w2d);

! get temporary transmission matrices
map11Max = SUM(numEachFreq1**2); map22Max = SUM(numEachFreq2**2)
PRINT *, 'maps', map11Max, map22Max
ALLOCATE(R12_dt(map11Max), T12_dt(MAX(map11Max,map22Max)*2), T21_dt(MAX(map11Max,map22Max)*2), &
		 R21_dt(map22Max), G12_dt(map11Max), G21_dt(map22Max), &
		 R12_st(map11Max), T12_st(MAX(map11Max,map22Max)), T21_st(MAX(map11Max,map22Max)), &
		 R21_st(map22Max), G12_st(map11Max), G21_st(map22Max))
		! *2 in T12_dt, T21_dt defn to allow for cross over of frequencies

ALLOCATE(spec1(numW(1),1),spec2(numW(2),1))
ALLOCATE(fBoth(w1d*w2d,2),weightsMatch(w1d*w2d))


WRITE(TempString, '(i4)') T
RT_fileRoot = 'RT_' // TRIM(matlParams(matlsInds(1))%matlName) // TRIM(matlParams(matlsInds(2))%matlName) &
								// '_' // TRIM(ADJUSTL(TempString)) // '_'

INQUIRE(FILE = TRIM(RT_fileRoot) // 'R12s.dat', EXIST=filesExist)
! filesExist = .FALSE.
IF (matlsInds(1)==matlsInds(2)) filesExist = .FALSE.

IF (filesExist) THEN		! read the files

	! !$OMP PARALLEL PRIVATE(m,dummy1, dummy2)

	! !$OMP SINGLE
	WRITE (*,*) 'Parallel part num threads: ', OMP_GET_NUM_THREADS()
	! !$OMP END SINGLE

	! !$OMP SECTIONS
	! !$OMP SECTION

		OPEN(1,FILE = TRIM(RT_fileRoot) // 'R12s.dat')
		READ(1,*) dummy1, dummy2, numR12s
		DO m = 1,numR12s
			READ(1,*) R12_st(m)%row, R12_st(m)%col, R12_st(m)%indVal
		END DO
		CLOSE(1)
		
	! !$OMP SECTION

		OPEN(2,FILE = TRIM(RT_fileRoot) // 'R12d.dat')
		READ(2,*) dummy1, dummy2, numR12d
		DO m = 1,numR12d
			READ(2,*) R12_dt(m)%row, R12_dt(m)%col, R12_dt(m)%indVal
		END DO
		CLOSE(2)
		
	! !$OMP SECTION

		OPEN(3,FILE = TRIM(RT_fileRoot) // 'T12s.dat')
		READ(3,*) dummy1, dummy2, numT12s
		DO m = 1,numT12s
			READ(3,*) T12_st(m)%row, T12_st(m)%col, T12_st(m)%indVal
		END DO
		CLOSE(3)
		
	! !$OMP SECTION
		OPEN(4,FILE = TRIM(RT_fileRoot) // 'T12d.dat')
		READ(4,*) dummy1, dummy2, numT12d
		DO m = 1,numT12d
			READ(4,*) T12_dt(m)%row, T12_dt(m)%col, T12_dt(m)%indVal
		END DO
		CLOSE(4)
		
	! !$OMP SECTION
		OPEN(5,FILE = TRIM(RT_fileRoot) // 'G12s.dat')
		READ(5,*) dummy1, dummy2, numG12s
		DO m = 1,numG12s
			READ(5,*) G12_st(m)%row, G12_st(m)%col, G12_st(m)%indVal
		END DO
		CLOSE(5)
		
	! !$OMP SECTION
		OPEN(7,FILE = TRIM(RT_fileRoot) // 'G12d.dat')
		READ(7,*) dummy1, dummy2, numG12d
		DO m = 1,numG12d
			READ(7,*) G12_dt(m)%row, G12_dt(m)%col, G12_dt(m)%indVal
		END DO
		CLOSE(7)
		
	! !$OMP SECTION
		OPEN(8,FILE = TRIM(RT_fileRoot) // 'R21s.dat')
		READ(8,*) dummy1, dummy2, numR21s
		DO m = 1,numR21s
			READ(8,*) R21_st(m)%row, R21_st(m)%col, R21_st(m)%indVal
		END DO
		CLOSE(8)
		
	! !$OMP SECTION
		OPEN(9,FILE = TRIM(RT_fileRoot) // 'R21d.dat')
		READ(9,*) dummy1, dummy2, numR21d
		DO m = 1,numR21d
			READ(9,*) R21_dt(m)%row, R21_dt(m)%col, R21_dt(m)%indVal
		END DO
		CLOSE(9)
		
	! !$OMP SECTION
		OPEN(10,FILE = TRIM(RT_fileRoot) // 'T21s.dat')
		READ(10,*) dummy1, dummy2, numT21s
		DO m = 1,numT21s
			READ(10,*) T21_st(m)%row, T21_st(m)%col, T21_st(m)%indVal
		END DO
		CLOSE(10)
		
	! !$OMP SECTION
		OPEN(11,FILE = TRIM(RT_fileRoot) // 'T21d.dat')
		READ(11,*) dummy1, dummy2, numT21d
		DO m = 1,numT21d
			READ(11,*) T21_dt(m)%row, T21_dt(m)%col, T21_dt(m)%indVal
		END DO
		CLOSE(11)
		
	! !$OMP SECTION
		OPEN(12,FILE = TRIM(RT_fileRoot) // 'G21s.dat')
		READ(12,*) dummy1, dummy2, numG21s
		DO m = 1,numG21s
			READ(12,*) G21_st(m)%row, G21_st(m)%col, G21_st(m)%indVal
		END DO
		CLOSE(12)
		
	! !$OMP SECTION
		OPEN(13,FILE = TRIM(RT_fileRoot) // 'G21d.dat')
		READ(13,*) dummy1, dummy2, numG21d
		DO m = 1,numG21d
			READ(13,*) G21_dt(m)%row, G21_dt(m)%col, G21_dt(m)%indVal
		END DO
		CLOSE(13)

	! !$OMP SECTION
		OPEN(14,FILE = TRIM(RT_fileRoot) // 'spec1.dat')
		! DO m = 1,w1
			READ(14,*) spec1(:,1)
		! END DO
		CLOSE(14)
		
	! !$OMP SECTION
		OPEN(15,FILE = TRIM(RT_fileRoot) // 'spec2.dat')
		! DO m = 1,w2
			READ(15,*) spec2(:,1)
		! END DO
		CLOSE(15)
		
	! !$OMP SECTION
		OPEN(16,FILE = TRIM(RT_fileRoot) // 'fMatch.dat')
		READ(16,*) numFBoth
		DO m = 1,numFBoth
			READ(16,*) fBoth(m,1), fBoth(m,2), weightsMatch(m)
		END DO
		CLOSE(16)
		
	! !$OMP END SECTIONS
	! !$OMP END PARALLEL
	
ELSE 
	PRINT *, 'solving for matrices...'

	
	CALL GetRTMats((/ matlParams(matlsInds(1)),matlParams(matlsInds(2)) /), Ls, RT_fileRoot, &
					R12_dt,T12_dt,R21_dt,T21_dt,G12_dt,G21_dt, &
					R12_st,T12_st,R21_st,T21_st,G12_st,G21_st, &
					numG12d,numG21d,numR12d,numR21d,numT12d,numT21d,&
					numG12s,numG21s,numR12s,numR21s,numT12s,numT21s,fBoth,numFBoth,weightsMatch,spec1(:,1),spec2(:,1))
					
END IF


ALLOCATE(R12_d(numR12d), T12_d(numT12d), T21_d(numT21d), R21_d(numR21d), G12_d(numG12d), G21_d(numG21d), &
		 R12_s(numR12s), T12_s(numT12s), T21_s(numT21s), R21_s(numR21s), G12_s(numG12s), G21_s(numG21s))

R12_d = R12_dt(1:numR12d); T12_d = T12_dt(1:numT12d); T21_d = T21_dt(1:numT21d); R21_d = R21_dt(1:numR21d)
G12_d = G12_dt(1:numG12d); G21_d = G21_dt(1:numG21d)
R12_s = R12_st(1:numR12s); T12_s = T12_st(1:numT12s); T21_s = T21_st(1:numT21s); R21_s = R21_st(1:numR21s)
G12_s = G12_st(1:numG12s); G21_s = G21_st(1:numG21s)
				
DEALLOCATE(R12_dt,T12_dt,T21_dt,R21_dt,R12_st,T12_st,T21_st,R21_st)

PRINT *, 'got RTMats'

DEALLOCATE(numEachFreq1,numEachFreq2,uniqueFreqs1,uniqueFreqs2)

PRINT *, 'successful RTMat call' 




!!!!!!!!!! solving for B1P2s !!!!!!!!!!!!!!!!!!111



ALLOCATE(B1P2s((wMax+wMaxd)*2,1,2),B1wvecs_s(wMax,1,2),P2wvecs_s(wMax,1,2), &
		 tempB1_d(wMaxd,1,2), tempP2_d(wMaxd,1,2), B1wvecs_d(wMax,1,2), P2wvecs_d(wMax,1,2), &
		 B1wvecs(wMax,1,2), P2wvecs(wMax,1,2))
		 ! B1wvecs_d, P2wvecs_d contain data from 1:w1, 1:w2
		 ! it's expanded from tempB1_d, tempP2_d from 1:w1d, 1:w2d

B1P2s = 0.0; B1wvecs_s = 0.0; P2wvecs_s = 0.0; tempB1_d = 0.0; tempP2_d = 0.0; 
B1wvecs_d = 0.0; P2wvecs_d = 0.0; B1wvecs = 0.0; P2wvecs = 0.0;




!!!!!!!!!!! get MatABCD, MatEF, MatGH, VecJK !!!!!!!!!!!!!!!!!!


! layer 1/2
! get block matrices for first interface (interface 1)
w1 = numW(1); w1d = numWd(1)
w2 = numW(2); w2d = numWd(2)

numRows = w1+w2+w1d+w2d;	! number of rows for given interface 
numRows0 = w1+w1d;	! number of nonzero cols for MatGH
numRows1 = w2+w2d;	! number of nonzero cols for MatEF

PRINT *, 'Calling block Mats'
CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'start time block Mats: ', timeT
PRINT *, 'test', EXP(-Ls(1)/skD)

! ALLOCATE(MatABCD_t1(numRows,numRows), MatEF_t1(numRows,numRows1), MatGH_t1(numRows,numRows0), VecJK_t1(numRows,1), &
		 ! MatABCD_t0(numRows,numRows), MatEF_t0(numRows,numRows0), MatGH_t0(numRows,numRows1), VecJK_t0(numRows,1), &
		 ! AmatX(N+1,N+1,2),fX1Coeffs(N+1,wMax,2),fX2Coeffs(N+1,wMax,2),fX3Coeffs(N+1,1,2))
ALLOCATE(MatABCD_t1(numRows,numRows), MatEFGH_t1(numRows,numRows), VecJK_t1(numRows,1), &
		 MatABCD_t0(numRows,numRows), MatEFGH_t0(numRows,numRows), VecJK_t0(numRows,1), &
		 AmatX(N+1,N+1,2),fX1Coeffs(N+1,wMax,2),fX2Coeffs(N+1,wMax,2),fX3Coeffs(N+1,1,2))
		 
MatABCD_t1 = 0; MatEFGH_t1 = 0; VecJK_t1 = 0;
MatABCD_t0 = 0; MatEFGH_t0 = 0; VecJK_t0 = 0;
AmatX = 0; fX1Coeffs = 0; fX2Coeffs = 0; fx3Coeffs = 0;


! layer 1/2
CALL GetInterfaceBM_PBC((/ matlParams(matlsInds(1)),matlParams(matlsInds(2)) /),(/Ls(1),Ls(1)/),Qs(1:2),skD, spec1, spec2, 0.0, &
		G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
		MatABCD_t1, MatEFGH_t1, VecJK_t1, &
		AmatX(1:N+1,1:N+1,1), fX1Coeffs(1:N+1,1:w1,1), fX2Coeffs(1:N+1,1:w1,1), fX3Coeffs(1:N+1,:,1), &
		AmatX(1:N+1,1:N+1,2), fX1Coeffs(1:N+1,1:w2,2), fX2Coeffs(1:N+1,1:w2,2), fX3Coeffs(1:N+1,:,2))

! layer 2/3
CALL GetInterfaceBM_PBC((/ matlParams(matlsInds(2)),matlParams(matlsInds(1)) /),(/Ls(1),Ls(1)/),Qs(2:1:-1),skD, spec2, spec1, dT2-dT1, &
		G21_s, R21_s, T21_s, G12_s, R12_s, T12_s, G21_d, R21_d, T21_d, G12_d, R12_d, T12_d, &
		MatABCD_t0, MatEFGH_t0, VecJK_t0, &
		AmatX(1:N+1,1:N+1,2), fX1Coeffs(1:N+1,1:w2,2), fX2Coeffs(1:N+1,1:w2,2), fX3Coeffs(1:N+1,:,2), &
		AmatX(1:N+1,1:N+1,1), fX1Coeffs(1:N+1,1:w1,1), fX2Coeffs(1:N+1,1:w1,1), fX3Coeffs(1:N+1,:,1))	

		
CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'end block mats call time: ', timeT


PRINT *, 'Defining B1P2s'
ALLOCATE(temp1(numRows,numRows),temp2(numRows,1))
temp1 = MatABCD_t1 - MATMUL(MatEFGH_t1, cAXB(MatABCD_t0,MatEFGH_t0));
temp2 = MATMUL(MatEFGH_t1,cAXB(MatABCD_t0,VecJK_t0))+VecJK_t1;
B1P2s(:,:,1) = cAXB(temp1,temp2)
	! 1/2 interface
B1P2s(:,:,2) = cAXB(MatABCD_t0,(MATMUL(MatEFGH_t0,B1P2s(:,:,1)) + VecJK_t0));
	! 2/3 (or 0/1) interface
DEALLOCATE(temp1,temp2)


DO m = 1,nL
	m1 = 2-MODULO(m,2)		! will only have two layers
	m2 = 2-MODULO(m+1,2)
	! 1 = mat1; 2 = mat2;

	w1 = numW(m1); w1d = numWd(m1); 
	w2 = numW(m2); w2d = numWd(m2);
		
	 
	B1wvecs_s(1:w1,1,m) = B1P2s(1:w1,1,m); 
	PRINT *, 'B1wvecs_s done'
	P2wvecs_s(1:w2,1,m) = B1P2s((w1+1):(w1+w2),1,m);
	PRINT *, 'P2wvecs_s done'
	tempB1_d(1:w1d,1,m) = B1P2s((w1+w2+1):(w1+w2+w1d),1,m);
	PRINT *, 'tempB1_d done'
	tempP2_d(1:w2d,1,m) = B1P2s((w1+w2+w1d+1):(w1+w1d+w2+w2d),1,m);
	PRINT *, 'tempP2_d done'

	B1wvecs_d(1:w1,1,m) = cDIFFVECtoFULL(tempB1_d(1:w1d,1,m), w1, matlParams(matlsInds(m1))%numEachFreq(1:w1d));
	P2wvecs_d(1:w2,1,m) = cDIFFVECtoFULL(tempP2_d(1:w2d,1,m), w2, matlParams(matlsInds(m2))%numEachFreq(1:w2d));

	B1wvecs(1:w1,1,m) = B1wvecs_s(1:w1,1,m) + B1wvecs_d(1:w1,1,m)
	P2wvecs(1:w2,1,m) = P2wvecs_s(1:w2,1,m) + P2wvecs_d(1:w2,1,m)
END DO

PRINT *, 'finished solving for b1p2'



!!!!!!!!!!!! calculate matrix error !!!!!!!!!!!!!!!	
ALLOCATE(matErr(numRows,1))

matErr = MATMUL(MatABCD_t1,B1P2s(:,:,1)) - MATMUL(MatEFGH_t1,B1P2s(:,:,2)) - VecJK_t1;

OPEN(1,FILE = TRIM(write_path) // 'output_matErr_' // fileNum // '.dat')
WRITE(1,"(40000(ES15.7,','))") REAL(matErr(:,1))
CLOSE(1)


matErr = MATMUL(MatABCD_t0,B1P2s(:,:,2)) - MATMUL(MatEFGH_t0,B1P2s(:,:,1)) - VecJK_t0;

! append to existing file
PRINT *, 'append?'
OPEN(1,FILE = TRIM(write_path) // 'output_matErr_' // fileNum // '.dat', POSITION='append')
WRITE(1,"(40000(ES15.7,','))") REAL(matErr(:,1))
CLOSE(1)


DEALLOCATE(matErr)
DEALLOCATE(MatABCD_t1,MatEFGH_t1,VecJK_t1, MatABCD_t0,MatEFGH_t0,VecJK_t0)
	


!!!!!!!!!!!!!!! Calculating dT !!!!!!!!!!!!!

PRINT *, 'Calculating dT'
CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'start time calc dT: ', timeT
	
! Calculation of deviational temperature distribution
ALLOCATE(dTCoeffs(N+1,nL), dTCoeffs_s(N+1,nL), dTCoeffs_d(N+1,nL), cosMat(Nx,N+1), xs(Nx))

PRINT *, 'size dTcoeffs', SHAPE(dTCoeffs)
dTCoeffs = 0;
DO m = 1,nL
	PRINT *, 'calcT loop', m
	w1 = numW(m)
	w1d = numWd(m)
	m1 = 2-MODULO(m-1,2);
	m2 = 2-MODULO(m,2);  ! m odd = layer 1 of unit cell, m even = layer 2 of unit cell
	dTCoeffs(:,m) = RESHAPE(cAXB(AmatX(:,:,m2), &
						(MATMUL(fX1Coeffs(:,1:w1,m2), P2wvecs(1:w1,:,m1)) + &
						 MATMUL(fX2Coeffs(:,1:w1,m2), B1wvecs(1:w1,:,m2)) + &
						fX3Coeffs(:,:,m2))),(/w1/))
	dTCoeffs_s(:,m) = RESHAPE(cAXB(AmatX(:,:,m2), &
						(MATMUL(fX1Coeffs(:,1:w1,m2), P2wvecs_s(1:w1,:,m1)) + &
						 MATMUL(fX2Coeffs(:,1:w1,m2), B1wvecs_s(1:w1,:,m2)) + &
						fX3Coeffs(:,:,m2))),(/w1/))
	dTCoeffs_d(:,m) = RESHAPE(cAXB(AmatX(:,:,m2), &
						(MATMUL(fX1Coeffs(:,1:w1,m2), P2wvecs_d(1:w1,:,m1)) + &
						 MATMUL(fX2Coeffs(:,1:w1,m2), B1wvecs_d(1:w1,:,m2)) + &
						fX3Coeffs(:,:,m2))),(/w1/))
END DO 



!!!!!!!!!!!!!!!! calculating dT !!!!!!!!!!!!!!!!!

xs = (/ (m*(1.0/(Nx-1)), m = 0,Nx-1) /)

DO m = 0,N
	cosMat(:,m+1) = COS(m*pi*xs)
END DO

ALLOCATE( dTtemp(SIZE(dT,1), SIZE(dT,2)) )
PRINT *, 'size dT, dTtemp', SHAPE(dTtemp)
dTtemp = MATMUL(cosMat,dTCoeffs)

CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'end time calc dT: ', timeT
PRINT *, 'writing file'

! write dT results to a txt file
OPEN(1,FILE = TRIM(write_path) // 'output_dT_' // fileNum // '.txt')
DO m = 1,SIZE(dTtemp,1)
	WRITE(1,"(10(ES15.7,','))") xs(m), REAL(dTtemp(m,:))
END DO
CLOSE(1)



CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'end time writing files: ', timeT
PRINT *, 'going on...'

PRINT *, 'calling CheckBCs'
PRINT *, Ls, skD
PRINT *, 'test', EXP(-Ls(1)/skD)


OPEN(2,FILE = TRIM(write_path) // 'output_bcErrS_' // fileNum // '.dat')
OPEN(3,FILE = TRIM(write_path) // 'output_bcErrD_' // fileNum // '.dat')
OPEN(4,FILE = TRIM(write_path) // 'output_bcErr_' // fileNum // '.dat')
OPEN(5,FILE = TRIM(write_path) // 'output_specErr_' // fileNum // '.dat')

DO m = 1,nL
	m1 = 2-MODULO(m,2)		! will only have two layers
	m2 = 2-MODULO(m+1,2)
	IF (m1==1) THEN	! mat 1 -> mat 2
		
		! bc for 1/2 interface		
		CALL CheckBCs_Pbc( fileNum, write_path, (/ matlParams(matlsInds(1)),matlParams(matlsInds(2)) /), &
							 0.0, Qs(1:2), skD, (/Ls(1),Ls(1)/), Ls(1)*2*numPer, numPer, &
							 fX1Coeffs(:,:,(/1,2/)), fX2Coeffs(:,:,(/1,2/)), &
							 B1wvecs_s(:,:,(/2,1,2/)), B1wvecs_d(:,:,(/2,1,2/)), P2wvecs_s(:,:,(/2,1,2/)), P2wvecs_d(:,:,(/2,1,2/)), &
							 G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, spec1, spec2, &
							 fBoth(:,(/1,2/)), numFBoth, weightsMatch, &
							 dTCoeffs_s(:,1:2), dTCoeffs_d(:,1:2), kCond(1:1))
	ELSE
		! bc for 2/3 (0/1) interface
		CALL CheckBCs_Pbc( fileNum, write_path, (/ matlParams(matlsInds(2)),matlParams(matlsInds(1)) /), &
							 dT2-dT1, Qs((/2,1/)), skD, (/Ls(1),Ls(1)/), Ls(1)*2*numPer, numPer, &
							 fX1Coeffs(:,:,(/2,1/)), fX2Coeffs(:,:,(/2,1/)), &
							 B1wvecs_s(:,:,(/1,2,1/)), B1wvecs_d(:,:,(/1,2,1/)), P2wvecs_s(:,:,(/1,2,1/)), P2wvecs_d(:,:,(/1,2,1/)), &
							 G21_s, R21_s, T21_s, G12_s, R12_s, T12_s, G21_d, R21_d, T21_d, G12_d, R12_d, T12_d, spec2, spec1, &
							 fBoth(:,(/2,1/)), numFBoth, weightsMatch, &
							 dTCoeffs_s(:,(/2,1/)), dTCoeffs_d(:,(/2,1/)), kCond(2:2))
	END IF
END DO		
				 
PRINT *, 'interface conductance', kCond			 

CLOSE(2)
CLOSE(3)
CLOSE(4)
CLOSE(5)

CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'end time checkbcs: ', timeT

DEALLOCATE(spec1,spec2, numW, numWd, Qs)
DEALLOCATE(fBoth,weightsMatch)
DEALLOCATE(B1P2s, B1wvecs, P2wvecs, tempB1_d, tempP2_d, B1wvecs_s, B1wvecs_d, P2wvecs_s, P2wvecs_d)
DEALLOCATE(AmatX, fX1Coeffs, fX2Coeffs, fX3Coeffs, dTCoeffs, dTCoeffs_s, dTCoeffs_d, cosMat, xs)

dT = dTtemp

DEALLOCATE(dTtemp)

END SUBROUTINE MultilayerT_PBC_fct_constL