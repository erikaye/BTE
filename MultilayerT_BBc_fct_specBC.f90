SUBROUTINE MultilayerT_BBC_fct_specBC(fileNum,DFT_path,write_path,matlsInds,Ls,Q0,skD,numPer, dT, kCond)
!Calculate deviational temperature and heat flux 
!through super lattice, assuming blackbody boundary conditions with dT1 and
!dT2 as temperatures at either end, and specularity spec at each interface

! assumes alternating layers of 2 (or 1) materials, but not necessarily periodic lengths

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
	
	SUBROUTINE GetRTMats(matlParams,Ls,write_path, &
							R12_d,T12_d,R21_d,T21_d,g1mCoeff_d,g2pCoeff_d, &
							R12_s,T12_s,R21_s,T21_s,g1mCoeff_s,g2pCoeff_s, &
							numG12d,numG21d,numR12d,numR21d,numT12d,numT21d, &
							numG12s,numG21s,numR12s,numR21s,numT12s,numT21s, fBoth, numFBoth, weightsMatch, spec1,spec2)
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
		REAL*8, INTENT(OUT) :: spec1(:), spec2(:), weightsMatch(:)
	END SUBROUTINE
	
	SUBROUTINE GetInterfaceBM_vSpec(matlParams,Ls,Qs,skD,bc, spec1,spec2, &
							G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
							MatABCD, MatEF, MatGH, VecJK, Amat1, f11Coeffs, f12Coeffs, f13Coeffs, Amat2, f21Coeffs, f22Coeffs, f23Coeffs)
		USE Constants
		USE omp_lib
		IMPLICIT NONE
		
		TYPE(material), DIMENSION(2), INTENT(IN) :: matlParams
					! material parameters from DFT for 2 adjacent blocks
		REAL, DIMENSION(2), INTENT(IN) :: Qs, Ls
					! heat source magnitude from laser at surface of slab, length of slabs
					! for the two slabs of interest
		CHARACTER(2), INTENT(IN) :: bc
		REAL, INTENT(IN) :: skD, spec1(:,:), spec2(:,:)
					! specularity for each mode
		TYPE(sIndVal), DIMENSION(:), INTENT(IN) :: R12_d, T12_d, T21_d, R21_d, R12_s, T12_s, T21_s, R21_s, G12_d, G12_s, G21_d, G21_s
		REAL*8, INTENT(OUT) :: MatABCD(:,:), MatEF(:,:), MatGH(:,:), VecJK(:,:), &
								Amat1(:,:), f11Coeffs(:,:), f12Coeffs(:,:), f13Coeffs(:,:), &
								Amat2(:,:), f21Coeffs(:,:), f22Coeffs(:,:), f23Coeffs(:,:)
	END SUBROUTINE
	
	SUBROUTINE CheckBCs_BBbc(fileNum, write_path, matlParams, Qs, skD, Ls, totL, numPer, &
						 fX1Coeffs, fX2Coeffs, B1wvecs_s, B1wvecs_d, P2wvecs_s, P2wvecs_d, &
						 G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
						 fBoth, numFBoth, weightsMatch, &
						 dTCoeffs_s, dTCoeffs_d, kCond)

		USE Constants
		IMPLICIT NONE
		
		CHARACTER(1), INTENT(IN) :: fileNum
		CHARACTER(longStr), INTENT(IN) :: write_path
		TYPE(material), DIMENSION(2), INTENT(IN) :: matlParams
		REAL, INTENT(IN) :: Ls(:), Qs(:), skD, totL
		REAL*8, DIMENSION(:,:,:), INTENT(IN) :: fX1Coeffs, fX2Coeffs
		REAL*8, DIMENSION(:,:,:), INTENT(IN) :: B1wvecs_s, B1wvecs_d, P2wvecs_s, P2wvecs_d
		REAL*8, DIMENSION(:,:), INTENT(IN) :: dTCoeffs_d, dTCoeffs_s
		TYPE(sIndVal), DIMENSION(:), INTENT(IN) :: R12_d, T12_d, T21_d, R21_d, R12_s, T12_s, T21_s, R21_s, G12_d, G12_s, G21_d, G21_s
		INTEGER, INTENT(IN) :: fBoth(:,:), numFBoth, numPer
		REAL*8, DIMENSION(:), INTENT(IN) :: weightsMatch
		REAL*8, DIMENSION(:), INTENT(OUT) :: kCond
	END SUBROUTINE
END INTERFACE



INTEGER, INTENT(IN) :: matlsInds(:), numPer
CHARACTER(1), INTENT(IN) :: fileNum
CHARACTER(longStr), INTENT(IN) :: DFT_path, write_path
REAL, INTENT(IN) :: Q0, skD, Ls(:)
REAL*8, DIMENSION(:) :: kCond

INTEGER(KIND=OMP_LOCK_KIND) LCK

TYPE(material), DIMENSION(numMatls) :: matlParams

INTEGER :: nL, maxNumRows, m, ierr, dummy1, dummy2
REAL, DIMENSION(:), ALLOCATABLE :: Qs
INTEGER, DIMENSION(:), ALLOCATABLE :: numW, numWd
REAL, DIMENSION(:,:), ALLOCATABLE :: spec1, spec2, weights1, weights2

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

REAL*8, DIMENSION(:,:), ALLOCATABLE :: MatABCD_t, MatEF_t, MatGH_t, VecJK_t, tempMat, tempMatFull
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: AmatX, fX1Coeffs, fX2Coeffs, fX3Coeffs, fX1Coeffs_d, fX2Coeffs_d
										
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: MatABCD1, MatEF1, MatGH1, VecJK1


REAL*8, DIMENSION(:,:), ALLOCATABLE :: dTCoeffs, dTCoeffs_s, dTCoeffs_d
REAL*8, DIMENSION(:), ALLOCATABLE :: uniqueFreqs1, uniqueFreqs2
INTEGER, DIMENSION(:), ALLOCATABLE :: numEachFreq1, numEachFreq2
REAL*8, ALLOCATABLE :: bcErr(:,:), intGs(:,:,:), intQs(:,:,:)
REAL, ALLOCATABLE :: cosMat(:,:), xs(:)
REAL*8, DIMENSION(:,:), ALLOCATABLE :: dTtemp
CHARACTER(10) :: timeT


REAL*8, DIMENSION(Nx,SIZE(matlsInds)), INTENT(OUT) :: dT


PRINT *, 'inside multilayer T fct'
PRINT *, 'T1', dT1, 'T2', dT2
PRINT *, 'Ls', Ls
PRINT *, 'number of periods', numPer
PRINT *, 'matlInds', matlsInds
PRINT *, 'rms roughness', roughness
PRINT *, 'fileNum', fileNum

CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'start time: ', timeT
nL = SIZE(matlsInds,1)*numPer  ! numLayers in unit cell * number of periods

ALLOCATE(Qs(1:nL), numW(0:nL+1), numWd(0:nL+1))
Qs(1) = Q0;

matlParams = GetMatlParams(DFT_path);

DO m = 1,nL
	Qs(m) = Q0*EXP(-SUM(Ls(1:(m-1)))/skD)
	numW(m) = matlParams(matlsInds(m))%wTot;
	numWd(m) = matlParams(matlsInds(m))%numUniqueFreq;
END DO
numW(0) = 0; numWd(0) = 0;
numW(nL+1) = 0; numWd(nL+1) = 0;

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
ALLOCATE(R12_dt(map11Max), T12_dt(MAX(map11Max,map22Max)*2), T21_dt(MAX(map11Max,map22Max)*2), &
		 R21_dt(map22Max), G12_dt(map11Max), G21_dt(map22Max), &
		 R12_st(map11Max), T12_st(MAX(map11Max,map22Max)), T21_st(MAX(map11Max,map22Max)), &
		 R21_st(map22Max), G12_st(map11Max), G21_st(map22Max))
		! *2 in T12_dt, T21_dt defn to allow for cross over of frequencies

PRINT *, 'max11', map11Max, 'max22', map22Max
		
ALLOCATE(spec1(numW(1),1),spec2(numW(2),1))
ALLOCATE(fBoth(w1d*w2d,2),weightsMatch(w1d*w2d))
! spec1(:,1) = matlParams(matlsInds(1))%spec(1:numW(1))
! spec2(:,1) = matlParams(matlsInds(nL))%spec(1:numW(nL))

WRITE(TempString, '(i4)') T
RT_fileRoot = 'RT_' // TRIM(matlParams(matlsInds(1))%matlName) // TRIM(matlParams(matlsInds(2))%matlName) &
								// '_' // TRIM(ADJUSTL(TempString)) // '_'

INQUIRE(FILE = TRIM(RT_fileRoot) // 'R12s.dat', EXIST = filesExist)
! filesExist = .FALSE.

IF (filesExist) THEN		! read the files

	CALL OMP_INIT_LOCK(LCK)	

	! !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(m,dummy1, dummy2)

	! !$OMP SINGLE
	! WRITE (*,*) 'Parallel part num threads: ', OMP_GET_NUM_THREADS()
	! !$OMP END SINGLE

	! !$OMP SECTIONS
	! !$OMP SECTION

		! CALL OMP_SET_LOCK(LCK)
		! PRINT *, 'thread num', OMP_GET_THREAD_NUM()
		! CALL OMP_UNSET_LOCK(LCK)

		OPEN(1,FILE = TRIM(RT_fileRoot) // 'R12s.dat')
		READ(1,*) dummy1, dummy2, numR12s
		DO m = 1,numR12s
			READ(1,*) R12_st(m)%row, R12_st(m)%col, R12_st(m)%indVal
		END DO
		CLOSE(1)
		
	! !$OMP SECTION

		! CALL OMP_SET_LOCK(LCK)
		! PRINT *, 'thread num', OMP_GET_THREAD_NUM()
		! CALL OMP_UNSET_LOCK(LCK)

		OPEN(2,FILE = TRIM(RT_fileRoot) // 'R12d.dat')
		READ(2,*) dummy1, dummy2, numR12d
		DO m = 1,numR12d
			READ(2,*) R12_dt(m)%row, R12_dt(m)%col, R12_dt(m)%indVal
		END DO
		CLOSE(2)
		
	! !$OMP SECTION

		! CALL OMP_SET_LOCK(LCK)
		! PRINT *, 'thread num', OMP_GET_THREAD_NUM()
		! CALL OMP_UNSET_LOCK(LCK)

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
		PRINT *, 'numG12d', numG12d
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
		PRINT *, 'numG21d', numG21d
		DO m = 1,numG21d
			READ(13,*) G21_dt(m)%row, G21_dt(m)%col, G21_dt(m)%indVal
		END DO
		CLOSE(13)

	! !$OMP SECTION
		OPEN(14,FILE = TRIM(RT_fileRoot) // 'spec1.dat')
		! DO m = 1,w1
		READ(14,*) spec1(1:w1,1)
		! END DO
		CLOSE(14)
		
	! !$OMP SECTION
		OPEN(15,FILE = TRIM(RT_fileRoot) // 'spec2.dat')
		! DO m = 1,w2
		READ(15,*) spec2(1:w2,1)
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
		
	CALL OMP_DESTROY_LOCK(LCK)
	
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





!!!!!!!!!!!!!!!!!! definition of boundary conditions !!!!!!!!!!!!

ALLOCATE(P1wvec(numW(1),1),BEndwvec(numW(nL),1))
P1wvec(:,1) = matlParams(matlsInds(1))%C(1:numW(1)) *dT1;
BEndwvec(:,1) = matlParams(matlsInds(nL))%C(1:numW(nL)) *dT2;

ALLOCATE(B1P2s((wMax+wMaxd)*2,1,0:nL),B1wvecs_s(wMax,1,0:nL),P2wvecs_s(wMax,1,0:nL), &
		 tempB1_d(wMaxd,1,0:nL), tempP2_d(wMaxd,1,0:nL), B1wvecs_d(wMax,1,0:nL), P2wvecs_d(wMax,1,0:nL), &
		 B1wvecs(wMax,1,0:nL), P2wvecs(wMax,1,0:nL))
		 ! B1wvecs_d, P2wvecs_d contain data from 1:w1, 1:w2
		 ! it's expanded from tempB1_d, tempP2_d from 1:w1d, 1:w2d

B1P2s = 0.0; B1wvecs_s = 0.0; P2wvecs_s = 0.0; tempB1_d = 0.0; tempP2_d = 0.0; 
B1wvecs_d = 0.0; P2wvecs_d = 0.0; B1wvecs = 0.0; P2wvecs = 0.0;


!!!!!! defining B0P1, BNPN+1 !!!!!!!!!


! B0P1 = B1P2s(0)
w1 = numW(1); w1d = numWd(1);  ! current values of w of interest
w2 = numW(nL); w2d = numWd(nL); 

ALLOCATE(weights1(numW(1),1),weights2(numW(nL),1),numEachFreq1(w1d), numEachFreq2(w2d))
weights1(:,1) = matlParams(matlsInds(1))%weights(1:w1);
weights2(:,1) = matlParams(matlsInds(nL))%weights(1:w2);
numEachFreq1 = matlParams(matlsInds(1))%numEachFreq(1:w1d);
numEachFreq2 = matlParams(matlsInds(nL))%numEachFreq(1:w2d);


PRINT *, 'defining B0P1'
PRINT *, 'w1', w1, 'w1d', w1d


! G12_d = summation over mu, returns same number for each element for mu for a given frequency.
! cblocksum extracts that number.

B1wvecs_s(:,:,0) = 0.0; B1wvecs_d(:,:,0) = 0.0; 
P2wvecs_s(1:w1,:,0) = P1wvec*spec1;
P2wvecs_d(1:w1,1,0) = cVECDIFF(P1wvec(:,1)*(1-spec1(:,1)),weights1(:,1),numEachFreq1);
! P2wvecs_d(1:w1,:,0) = sSPARSE_MATMUL_L(G12_d,P1wvec*(1-spec1),w1,1);
! same value for each mu for a given frequency.  value = summation over mu (G12_d * ...)

tempB1_d(:,:,0) = 0.0;
! tempP2_d(1:w1d,:,0) = cBLOCKSUM_avg(sSPARSE_MATMUL_L(G12_d,P1wvec*(1-spec1),w1,1),(/w1d,1/),numEachFreq1);
tempP2_d(1:w1d,:,0) = cBLOCKSUM_avg(P2wvecs_d(1:w1,:,0),(/w1d,1/),numEachFreq1);

CALL cVECJOIN(P2wvecs_s(1:w1,1,0), tempP2_d(1:w1d,1,0), B1P2s(1:(w1+w1d),1,0))
	! size of B1wvecs_s, B1wvecs_d is 0

B1wvecs(:,1,0) = B1wvecs_s(:,1,0) + B1wvecs_d(:,1,0)
P2wvecs(1:w1,1,0) = P2wvecs_s(1:w1,1,0) + P2wvecs_d(1:w1,1,0)
	
PRINT *, 'blocksum diffvec check', SUM(P2wvecs_d(1:w2,:,0) - sSPARSE_MATMUL_L(G12_d,P1wvec*(1-spec1),w1,1))



PRINT *, 'defining BNPN+1'
PRINT *, 'wN', w2, 'wNd', w2d

P2wvecs_s(:,:,nL) = 0.0; P2wvecs_d(:,:,nL) = 0.0; 
B1wvecs_s(1:w2,:,nL) = BEndwvec*spec2;
B1wvecs_d(1:w2,1,nL) = cVECDIFF(BEndwvec(:,1)*(1-spec2(:,1)),weights2(:,1),numEachFreq2)
! B1wvecs_d(1:w2,:,nL) = sSPARSE_MATMUL_L(G21_d,BEndwvec*(1-spec2),w2,1)

tempP2_d(:,:,nL) = 0.0;
! tempB1_d(1:w2d,:,nL) = cBLOCKSUM_avg(sSPARSE_MATMUL_L(G21_d,BEndwvec*(1-spec2),w2,1),(/w2d,1/),numEachFreq2);
tempB1_d(1:w2d,:,nL) = cBLOCKSUM_avg(B1wvecs_d(1:w2,:,nL),(/w2d,1/),numEachFreq2);

CALL cVECJOIN(B1wvecs_s(1:w2,1,nL), tempB1_d(1:w2d,1,nL), B1P2s(1:(w2+w2d),1,nL))

B1wvecs(1:w2,1,nL) = B1wvecs_s(1:w2,1,nL) + B1wvecs_d(1:w2,1,nL)
P2wvecs(:,1,nL) = P2wvecs_s(:,1,nL) + P2wvecs_d(:,1,nL)

DEALLOCATE(weights1,weights2,numEachFreq1,numEachFreq2)






!!!!!!!!!!! get MatABCD, MatEF, MatGH, VecJK !!!!!!!!!!!!!!!!!!

! maxNumRows = MAXVAL(numW+numWd) * 2;
maxNumRows = MAXVAL((numW(1:nL-1)+numWd(1:nL-1)) + (numW(2:nL)+numWd(2:nL)))
ALLOCATE(MatABCD1(maxNumRows,maxNumRows,nL-1), MatEF1(maxNumRows,MAXVAL(numW(1:nL)+numWd(1:nL)),nL-1), VecJK1(maxNumRows,1,nL-1))

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

ALLOCATE(MatABCD_t(numRows,numRows), MatEF_t(numRows,numRows1), MatGH_t(numRows,numRows0), VecJK_t(numRows,1), &
		 AmatX(N+1,N+1,nL),fX1Coeffs(N+1,wMax,nL),fX2Coeffs(N+1,wMax,nL),fX3Coeffs(N+1,1,nL))
		 
MatABCD_t = 0; MatEF_t = 0; MatGH_t = 0; VecJK_t = 0;
AmatX = 0; fX1Coeffs = 0; fX2Coeffs = 0; fx3Coeffs = 0;


! assume first layer is specularly reflecting boundary
CALL GetInterfaceBM_vSpec((/ matlParams(matlsInds(1)),matlParams(matlsInds(2)) /),Ls(1:2),Qs(1:2),skD,'sb', spec1, spec2, &
		G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
		MatABCD_t, MatEF_t, MatGH_t, VecJK_t, &
		AmatX(1:N+1,1:N+1,1), fX1Coeffs(1:N+1,1:w1,1), fX2Coeffs(1:N+1,1:w1,1), fX3Coeffs(1:N+1,:,1), &
		AmatX(1:N+1,1:N+1,2), fX1Coeffs(1:N+1,1:w2,2), fX2Coeffs(1:N+1,1:w2,2), fX3Coeffs(1:N+1,:,2))


CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'end block mats call time: ', timeT

PRINT *, 'Defining new block Mats'
PRINT *, 'size MatABCD1', numRows, SHAPE(MatABCD1)
MatABCD1(1:numRows,1:numRows,1) = MatABCD_t;
MatEF1(1:numRows,1:numRows1,1) = MatEF_t;
VecJK1(1:numRows,:,1) = VecJK_t + MATMUL(MatGH_t,B1P2s(1:(w1+w1d),:,0))  ! w1 = numFreq first layer; w0 = 0
PRINT *, 'mats layer 1 done'


DO m = 2,nL-1
	w0 = numW(m-1); w0d = numWd(m-1)
	w1 = numW(m); w1d = numWd(m)
	w2 = numW(m+1); w2d = numWd(m+1)
	w3 = numW(m+2); w3d = numWd(m+2)
	
	numRows = w1+w2+w1d+w2d;	! number of rows for given interface 
	numRows0 = w1+w1d+w0+w0d;	! number of rows for previous interface
	numRows1 = w2+w2d+w3+w3d;	! number of rows for next interface
	
	DEALLOCATE(MatABCD_t,MatEF_t,MatGH_t,VecJK_t)
	ALLOCATE(MatABCD_t(numRows,numRows), MatEF_t(numRows,w2+w2d), MatGH_t(numRows,w1+w1d), VecJK_t(numRows,1))
	
	IF (MODULO(m,2)==0) THEN  ! if even layer;  mat 2 -> mat 1
	
		CALL GetInterfaceBM_vSpec((/ matlParams(matlsInds(m)),matlParams(matlsInds(m+1)) /),Ls(m:m+1),Qs(m:m+1),skD,'bb', spec1, spec2, &
				G21_s, R21_s, T21_s, G12_s, R12_s, T12_s, G21_d, R21_d, T21_d, G12_d, R12_d, T12_d, &
				MatABCD_t, MatEF_t, MatGH_t, VecJK_t, &
				AmatX(1:N+1,1:N+1,m), fX1Coeffs(1:N+1,1:w1,m), fX2Coeffs(1:N+1,1:w1,m), fX3Coeffs(1:N+1,:,m), &
				AmatX(1:N+1,1:N+1,m+1), fX1Coeffs(1:N+1,1:w2,m+1), fX2Coeffs(1:N+1,1:w2,m+1), fX3Coeffs(1:N+1,:,m+1))
				
	ELSE  ! if odd layer;  mat 1 -> mat 2
	
		CALL GetInterfaceBM_vSpec((/ matlParams(matlsInds(m)),matlParams(matlsInds(m+1)) /),Ls(m:m+1),Qs(m:m+1),skD,'bb', spec1, spec2, &
				G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
				MatABCD_t, MatEF_t, MatGH_t, VecJK_t, &
				AmatX(1:N+1,1:N+1,m), fX1Coeffs(1:N+1,1:w1,m), fX2Coeffs(1:N+1,1:w1,m), fX3Coeffs(1:N+1,:,m), &
				AmatX(1:N+1,1:N+1,m+1), fX1Coeffs(1:N+1,1:w2,m+1), fX2Coeffs(1:N+1,1:w2,m+1), fX3Coeffs(1:N+1,:,m+1))
				
	END IF
	
	
	ALLOCATE(tempMat(numRows,w1+w1d),tempMatFull(numRows,numRows))
	tempMat = 0.0; tempMatFull = 0.0;
	tempMat = cBLOCK_MATMUL(MatGH_t, & 
						cAXB(MatABCD1(1:numRows0,1:numRows0,m-1),MatEF1(1:numRows0,1:w1+w1d,m-1)),(/0,w1,0,w1d/),(/w0,w1,w0d,w1d/));
			! full matrix has columns of 0 where MatEF has columns of 0
	tempMatFull(1:numRows,1:w1) = tempMat(1:numRows,1:w1)
	tempMatFull(1:numRows,w1+w2+1:w1+w2+w1d) = tempMat(1:numRows,w1+1:w1+w1d)
	
	MatABCD1(1:numRows,1:numRows,m) = MatABCD_t - tempMatFull;
	! MatABCD1(1:numRows,1:numRows,m) = MatABCD_t - cBLOCK_MATMUL(MatGH_t, & 
						! cAXB(MatABCD1(1:numRows0,1:numRows0,m-1),MatEF1(1:numRows0,1:w1+w1d,m-1)),(/0,w1,0,w1d/),(/w0,w1,w0d,w1d/));
	MatEF1(1:numRows,1:w2+w2d,m) = MatEF_t;
	VecJK1(1:numRows,:,m) = VecJK_t + cBLOCK_MATMUL(MatGH_t, &
						cAXB(MatABCD1(1:numRows0,1:numRows0,m-1),VecJK1(1:numRows0,:,m-1)),(/0,w1,0,w1d/),(/w0,w1,w0d,w1d/));
	
	DEALLOCATE(tempMat,tempMatFull)
	
	PRINT *, 'layer ', m, ' done'
	
END DO
	
DEALLOCATE(MatABCD_t,MatEF_t,MatGH_t,VecJK_t)

! solve explicitly for B, P vecs

PRINT *, 'solving for b1p2'
CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'start time cAXB MatABCD: ', timeT

OPEN(1,FILE = TRIM(write_path) // 'output_matErr_' // fileNum // '.dat')

PRINT *, 'updated b1p2 calc'
DO m = nL-1,1,-1
	w0 = numW(m-1); w0d = numWd(m-1);
	w1 = numW(m); w1d = numWd(m); 
	w2 = numW(m+1); w2d = numWd(m+1);
	w3 = numW(m+2); w3d = numWd(m+2);
		
	numRows0 = w0+w0d+w1+w1d;	
	numRows = w1+w1d+w2+w2d;
	numRows1 = w2+w2d+w3+w3d;  ! w3, w3d = 0 for a 2 layer interface

	B1P2s(1:numRows,:,m) = cAXB(MatABCD1(1:numRows,1:numRows,m), &  ! square matrix
									cBLOCK_MATMUL(MatEF1(1:numRows,1:w2+w2d,m),B1P2s(1:numRows1,:,m+1),(/0,w2,0,w2d/),(/0,w2,w3,w2d/)) & 
													+ VecJK1(1:numRows,:,m))
	! ! B1P2s(1:numRows,:,1) = cAXB(MatABCD1(1:numRows,1:numRows,1), &  ! square matrix
					! ! MATMUL(MatEF1(1:numRows,1:w2+w2d,1),B1P2s(1:(w2+w2d),:,nL)) + VecJK1(1:numRows))

	! ! specular part
	! B1P2s(1:(w1+w2),:,1) = cAXB(MatABCD1(1:(w1+w2),1:(w1+w2),m), &  ! square matrix
					! MATMUL(MatEF1(1:(w1+w2),1:w2,m),B1P2s(1:w2,:,m+1)) + VecJK1(1:(w1+w2),:,m))
					
	! ! diffuse part
	! B1P2s((w1+w2+1):(numRows),:,1) = cAXB(MatABCD1((w1+w2+1):(numRows),(w1+w2+1):(numRows),m), &  ! square matrix
					! MATMUL(MatEF1((w1+w2+1):(numRows),(w2+1):(w2+w2d),m),B1P2s((w2+1):(w2+w2d),:,m+1)) + VecJK1((w1+w2+1):(numRows),:,m))

					
	PRINT *, 'finished solving for b1p2'
	CALL DATE_AND_TIME(TIME=timeT)
	PRINT *, 'end time cAXB MatABCD: ', timeT
	 
	B1wvecs_s(1:w1,1,m) = B1P2s(1:w1,1,m); 
	PRINT *, 'B1wvecs_s done'
	P2wvecs_s(1:w2,1,m) = B1P2s((w1+1):(w1+w2),1,m);
	PRINT *, 'P2wvecs_s done'
	tempB1_d(1:w1d,1,m) = B1P2s((w1+w2+1):(w1+w2+w1d),1,m);
	PRINT *, 'tempB1_d done'
	tempP2_d(1:w2d,1,m) = B1P2s((w1+w2+w1d+1):(w1+w1d+w2+w2d),1,m);
	PRINT *, 'tempP2_d done'

	B1wvecs_d(1:w1,1,m) = cDIFFVECtoFULL(tempB1_d(1:w1d,1,m), w1, matlParams(matlsInds(m))%numEachFreq(1:w1d));
	P2wvecs_d(1:w2,1,m) = cDIFFVECtoFULL(tempP2_d(1:w2d,1,m), w2, matlParams(matlsInds(m+1))%numEachFreq(1:w2d));
	
	B1wvecs(1:w1,1,m) = B1wvecs_s(1:w1,1,m) + B1wvecs_d(1:w1,1,m)
	P2wvecs(1:w2,1,m) = P2wvecs_s(1:w2,1,m) + P2wvecs_d(1:w2,1,m)
	
	
	ALLOCATE(matErr(numRows,1))
	! matErr = MATMUL(MatABCD1,B1P2s(1:numRows,:,1)) - MATMUL(MatEF1,B1P2s(1:(w2+w2d),:,2)) - VecJK1
	matErr = MATMUL(MatABCD1(1:numRows,1:numRows,m),B1P2s(1:numRows,:,m)) - &
				cBLOCK_MATMUL(MatEF1(1:numRows,1:(w2+w2d),m),B1P2s(1:numRows,:,m+1),(/0,w2,0,w2d/),(/0,w2,w3,w2d/)) - VecJK1(1:numRows,:,m)

	! append to existing file?
	! IF (m == 1) THEN
		! ! create a new file
		! OPEN(1,FILE = TRIM(write_path) // 'output_matErr_' // fileNum // '.dat')
		! WRITE(1,"(40000(ES15.7,','))") REAL(matErr(:,1))
		! CLOSE(1)
	! ELSE
		! ! append to existing file
		! PRINT *, 'appended matErr?'
		! OPEN(1,FILE = TRIM(write_path) // 'output_matErr_' // fileNum // '.dat', POSITION='append')
		! WRITE(1,"(40000(ES15.7,','))") REAL(matErr(:,1))
		! CLOSE(1)
	! END IF

	WRITE(1,"(40000(ES15.7,','))") REAL(matErr(:,1))
	
	DEALLOCATE(matErr)

END DO
	
CLOSE(1)



!!!!!!!!!!!!!!! Calculating dT !!!!!!!!!!!!!

PRINT *, 'Calculating dT'
CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'start time calc dT: ', timeT
	
! Calculation of deviational temperature distribution
ALLOCATE(dTCoeffs(N+1,nL), dTCoeffs_s(N+1,nL), dTCoeffs_d(N+1,nL), cosMat(Nx,N+1), xs(Nx))

dTCoeffs = 0;
DO m = 1,nL
	w1 = numW(m)
	w1d = numWd(m)
	dTCoeffs(:,m) = RESHAPE(cAXB(AmatX(:,:,m), &
						(MATMUL(fX1Coeffs(:,1:w1,m), P2wvecs(1:w1,:,m-1)) + &
						 MATMUL(fX2Coeffs(:,1:w1,m), B1wvecs(1:w1,:,m)) + &
						fX3Coeffs(:,:,m))),(/w1/))
	dTCoeffs_s(:,m) = RESHAPE(cAXB(AmatX(:,:,m), &
						(MATMUL(fX1Coeffs(:,1:w1,m), P2wvecs_s(1:w1,:,m-1)) + &
						 MATMUL(fX2Coeffs(:,1:w1,m), B1wvecs_s(1:w1,:,m)) + &
						fX3Coeffs(:,:,m))),(/w1/))
	dTCoeffs_d(:,m) = RESHAPE(cAXB(AmatX(:,:,m), &
						(MATMUL(fX1Coeffs(:,1:w1,m), P2wvecs_d(1:w1,:,m-1)) + &
						 MATMUL(fX2Coeffs(:,1:w1,m), B1wvecs_d(1:w1,:,m)) + &
						fX3Coeffs(:,:,m))),(/w1/))
END DO 




!!!!!!!!!!!!!!!! calculating dT !!!!!!!!!!!!!!!!!

xs = (/ (m*(1.0/(Nx-1)), m = 0,Nx-1) /)

DO m = 0,N
	cosMat(:,m+1) = COS(m*pi*xs)
END DO

ALLOCATE( dTtemp(SIZE(dT,1), SIZE(dT,2)) )
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

PRINT *, 'size fboth', SHAPE(fBoth)

OPEN(2,FILE = TRIM(write_path) // 'output_bcErrS_' // fileNum // '.dat')
OPEN(3,FILE = TRIM(write_path) // 'output_bcErrD_' // fileNum // '.dat')
OPEN(4,FILE = TRIM(write_path) // 'output_bcErr_' // fileNum // '.dat')

w1 = numW(1); w2 = numW(2)
PRINT *, 'w1', w1, 'w2', w2

! CALL CheckBCs_BBbc( fileNum, write_path, (/ matlParams(matlsInds(1)),matlParams(matlsInds(2)) /), Qs, skD, Ls, numPer, &
					 ! fX2Coeffs(:,1:w1,1), fX1Coeffs(:,1:w2,2), &
					 ! B1wvecs_s(:,:,0:nL), B1wvecs_d(:,:,0:nL), P2wvecs_s(:,:,0:nL), P2wvecs_d(:,:,0:nL), &
					 ! G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
					 ! fBoth, numFBoth, weightsMatch, &
					 ! dTCoeffs_s, dTCoeffs_d, kCond)

DO m = 1,nL-1
						 
	IF (MODULO(m,2)==1) THEN	! mat 1 -> mat 2
		CALL CheckBCs_BBbc( fileNum, write_path, (/ matlParams(matlsInds(m)),matlParams(matlsInds(m+1)) /), &
							 Qs(m:m+1), skD, Ls(m:m+1), SUM(Ls), numPer, &
							 fX1Coeffs(:,:,m:m+1), fX2Coeffs(:,:,m:m+1), &
							 B1wvecs_s(:,:,m-1:m+1), B1wvecs_d(:,:,m-1:m+1), P2wvecs_s(:,:,m-1:m+1), P2wvecs_d(:,:,m-1:m+1), &
							 G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
							 fBoth(:,1:2), numFBoth, weightsMatch, &
							 dTCoeffs_s(:,m:m+1), dTCoeffs_d(:,m:m+1), kCond(m:m))
	ELSE
		CALL CheckBCs_BBbc( fileNum, write_path, (/ matlParams(matlsInds(m)),matlParams(matlsInds(m+1)) /), &
							 Qs(m:m+1), skD, Ls(m:m+1), SUM(Ls), numPer, &
							 fX1Coeffs(:,:,m:m+1), fX2Coeffs(:,:,m:m+1), &
							 B1wvecs_s(:,:,m-1:m+1), B1wvecs_d(:,:,m-1:m+1), P2wvecs_s(:,:,m-1:m+1), P2wvecs_d(:,:,m-1:m+1), &
							 G21_s, R21_s, T21_s, G12_s, R12_s, T12_s, G21_d, R21_d, T21_d, G12_d, R12_d, T12_d, &
							 fBoth(:,2:1:-1), numFBoth, weightsMatch, &
							 dTCoeffs_s(:,m:m+1), dTCoeffs_d(:,m:m+1), kCond(m:m))	
	END IF
END DO

CLOSE(2)
CLOSE(3)
CLOSE(4)

PRINT *, 'interface conductance', kCond			 

CALL DATE_AND_TIME(TIME=timeT)
PRINT *, 'end time checkbcs: ', timeT

DEALLOCATE(spec1,spec2, numW, numWd, Qs)
DEALLOCATE(fBoth, weightsMatch)
DEALLOCATE(P1wvec, BEndwvec, B1P2s, B1wvecs, P2wvecs, tempB1_d, tempP2_d, &
		   B1wvecs_s, B1wvecs_d, P2wvecs_s, P2wvecs_d)
DEALLOCATE(AmatX, fX1Coeffs, fX2Coeffs, fX3Coeffs, MatABCD1, MatEF1, VecJK1, dTCoeffs, dTCoeffs_s, dTCoeffs_d, cosMat, xs)

dT = dTtemp

DEALLOCATE(dTtemp)
PRINT *, 'leaving multilayerT'

END SUBROUTINE MultilayerT_BBC_fct_specBC