MODULE constants

USE IFPORT

IMPLICIT NONE

! CHARACTER(1), PARAMETER :: fileNum = '1';
REAL, PARAMETER :: roughness = 0.8e-10
! REAL, PARAMETER :: roughness = 0.2e-10

REAL*8, PARAMETER :: pi = 3.1415927, kB = 1.3806488e-23, h = 6.62606957e-34
REAL*8, PARAMETER :: hbar = h/2/pi
INTEGER, PARAMETER :: T = 200       ! ambient temperature
INTEGER, PARAMETER :: maxLen = 16384		! max side length of a square matrix

! parameters for extracting material parameters from DFT calculations
INTEGER,PARAMETER :: maxPts = 32769, maxModes = 6, maxLines = MIN(maxPts*maxModes,200000)
						! max kpts, max num modes, max lines in file
INTEGER,PARAMETER :: shortStr = 4, medStr = 32, longStr = 128

! parameters for matching k's, freqs
REAL, PARAMETER :: ktol = 1e-7, fgs_tol = 0.01
! tolerance (fraction) for error when matching k's, freqs, freq comparision with gaussian smearing

! parameters describing the superlattice stack
INTEGER, PARAMETER :: numMatls = 2  ! number of unique materials
! INTEGER, PARAMETER :: numLayers = 2 ! number of layers in stack
CHARACTER(shortStr), DIMENSION(numMatls) :: matlist = (/'Si','SiX'/)   ! list of unique materials
! INTEGER, DIMENSION (numLayers):: matStack = (/1,1/)  ! numbers correspond to position in matlist
! REAL, DIMENSION(numLayers) :: LStack = (/1E-7,1E-7/)    ! thickness of each layer in stack
! REAL, PARAMETER :: Q0 = 0 ! amplitude of heat source at surface

! parameters for cosine expansion coefficients
INTEGER, PARAMETER :: N = 40        ! number of elements in cos expansion
REAL, PARAMETER :: eta = 0, kr = 0  ! fourier transform component in time, radial hankel transform component
INTEGER, PARAMETER :: Nx = 1000		! number of xpoints in final dT profile

! parameters for boundary conditions
REAL, PARAMETER :: dT1 = 1, dT2 = -1   ! temperatures at edges of SL

TYPE, PUBLIC :: material
    CHARACTER(shortStr) :: matlName ! material name
    CHARACTER(longStr) :: fileLoc   ! file location
    ! INTEGER :: numPts, numModes     ! number of k points, number of modes
	INTEGER :: numModes, wTot					! total number of modes
    REAL, DIMENSION(maxPts*maxModes) :: kx, ky, kz
	REAL, DIMENSION(maxPts*maxModes) :: dkx, dky, dkz
        ! relevant dimension:  length = numPts
    REAL, DIMENSION(maxPts*maxModes) :: vgx, vgy, vgz, tau, freq, C, MFP, weights, spec
        ! relevant dimensions:  numPts,numModes
		! weights sum over 1/2 Brillouin Zone
    REAL*8, DIMENSION(maxPts*maxModes) :: gamma
        ! relevant dimensions:  numPts,numModes
		
	! ordering by frequency information
	INTEGER :: numUniqueFreq
	INTEGER, DIMENSION(maxPts*maxModes) :: numEachFreq
	REAL, DIMENSION(maxPts*maxModes) :: uniqueFreqs
	REAL :: df, rho
END TYPE material


TYPE, PUBLIC :: sIndVal
	! for sparse matrices, list elements as (row,col,val)
	INTEGER :: row, col
	REAL :: indVal
END TYPE sIndVal



CONTAINS
	FUNCTION numLIMVAL(inNum, nanVal, infVal) RESULT(outNum)
		! gets rid of NaN's (set to nanVal), and Infs (set to infVal)
		REAL*8, INTENT(IN) :: inNum
		REAL*8, INTENT(IN) :: nanVal, infVal
		REAL*8 :: outNum

		outNum = inNum;
		
		! replace NaN's and Infs
		IF(ISNAN(inNum)) outNum = nanVal;
		IF(ABS(inNum) > 1.0E60) outNum = infVal;
		! IF(inNum > 1.0E60) outNum = infVal;
		! IF(inNum < -1.0E60) outNum = -1.0*infVal;
			! could be #/small neg number, which is just from a numerical error
		
		! also replace values that should be approx 0, 1, but are not due to manually introduced small errors
		IF(ABS(inNum - 0.0) < 1.0E-60) outNum = nanVal
		IF(ABS(inNum - 1.0) < 1.0E-60) outNum = infVal
		IF(ABS(inNum - (-1.0)) < 1.0E-60) outNum = -1.0*infVal	
	END FUNCTION

	FUNCTION vecLIMVAL(inVec, nanVal, infVal) RESULT(outVec)
		! gets rid of NaN's (set to nanVal), and Infs (set to infVal)
		REAL*8, DIMENSION(:), INTENT(IN) :: inVec
		REAL*8, INTENT(IN) :: nanVal, infVal
		REAL*8, DIMENSION(SIZE(inVec)) :: outVec

		outVec = inVec;
		
		! replace NaN's and Infs
		WHERE(ISNAN(inVec)) outVec = nanVal;
		WHERE(ABS(inVec) > 1.0E60) outVec = infVal;
		! WHERE(inVec > 1.0E60) outVec = infVal;
		! WHERE(inVec < -1.0E60) outVec = -1.0*infVal;
		
		! also replace values that should be approx nanVal, infVal, but are not due to manually introduced small errors
		WHERE(ABS(inVec - nanVal) < 1.0E-60) outVec = nanVal
		WHERE(ABS(infVal - inVec) < 1.0E-60) outVec = infVal
		WHERE(ABS(inVec - (-1.0)*infVal) < 1.0E-60) outVec = -1.0*infVal
		
	END FUNCTION

	SUBROUTINE cVECJOIN(rVec, lVec, totVec)
	! connect two vectors
		REAL*8, DIMENSION(:), INTENT(IN) :: rVec, lVec
		REAL*8, DIMENSION(:), INTENT(OUT) :: totVec
		INTEGER :: rLen, lLen, totLen
		
		rLen = SIZE(rVec); lLen = SIZE(lVec); totLen = SIZE(totVec)
		IF (rLen + lLen == totLen) THEN
			totVec(1:rLen) = rVec
			totVec((rLen+1):(rLen+lLen)) = lVec
		ELSE
			STOP ' vecjoin: dimension mismatch'
		END IF
	END SUBROUTINE cVECJOIN

	FUNCTION cDIAG(vec) RESULT(mat)
		REAL*8, DIMENSION(:), INTENT(IN) :: vec
		REAL*8, DIMENSION(SIZE(vec,1),SIZE(vec,1)) :: mat
		REAL*8, DIMENSION(SIZE(vec,1)**2) :: temp
		INTEGER :: vLen
		
		temp = 0;
		mat = 0;
		
		vLen = SIZE(vec,1)
		temp(1:vLen**2:(vLen+1)) = vec;
		mat = RESHAPE(temp,(/vLen,vLen/))
	END FUNCTION cDIAG

	FUNCTION cBLOCK_MATMUL(matInL,matInR,numElsL,numElsR) RESULT(matOut)
		REAL*8, INTENT(IN) :: matInL(:,:), matInR(:,:)
		INTEGER, INTENT(IN) :: numElsL(:), numElsR(:)  
			! (#cols skipped, #col in first block, #cols skipped, #col in next block, etc)
			! must be even in length
		REAL*8, DIMENSION(SIZE(matInL,1),SIZE(matInR,2)) :: matOut
		INTEGER :: m, ind1L, ind2L, ind1R, ind2R
		! matInL has columnar blocks, or matInR has blocks along each row
		matOut = 0;
		
		ind1L = 0; ind2L = 0; ind1R = 0; ind2R = 0;
		IF (SIZE(numElsR) .NEQV. SIZE(numElsL)) THEN
			 PRINT *, 'error in numElsDef'
		END IF
		DO m = 1,SIZE(numElsL),2
			ind1L = SUM(numElsL(1:m)) + 1
			ind2L = SUM(numElsL(1:m+1))
			ind1R = SUM(numElsR(1:m)) + 1
			ind2R = SUM(numElsR(1:m+1))
			matOut = matOut + MATMUL(matInL(:,ind1L:ind2L),matInR(ind1R:ind2R,:))
		END DO
		
	END FUNCTION
	
	FUNCTION sDIAG_MATMUL_L(diagVec, matIn) RESULT(matOut)
		REAL*8, INTENT(IN) :: diagVec(:), matIn(:,:)
		REAL*8 :: matOut(SIZE(matIn,1),SIZE(matIn,2))
		INTEGER :: i, nRows
		
		matOut = 0;
		
		nRows = SIZE(matIn,1) ! = SIZE(diagVec)
		
		IF (SIZE(diagVec)/=SIZE(matIn,1)) THEN
			PRINT *, 'error in diag_matmul_l'
		END IF
		
		DO i = 1,nRows
			matOut(i,:) = diagVec(i)*matIn(i,:)
		END DO	
	END FUNCTION
	
	SUBROUTINE sDIAG_MATMUL_L_SUB(diagVec, matIn)
		REAL*8, INTENT(IN) :: diagVec(:)
		REAL*8 :: matIn(:,:)
		REAL*8, ALLOCATABLE :: matOut(:,:)
		INTEGER :: i, nRows
		
		ALLOCATE(matOut(SIZE(matIn,1),SIZE(matIn,2)))
		matOut = 0;
		
		nRows = SIZE(matIn,1) ! = SIZE(diagVec)
		
		IF (SIZE(diagVec)/=SIZE(matIn,1)) THEN
			PRINT *, 'error in diag_matmul_l'
		END IF
		
		DO i = 1,nRows
			matOut(i,:) = diagVec(i)*matIn(i,:)
		END DO	
		
		matIn = matOut;
		DEALLOCATE(matOut)
	END SUBROUTINE
	
	SUBROUTINE sDIAG_MATMUL_R_SUB(matIn, diagVec)
		REAL*8, INTENT(IN) :: diagVec(:)
		REAL*8 :: matIn(:,:)
		REAL*8, ALLOCATABLE :: matOut(:,:)
		INTEGER :: i, nCols
		
		ALLOCATE(matOut(SIZE(matIn,1),SIZE(matIn,2)))
		matOut = 0;
		
		nCols = SIZE(matIn,2) ! = SIZE(diagVec)
		
		IF (SIZE(diagVec)/=SIZE(matIn,2)) THEN
			PRINT *, 'error in diag_matmul_r'
		END IF
		
		DO i = 1,nCols
			matOut(:,i) = matIn(:,i)*diagVec(i);
		END DO	
		
		matIn = matOut;
		DEALLOCATE(matOut)
	END SUBROUTINE
	
	FUNCTION sDIAG_MATMUL_R(matIn, diagVec) RESULT(matOut)
		REAL, INTENT(IN) :: diagVec(:)
		REAL*8, INTENT(IN) :: matIn(:,:)
		REAL*8 :: matOut(SIZE(matIn,1),SIZE(matIn,2))
		INTEGER :: i, nCols
		
		! PRINT *, 'in diagMatmulR'
		
		matOut = 0;
		
		nCols = SIZE(matIn,2) ! = SIZE(diagVec)
		! PRINT *, nCols
		
		DO i = 1,nCols
			matOut(:,i) = diagVec(i)*matIn(:,i)
		END DO	
		
		! PRINT *, 'end diagMatmulR'
	END FUNCTION

	FUNCTION cBLOCKSUM(block, outSize, numEachFreq) RESULT(blockOut)
	! sum together terms of the same frequency
	! R/T mats (with weighting per omega) is usually included
		REAL*8, DIMENSION(:,:), INTENT(IN) :: block
		INTEGER, DIMENSION(2), INTENT(IN) :: outSize  ! (/numrows, numcols/)
		INTEGER, DIMENSION(:), INTENT(IN) :: numEachFreq
				! number of modes for certain frequency
		REAL*8, DIMENSION(outSize(1),outSize(2)) :: blockOut
		INTEGER :: i, j, nFreqs, order, indNum
		! REAL*8, DIMENSION(:,:), ALLOCATABLE :: temp
		
		nFreqs = SIZE(numEachFreq)
		blockOut = 0;
		
		IF ((SIZE(block,1) == outSize(1)).AND.(SIZE(block,2) == outSize(2))) THEN
			! no summation required
			blockOut = block;
			
		ELSE
			! summation to reduce block size
			IF (SIZE(block,1) == outSize(1)) THEN
				order = 1;
			ELSEIF (SIZE(block,2) == outSize(2)) THEN
				order = 2;
			ELSE
				STOP 'error in blockSum'
			END IF
			
			! PRINT *, 'start blockSum', order
			! PRINT *, 'nFreqs', nFreqs
			! PRINT *, 'w', SIZE(block,1)
			
			IF (order == 1) THEN   ! sum across the columns
				! PRINT *, SIZE(numEachFreq), SUM(numEachFreq), SIZE(block,2)
				DO i = 1,nFreqs  ! number of unique frequencies
					indNum = SUM(numEachFreq(1:i-1));
					blockOut(:,i) = SUM(block(:,(indNum+1):(indNum+numEachFreq(i))),2);
				END DO
			ELSE
				! sum across the rows
				DO i = 1,nFreqs  ! number of unique frequencies
					indNum = SUM(numEachFreq(1:i-1));
					blockOut(i,:) = SUM(block((indNum+1):(indNum+numEachFreq(i)),:),1);
				END DO
			END IF
		END IF
		! PRINT *, 'end blockSum'
		
	END FUNCTION cBLOCKSUM

	FUNCTION cBLOCKSUM_avg(block, outSize, numEachFreq) RESULT(blockOut)
	! old:  sum together terms of the same frequency
	! now:  find the average of terms of the same frequency
		REAL*8, DIMENSION(:,:), INTENT(IN) :: block
		INTEGER, DIMENSION(2), INTENT(IN) :: outSize  ! (/numrows, numcols/)
		INTEGER, DIMENSION(:), INTENT(IN) :: numEachFreq
				! number of modes for certain frequency
		REAL*8, DIMENSION(outSize(1),outSize(2)) :: blockOut
		INTEGER :: i, j, nFreqs, order, indNum
		
		nFreqs = SIZE(numEachFreq)
		blockOut = 0;
		
		IF ((SIZE(block,1) == outSize(1)).AND.(SIZE(block,2) == outSize(2))) THEN
			! no summation required
			blockOut = block;
			
		ELSE
			! summation to reduce block size
			IF (SIZE(block,1) == outSize(1)) THEN
				order = 1;
			ELSEIF (SIZE(block,2) == outSize(2)) THEN
				order = 2;
			ELSE
				STOP 'error in blockSum'
			END IF
			
			IF (order == 1) THEN   ! sum across the columns
				DO i = 1,nFreqs  ! number of unique frequencies
					indNum = SUM(numEachFreq(1:i-1));
					blockOut(1,i) = SUM(block(1,(indNum+1):(indNum+numEachFreq(i))))/numEachFreq(i)
				END DO
			ELSE
				! sum across the rows
				DO i = 1,nFreqs  ! number of unique frequencies
					indNum = SUM(numEachFreq(1:i-1));
					blockOut(i,1) = SUM(block((indNum+1):(indNum+numEachFreq(i)),1))/numEachFreq(i)
				END DO
			END IF
		END IF
		! PRINT *, 'end blockSum'
		
	END FUNCTION cBLOCKSUM_avg
	
	! ! delete?
	! FUNCTION sVECDIFF(vec, weights, numEachFreq) RESULT(vecOut)
	! ! make vector s.t elements corresponding to modes with the same frequency have the same
	! ! value (value is the mean (of vec) over angles for given frequency)
		! REAL, DIMENSION(:), INTENT(IN) :: vec, weights  ! will be 1-dimensional
		! INTEGER, DIMENSION(:), INTENT(IN) :: numEachFreq  ! same size as vec
		! INTEGER :: nFreqs, i, indNum
		
		! REAL, ALLOCATABLE :: vecSum(:)
		! REAL :: vecOut(SIZE(vec))

		! nFreqs = SIZE(numEachFreq)
		
		! ALLOCATE(vecSum(nFreqs))
		
		! vecSum = 0;
		! DO i = 1,nFreqs
			! indNum = SUM(numEachFreq(1:i-1));
			! vecSum(i) = SUM(vec((indNum+1):(indNum+numEachFreq(i))) * weights((indNum+1):(indNum+numEachFreq(i)))) / &
								! SUM(weights((indNum+1):(indNum+numEachFreq(i))))
		! END DO
		
		! vecOut = 0;
		! DO i = 1,nFreqs
			! indNum = SUM(numEachFreq(1:i-1));
			! vecOut((indNum+1):(indNum+numEachFreq(i))) = vecSum(i) 
		! END DO
		
		! DEALLOCATE(vecSum)
		
	! END FUNCTION sVECDIFF

	FUNCTION cVECDIFF(vec, weights, numEachFreq) RESULT(vecOut)
	! make vector s.t elements corresponding to modes with the same frequency have the same
	! value (value is the mean (of vec) over angles for given frequency)

		REAL*8, DIMENSION(:), INTENT(IN) :: vec, weights  ! will be 1D
		INTEGER, DIMENSION(:), INTENT(IN) :: numEachFreq
		INTEGER :: nFreqs, i, indNum
		
		REAL*8 :: vecSum
		REAL*8 :: vecOut(SIZE(vec))
		
		nFreqs = SIZE(numEachFreq)
		
		
		vecSum = 0;
		DO i = 1,nFreqs  ! number of masks
			indNum = SUM(numEachFreq(1:i-1));
			! vecSum(i) = SUM(vec((indNum+1):(indNum+numEachFreq(i))))/numEachFreq(i)
			vecSum = SUM(vec((indNum+1):(indNum+numEachFreq(i)))*weights((indNum+1):(indNum+numEachFreq(i)))) / &
																		SUM(weights((indNum+1):(indNum+numEachFreq(i))))
							! weighting to match defn in R/T matrices (RT normalized by weights/sum(weights))
			! mean or sum? MEAN bc integrating over all angles, but spreading out again
			
			! want this to be weighted average of values in vec. which it should be.
			
			vecOut((indNum+1):(indNum+numEachFreq(i))) = vecSum
		END DO
		
	END FUNCTION cVECDIFF
	
	FUNCTION cVECAVG(vec, numEachFreq) RESULT(vecOut)
	! make vector s.t elements corresponding to modes with the same frequency have the same
	! value (value is the avg (no weights) over angles for given frequency)

		REAL*8, DIMENSION(:), INTENT(IN) :: vec ! will be 1D
		INTEGER, DIMENSION(:), INTENT(IN) :: numEachFreq
		INTEGER :: nFreqs, i, indNum
		
		REAL*8, ALLOCATABLE :: vecSum(:)
		REAL*8 :: vecOut(SIZE(vec))
		
		nFreqs = SIZE(numEachFreq)
		
		! ALLOCATE(vecSum(nFreqs))
		
		! vecSum = 0;
		! DO i = 1,nFreqs  ! number of masks
			! indNum = SUM(numEachFreq(1:i-1));
			! vecSum(i) = SUM(vec((indNum+1):(indNum+numEachFreq(i))))/numEachFreq(i)
		! END DO
		
		vecOut = 0;
		DO i = 1,nFreqs  ! number of masks
			indNum = SUM(numEachFreq(1:(i-1)));
			! vecOut((indNum+1):(indNum+numEachFreq(i))) = vecSum(i)
			vecOut((indNum+1):(indNum+numEachFreq(i))) = SUM(vec((indNum+1):(indNum+numEachFreq(i))))/numEachFreq(i)
		END DO
		
		! vecOut = vec;
		
		! DEALLOCATE(vecSum)
	END FUNCTION cVECAVG

	FUNCTION cVECNORM(vec, numEachFreq) RESULT(vecOut)
	! make vector s.t elements corresponding to modes with the same frequency sum to 1 but have the same relative weighting

		REAL*8, DIMENSION(:), INTENT(IN) :: vec ! will be 1D
		INTEGER, DIMENSION(:), INTENT(IN) :: numEachFreq
		INTEGER :: nFreqs, i, indNum
		
		REAL*8, ALLOCATABLE :: vecSum(:)
		REAL*8 :: vecOut(SIZE(vec))
		
		nFreqs = SIZE(numEachFreq)
		
		ALLOCATE(vecSum(nFreqs))
		
		vecSum = 0;
		! DO i = 1,nFreqs  ! number of masks
			! indNum = SUM(numEachFreq(1:i-1));
			! vecSum(i) = SUM(vec((indNum+1):(indNum+numEachFreq(i))))/numEachFreq(i)
		! END DO
		
		vecOut = 0;
		DO i = 1,nFreqs  ! number of masks
			indNum = SUM(numEachFreq(1:(i-1)));
			vecOut((indNum+1):(indNum+numEachFreq(i))) = vec((indNum+1):(indNum+numEachFreq(i)))/SUM(vec((indNum+1):(indNum+numEachFreq(i))))
		END DO
		
		DEALLOCATE(vecSum)
	END FUNCTION cVECNORM
	
	FUNCTION cDIFFVECtoFULL(vecIn, outLen, numEachFreq) RESULT (vecOut)
		! should be inverse of cBLOCKSUM_avg
		! takes a vec (of length wxd) and repeats it so that it's of length wx
		
		! now inverse of blocksum along col (summing rows, which are not weighted)
		REAL*8, INTENT(IN) :: vecIn(:)
		INTEGER, INTENT(IN) :: outLen
		INTEGER, INTENT(IN) :: numEachFreq(:)
		REAL*8 :: vecOut(outLen)
		INTEGER :: i, indNum
		
		vecOut = 0;
		! PRINT *, SIZE(vecIn)
		! PRINT *, SUM(numEachFreq)
		! PRINT *, SIZE(weights)
		
		DO i = 1, SIZE(vecIn)
			indNum = SUM(numEachFreq(1:(i-1)))
			! IF (i==1) THEN
				! PRINT *, indNum
			! END IF
			vecOut((indNum+1):indNum+numEachFreq(i)) = vecIn(i)!/numEachFreq(i)
			! vecOut((indNum+1):indNum+numEachFreq(i)) = vecIn(i)/SUM(weights(indNum+1:indNum+numEachFreq(i)))
		END DO 
	END FUNCTION 
	
	FUNCTION sSPARSE_MATMUL_L(spList, matIn, nRows, nCols) RESULT(matOut)
		TYPE(sIndVal), INTENT(IN) :: spList(:)
		REAL, INTENT(IN) :: matIn(:,:)
		INTEGER, iNTENT(IN) :: nRows, nCols
		REAL :: matOut(nRows,nCols)
		INTEGER :: i,j,m
		
		matOut = 0;
		
		! sparse mat * normal mat
		DO m = 1,SIZE(spList)			
			i = spList(m)%row; j = spList(m)%col;
			! IF (i > nRows) THEN
				! PRINT *, 'i>nRows', i, nRows, m
			! END IF
			! IF (j > nCols) THEN
				! PRINT *, 'j>nCols', j, nCols, m
			! END IF
			matOut(i,:) = matOut(i,:) + (spList(m)%indVal)*matIn(j,:)
		END DO
	END FUNCTION
	
	FUNCTION scSPARSE_MATMUL_L(spList, matIn, nRows, nCols) RESULT(matOut)
		TYPE(sIndVal), INTENT(IN) :: spList(:)
		REAL*8, INTENT(IN) :: matIn(:,:)
		INTEGER, iNTENT(IN) :: nRows, nCols
		REAL*8 :: matOut(nRows,nCols)
		INTEGER :: i,j,m

		matOut = 0;
		
		! sparse mat * normal mat
		DO m = 1,SIZE(spList)			
			i = spList(m)%row; j = spList(m)%col;
			matOut(i,:) = matOut(i,:) + (spList(m)%indVal)*matIn(j,:)
		END DO
	END FUNCTION
	
	FUNCTION sSPARSE_MATMUL_R(matIn, spList, nRows, nCols) RESULT(matOut)
		TYPE(sIndVal), INTENT(IN) :: spList(:)
		REAL, INTENT(IN) :: matIn(:,:)
		INTEGER, INTENT(IN) :: nRows, nCols
		REAL :: matOut(nRows,nCols)
		INTEGER :: i,j,m
		
		matOut = 0;
		
		! normal mat * sparse mat
		DO m = 1,SIZE(spList)			
			i = spList(m)%row; j = spList(m)%col;
			matOut(:,j) = matOut(:,j) + (spList(m)%indVal)*matIn(:,i)
		END DO
	END FUNCTION
	
	FUNCTION cTRANSPOSE_SPARSE(spList) RESULT(spOut)
		TYPE(sIndVal), INTENT(IN) :: spList(:)
		INTEGER :: m
		TYPE(sIndVal), DIMENSION(SIZE(spList)) :: spOut
		
		DO m = 1,SIZE(spList)
			spOut(m) = sIndVal(spList(m)%col,spList(m)%row,spList(m)%indVal)
		END DO
	END FUNCTION
			
	
	FUNCTION sSPARSE_DIAG(spIn, diagVec) RESULT(spOut)
		! sparse * diag
		! elements along cols in spIn matrix are multiplied by the same number from diagVec
		! col index in spIn element is relevant
		TYPE(sIndVal), INTENT(IN) :: spIn(:)
		REAL, INTENT(IN) :: diagVec(:)
		TYPE(sIndVal) :: spOut(SIZE(spIn))
		INTEGER :: m, i, j, nPts
		
		nPts = SIZE(spIn)
		
		DO m = 1,nPts
			i = spIn(m)%row; j = spIn(m)%col
			spOut(m) = sIndVal(i, j, spIn(m)%indVal*diagVec(j))
		END DO	
	END FUNCTION
	
	FUNCTION sDIAG_SPARSE(diagVec, spIn) RESULT(spOut)
		! diag * sparse
		! elements along rows in spIn matrix are multiplied by the same number from diagVec
		! row index in spIn element is relevant
		TYPE(sIndVal), INTENT(IN) :: spIn(:)
		REAL, INTENT(IN) :: diagVec(:)
		TYPE(sIndVal) :: spOut(SIZE(spIn))
		INTEGER :: m, i, j, nPts
		
		nPts = SIZE(spIn)
		
		DO m = 1,nPts
			i = spIn(m)%row; j = spIn(m)%col
			spOut(m) = sIndVal(i, j, spIn(m)%indVal*diagVec(i))
		END DO	
	END FUNCTION
	
	FUNCTION SUM_SPARSE(spIn, matDim, sumDim) RESULT(vecOut)
		TYPE(sIndVal), INTENT(IN) :: spIn(:)
		INTEGER, INTENT(IN) :: matDim(2), sumDim
		REAL *8 :: vecOut(matDim(sumDim))
		
		INTEGER :: m, i, j, nPts
		
		nPts = SIZE(spIn)
		vecOut = 0
		
		IF (sumDim == 1) THEN  			! sum across rows (along column)
			DO m = 1,nPts
				i = spIn(m)%row;  j = spIn(m)%col;
				vecOut(j) = vecOut(j) + spIn(m)%indVal
			END DO
		ELSE IF (sumDim == 2) THEN		! sum across cols (along row)
			DO m = 1,nPts
				i = spIn(m)%row;  j = spIn(m)%col;
				vecOut(i) = vecOut(i) + spIn(m)%indVal
			END DO
		ELSE
			PRINT *, 'choose dimension of 1 (sum across rows), 2 (sum across cols)'
		END IF
	END FUNCTION
	
	! Note:  find LAPACK function explanations here: 
	! http://physics.oregonstate.edu/~landaur/nacphy/lapack/linear.html

	FUNCTION cINV(A) RESULT(Ainv)
	! taken from here: http://fortranwiki.org/fortran/show/Matrix+inversion
		REAL*8, DIMENSION(:,:), INTENT(IN) :: A
		REAL*8, DIMENSION(SIZE(A,1),SIZE(A,2)) :: Ainv
		
		REAL*8, DIMENSION(SIZE(A,1)) :: work ! work array for lapack
			! should be larger than N*NB, where NB is optimal blocksize returned by ILAENV
		INTEGER, DIMENSION(SIZE(A,1)) :: ipiv ! pivot indices
		INTEGER :: lda, info
		
		! External procedures defined in LAPACK
		EXTERNAL DGETRF
		EXTERNAL DGETRI
		
		! Store A in Ainv to prevent it from being overwritten
		Ainv = A;
		lda = SIZE(A,1)
		
		! CGETRF computes an LU factorization of a general M-by-N matrix A
		! using partial pivoting with row interchanges
		CALL DGETRF(lda, lda, Ainv, lda, ipiv, info)
		
		IF (info /= 0) THEN
			STOP 'Matrix is numerically singular'
		END IF
		
		! CGETRI computes the inverse of a matrix using the LU factorization
		! computed by CGETRF
		CALL DGETRI(lda, Ainv, lda, ipiv, work, lda, info)
		
		IF (info /= 0) THEN
			STOP 'Matrix inversion failed'
		END IF
	END FUNCTION cINV

	FUNCTION sINV(A) RESULT(Ainv)
		REAL, DIMENSION(:,:), INTENT(IN) :: A
		REAL, DIMENSION(SIZE(A,1),SIZE(A,2)) :: Ainv
		
		REAL, DIMENSION(SIZE(A,1)) :: work ! work array for lapack
		INTEGER, DIMENSION(SIZE(A,1)) :: ipiv ! pivot indices
		INTEGER :: lda, info
		
		! External procedures defined in LAPACK
		EXTERNAL DGETRF
		EXTERNAL DGETRI
		
		! Store A in Ainv to prevent it from being overwritten
		Ainv = A;
		lda = SIZE(A,1)
		
		! CGETRF computes an LU factorization of a general M-by-N matrix A
		! using partial pivoting with row interchanges
		CALL DGETRF(lda, lda, Ainv, lda, ipiv, info)
		
		IF (INFO /= 0) THEN
			STOP 'Matrix is numerically singular'
		END IF
		
		! CGETRI computes the inverse of a matrix using the LU factorization
		! computed by CGETRF
		CALL DGETRI(lda, Ainv, lda, ipiv, work, lda, info)
		
		IF (INFO /= 0) THEN
			STOP 'Matrix inversion failed'
		END IF
	END FUNCTION sINV

	FUNCTION sDET3(A) RESULT(det)
		REAL, DIMENSION(3,3) :: A
		REAL :: det
		
		det = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) &
				- A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) &
				+ A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))  
	END FUNCTION sDET3

	FUNCTION cAXB(A,B) RESULT(X)
		REAL*8, DIMENSION(:,:), INTENT(IN) :: A, B
		REAL*8, DIMENSION(SIZE(B,1),SIZE(B,2)) :: X
		
		REAL*8, DIMENSION(SIZE(A,1),SIZE(A,2)) :: A_LU  ! A after LU factorization
		! REAL, DIMENSION(SIZE(A,1)) :: work ! work array for lapack
		INTEGER, DIMENSION(SIZE(A,1)) :: ipiv ! pivot indices
		INTEGER :: lda, ldb, nrhs, info
		
		! External procedures defined in LAPACK
		EXTERNAL DGETRF
		EXTERNAL DGETRS
		
		! Store A,B in output matrices to prevent it from being overwritten
		A_LU = A;
		X = B;
		lda = SIZE(A,1)
		ldb = SIZE(B,1)
		nrhs = SIZE(B,2)
		
		! CGETRF computes an LU factorization of a general M-by-N matrix A
		! using partial pivoting with row interchanges
		CALL DGETRF(lda, lda, A_LU, lda, ipiv, info)
		
		! PRINT *, 'got pivots'
		
		IF (INFO > 0) THEN
			STOP 'Matrix is numerically singular'
		ELSE IF (INFO < 0) THEN
			PRINT *, 'element', INFO
			STOP 'element has illegal value'
		END IF
		
		! CGETRS computes A*X = B of a matrix using the LU factorization
		! computed by CGETRF
		CALL DGETRS('N', lda, nrhs, A_LU, lda, ipiv, X, ldb, info)
		! 'N' = no transpose, 'T' = transpose, 'C' = conjugate transpose
		
		IF (INFO /= 0) THEN
			STOP 'Solving AX=B failed'
		END IF
		! $OMP END PARALLEL
		
		! PRINT *, 'caxb done'
	END FUNCTION cAXB

	FUNCTION pseudoInv(A) RESULT(A_pinv)
			! get pseudoInverse using SVD: A = U Sigma V -> pinv(A) = V "SigmaInv" UT
		REAL*8, DIMENSION(:,:), INTENT(IN) :: A
		REAL*8, DIMENSION(SIZE(A,2),SIZE(A,1)) :: A_pinv
		
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: A_temp, SVD_U, SVD_UT, SVD_V, SVD_VT, SVD_SigmaInv
		REAL*8, DIMENSION(:), ALLOCATABLE :: sigmaVec, WORK
		
		INTEGER :: INFO, LWORK, r, mA, nA
		
		
		mA = SIZE(A,1)
		nA = SIZE(A,2)
		r = MIN(mA,nA)
		
		ALLOCATE(A_temp(mA,nA),SVD_U(mA,mA),SVD_UT(mA,mA),SVD_VT(nA,nA),SVD_V(nA,nA),SVD_SigmaInv(nA,mA),sigmaVec(r))
		
		LWORK = MAX(1,3*MIN(mA,nA)+MAX(mA,nA),5*MIN(mA,nA))*10
		ALLOCATE(WORK(MAX(1,LWORK)))
		
		A_temp = A;
		
		CALL DGESVD( 'a', 'a', mA, nA, A_temp, mA, sigmaVec, SVD_U, mA, SVD_VT, nA, WORK, LWORK, INFO )
		! see http://www.netlib.org/lapack/double/dgesvd.f/
		! A is an m x n matrix
		! arg1:  jobu = 'a' (all M columns of U are returned), 's' (first min(m,n) cols of U are returned (left singular vectors)), 
		!				'o' (first min(m,n) cols of U are overwritten on A), 'n' (no columns of U are computed)
		! 		 note that U(1:r) = colspace(A), U(r+1:m) = left nullspace(A), where r = rank
		! arg2:  jobv but for matrix V.  V(:,1:r) = rowspace(A), V(:,r+1:n) = nullspace(A), where r = rank  (what if r \= m,n?)
		! m:  number of rows in input matrix A
		! n:  number of cols in input matrix A
		! A:  input matrix
		! LDA:  leading dimension of array A (= #cols = rank)
		! S (output):  singular values of A, decreasing order
		! U (output):  dimension (LDU,UCOL)
		! LDU:  leading dimension of U (if JOBU = 's' or 'a', LDU >= M)
		! VT (output):  dimension(LDVT,N)
		! LDV:  leading dimension of array VT
		! work (workspace/output):  dimension(max(1,LWORK))
		! lwork:  dimension of array work
		! info:  = 0 -> successful,+ <0 = -i if illegal value at ith argument, >0 if didn't converge
		
		IF (INFO > 0) THEN
			STOP 'Bad convergence of SVD'
		ELSE IF (INFO < 0) THEN
			PRINT *, 'element', INFO
			STOP 'element has illegal value'
		END IF
		
		SVD_V = TRANSPOSE(SVD_VT)
		SVD_UT = TRANSPOSE(SVD_U)
		
		WHERE (ABS(sigmaVec/MAXVAL(sigmaVec)) > 1.0E-12) sigmaVec = 1.0/sigmaVec;
		SVD_SigmaInv = 0;
	
		SVD_SigmaInv(1:r,1:r) = cDIAG(sigmaVec(1:r))	! using the fact that zeros should be at the end
		
		A_pinv = MATMUL(SVD_V, MATMUL(SVD_SigmaInv,SVD_UT))
		
		DEALLOCATE(SVD_U,SVD_UT,SVD_VT,SVD_V,SVD_SigmaInv,sigmaVec)
		
	END FUNCTION

	FUNCTION getError(xGuess) RESULT(error)
		REAL*8, INTENT(IN) :: xGuess(:,:)
		REAL*8, DIMENSION(SIZE(xGuess,1),SIZE(xGuess,2)) :: error
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: errorX
		INTEGER :: m
		
		ALLOCATE(errorX(SIZE(xGuess,1),3))
		
		! find error of new guess
		DO m = 1,SIZE(xGuess,2)
			errorX(:,1) = xGuess(:,m) - 1.0;
			errorX(:,2) = 0.0 - xGuess(:,m);
			errorX(:,3) = 0.0;
			error(:,m) = MAXVAL(errorX,2);
		END DO
	
	END FUNCTION
	
	RECURSIVE FUNCTION minimizeError(xGuess, nullspace, cStep, errTol, maxNumIt) RESULT(X)
		REAL*8 :: xGuess(:,:), nullspace(:,:), X(SIZE(nullspace,1),1), cStep, errTol
		! xGuess = x_null + x_particular; assumed that it satisfies Ax= b
		! nullspace = vectors x_null along the columns
		! cstep = starting step size when evaluating jacobian
		
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: prevErrorX, errordX, newErrorX, dX, xNew, bestX, xParticular, jacobian, jacobian_pinv, &
												dCoeffs, prevErr, newErr, minErr, randCoeffs, mJacobian, mJacobian_pinv
		REAL*8 :: totPrevErr, totNewErr, totMinErr, cStepMax, cStepMin, &
					thresholdA, thresholdM, thresholdR, thresholdC, threshStrength
		INTEGER :: numIt, maxNumIt, nA, nullDim, m, maxErrInd, numItFound, randGuessIt, nullVecInd, numStep, maxNumStep
		LOGICAL :: ACCEPT

		
		! errTol = 1.0E-12
		threshStrength = 2.2
		
		nA = SIZE(nullspace,1)  ! # of rows in nullspace = # of elements in X ( = # of cols in A for AX=B)
		nullDim = SIZE(nullspace,2)  ! # of cols in nullspace
	
		ALLOCATE(prevErr(nA,1), errordX(nA,1), newErr(nA,1), minErr(nA,1), &
					dX(nA,1), xNew(nA,1), xParticular(nA,1), bestX(nA,1), &
					dCoeffs(nullDim,1), randCoeffs(nullDim,1), jacobian(nA,nullDim), jacobian_pinv(nullDim,nA))
		

		!!! algorithm: 
		!! get new error after least squares correction method
		!! probability of acceptance depends on change in error
		!!!! if error < tolerance, end loop
		!!!! if accept, generate new dCoeff for new x (smaller step size?)
		!!!! if reject, move only by a single step (dCoeff*cStep)  (cStep < 1)
		!!!!			find index of largest error
		!!!!			generate new dCoeff only changing coefficients of vectors that contribute the most to that index
		!! repeat.
		
		! assumes B is a 1D vector
	
		numIt = 0;  numStep = 0;
		maxNumStep = 2000;
		! maxNumIt = 20;	
		cStepMax = cStep;
		numItFound = 0;
		randGuessIt = 0;
		
		cStepMin = 1.0E-3;
		xParticular = xGuess;			! will define first guess to be the particular solution		

		prevErr = getError(xGuess)
		totPrevErr = SUM(prevErr)

		bestX = xGuess;
		minErr = prevErr
		totMinErr = totPrevErr
		
				
		! PRINT *, 'particular error'
		! PRINT *, MATMUL(mat1,xGuess)-mat2
		
		DO WHILE ((totPrevErr > errTol).AND.(numIt < maxNumIt).AND.(numStep < maxNumStep))
		
			prevErr = getError(xGuess)
			totPrevErr = SUM(prevErr)
			
			IF (totPrevErr < totMinErr) THEN
				totMinErr = totPrevErr
				minErr = prevErr;
				bestX = xGuess;
				numItFound = numIt;
			END IF
		
			! getting the Jacobian = -d residual/d coefficients
			DO m = 1, nullDim
				! looking a m-th column in nullSpace, and seeing the change in each element of dX
				dX(:,1) = nullSpace(:,m)*cStep
			
				! and the change in error associated with dX
				! errordX(:,1) = xGuess(:,1) + dX(:,1) - 1.0;
				! errordX(:,2) = 0.0 - (xGuess(:,1) + dX(:,1));
				! errordX(:,3) = 0.0;
				errordX = getError(xGuess + dX)
				jacobian(:,m) =  (errordX(:,1)-prevErr(:,1))/cStep;  ! nonlinear regression Jacobian. dresidual/dStep
				! PRINT *, 'delta error', MAXVAL(errordX,2)-prevErr(:,1)
					! change in error / step size for coeffs
					! diff => derivative
			END DO
			
			
			! desired dError = 0;  find dCoeffs that will get you there assuming linearity
			jacobian_pinv = pseudoInv(jacobian)
			dCoeffs = MATMUL(jacobian_pinv,0.0-prevErr);
			
			! then find new X and new error if used actually used those dCoeffs
			xNew = xGuess + MATMUL(nullspace,dCoeffs);	
				! don't need xParticular explicitly bc already included in initial xGuess
			
			newErr = getError(xNew);
			totNewErr = SUM(newErr)
			
			! PRINT *, 'new error', totNewErr, 'old error', totPrevErr
			
			
			
			IF (totNewErr < totMinErr) THEN
				totMinErr = totNewErr
				minErr = newErr;
				bestX = xNew;
				numItFound = numIt;
			END IF
			
			IF (totMinErr < errTol) EXIT
			
			! should sample random numbers from uniform? distribution to see if we accept or not (Ph129c)
			! continue wih smaller step size OR move on to a new point
			! with probability determined by change in error
			ACCEPT = .FALSE.;
			thresholdA = 1.0 - MINVAL((/1.0,REAL(0.8*EXP((totNewErr-totPrevErr)/(totNewErr+totPrevErr)*threshStrength))/));
				! when totNewErr < totPrevErr, threshold is high; want high probability of accepting result.
			
			! can be any value
			thresholdM = totNewErr/(totNewErr+totPrevErr);
				! when totNewErr >> totPrevErr, thresholdM is close to 1; want high probability of finer mesh
			! thresholdC = 1.0 - 50*cStep/cStepMax; ! - 1.0*(numIt-randGuessIt)/randGuessIt
				! ! want to make educated guess on guessing new x0 when cStep gets too small (thresholdC is high)
				! ! (and is not making much significant improvements)
			thresholdR = 7.0*(MINVAL((/totPrevErr,totNewErr/)) - totMinErr)/(MINVAL((/totPrevErr,totNewErr/)) + totMinErr);
				! change point if stepping leads to divergence... 
				! avoid switching too much when minErr is small
				! usually when totNewErr < totPrevErr, xNew is accepted so no random guessing
			
			! PRINT *, 'A', thresholdA, 'M', thresholdM, 'R', thresholdR
			
			IF (RAND() < thresholdA) ACCEPT = .TRUE.
			! IF (totNewErr < totPrevErr) ACCEPT = .TRUE.
			
			IF (ACCEPT) THEN
				! just update position; use xNew to find another xGuess
				! PRINT *, 'accepted: new error', totNewErr, 'old error', totPrevErr
				cStep = MAXVAL( (/ MINVAL((/cStep/2.0, totNewErr/)), cStepMin /) );
				xGuess = minimizeError(xNew, nullspace, cStep, errTol, maxNumIt-1)

			ELSE		! effectively ignoring xNew; guessing next xGuess
				IF ((thresholdM >= thresholdR).OR.(totNewErr > 100*totPrevErr)) THEN
				! sudden spike in error might mean that we're actually close to the completely reflecting solution...
					
					! move to spot closer to original guess
					xGuess = xGuess + MATMUL(nullspace,dCoeffs*cStep);	
					numStep = numStep + 1;
					
				ELSE ! completely random new xGuess
					! PRINT *, 'random xGuess: new error', totNewErr, 'old error', totPrevErr
					DO m = 1,nullDim
						dCoeffs(m,1) = RAND()		! scaling times random number
					END DO
					xGuess = xParticular + MATMUL(nullspace,dCoeffs)

					cStep = cStepMax
					
					randGuessIt = numIt
					numStep = numStep + 1;
					! PRINT *, 'after random xGuess: new error', totNewErr, 'old error', totPrevErr
				END IF
				
			END IF
			
			numIt = numIt + 1;
			
		END DO
		
		! PRINT *, 'out of do loop'
		X = bestX
		! PRINT *, 'END FCT: min error', totMinErr
		
		DEALLOCATE(errordX,prevErr,newErr,minErr, &
					dCoeffs,dX,xNew,xParticular,bestX, jacobian,jacobian_pinv, randCoeffs)
		
	END FUNCTION
	

	FUNCTION cAXB_null(A,B,minX,maxX) RESULT(X)
			! finds x, accounting for singular (underdefined) A
		REAL*8, DIMENSION(:,:), INTENT(IN) :: A, B
		REAL*8, DIMENSION(SIZE(A,2),SIZE(B,2)) :: X
		REAL*8, INTENT(IN) :: minX, maxX		! constraints on values in output vector
		
		REAL*8, DIMENSION(:,:), ALLOCATABLE :: A_temp, B_temp, A_pinv, SVD_V, SVD_VT, SVD_U, SVD_UT, SVD_SigmaInv, &
										xParticular, xGuess, nullSpace_pinv, nullSpace, &
										coeffMat, dCoeff, dError_pinv, bestCoeff, bestXGuess
		REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: errorX, errorX_new, dError, cCorr, currentErr, dX
		REAL*8 :: cStep, totError, error_tol, smallestErr
		! REAL, DIMENSION(SIZE(A,1)) :: work ! work array for lapack
		REAL*8, DIMENSION(:), ALLOCATABLE :: WORK, sigmaVec, nullSigmaVec
		INTEGER, DIMENSION(:), ALLOCATABLE :: fixed_pivs, free_pivs
		INTEGER :: mA,nA,mB,nB,r,nullDim, i,j,m, LWORK, INFO, numIt,maxNumIt
		
		
		! External procedures defined in LAPACK
		EXTERNAL DGESVD		! computes singular value decomposition of rectangular matrix
		EXTERNAL DGETRF

		
		mA = SIZE(A,1)
		nA = SIZE(A,2)
		mB = SIZE(B,1)
		nB = SIZE(B,2)
	
		
		!!!!! perform SVD to get colspace, nullspace, and rank
		! SVD:  A = U*Sigma*transpose(V)
		
		ALLOCATE(A_temp(mA,nA),B_temp(mB,nB),A_pinv(nA,mA),fixed_pivs(mA),free_pivs(nA),SVD_U(mA,mA),SVD_UT(mA,mA),&
							SVD_VT(nA,nA),SVD_V(nA,nA),SVD_SigmaInv(nA,mA),sigmaVec(MIN(mA,nA)),xParticular(nA,nB))
		
		A_temp = A;
		B_temp = B;
		
		LWORK = MAX(1,3*MIN(mA,nA)+MAX(mA,nA),5*MIN(mA,nA))*10
		ALLOCATE(WORK(MAX(1,LWORK)))
		
		CALL DGESVD( 'a', 'a', mA, nA, A_temp, mA, sigmaVec, SVD_U, mA, SVD_VT, nA, WORK, LWORK, INFO )
		! see http://www.netlib.org/lapack/double/dgesvd.f
		! A is an m x n matrix
		! arg1:  jobu = 'a' (all M columns of U are returned), 's' (first min(m,n) cols of U are returned (left singular vectors)), 
		!				'o' (first min(m,n) cols of U are overwritten on A), 'n' (no columns of U are computed)
		! 		 note that U(1:r) = colspace(A), U(r+1:m) = left nullspace(A), where r = rank
		! arg2:  jobv but for matrix V.  V(:,1:r) = rowspace(A), V(:,r+1:n) = nullspace(A), where r = rank  (what if r \= m,n?)
		! m:  number of rows in input matrix A
		! n:  number of cols in input matrix A
		! A:  input matrix
		! LDA:  leading dimension of array A (= #cols = rank)
		! S (output):  singular values of A, decreasing order
		! U (output):  dimension (LDU,UCOL)
		! LDU:  leading dimension of U (if JOBU = 's' or 'a', LDU >= M)
		! VT (output):  dimension(LDVT,N)
		! LDV:  leading dimension of array VT
		! work (workspace/output):  dimension(max(1,LWORK))
		! lwork:  dimension of array work
		! info:  = 0 -> successful,+ <0 = -i if illegal value at ith argument, >0 if didn't converge
		
		IF (INFO > 0) THEN
			STOP 'Bad convergence of SVD'
		ELSE IF (INFO < 0) THEN
			PRINT *, 'element', INFO
			STOP 'element has illegal value'
		END IF
		
		SVD_V = TRANSPOSE(SVD_VT)
		SVD_UT = TRANSPOSE(SVD_U)
		
		PRINT *, 'SVD success'
		
		! PRINT *, 'sigmas', sigmaVec
		
		! get dimension of nullspace (and hence rank)
		r = 0;
		nullDim = 0;
		DO i = 1, mA
			IF (ABS(sigmaVec(i)/MAXVAL(sigmaVec)) > 1.0E-12) THEN 
				r = r + 1;
				fixed_pivs(r) = i
			ELSE
				! note, all elements in sigma should be positive, and elements are in decreasing order.. so zeros should be last
				nullDim = nullDim + 1;
				free_pivs(nullDim) = i
			END IF
		END DO
		
		free_pivs(nullDim+1:nA-r) = (/ (m, m=mA,nA) /);	! cols of V (rows of VT) corresponding to 0s in svd_sigma
		nullDim = nA - r;
		PRINT *, 'nullDim', nullDim, 'rank', r
		
		ALLOCATE(nullSpace(nA,nullDim))
		! nullSpace = SVD_V(:,free_pivs(1:nullDim))
		nullSpace = SVD_V(:,r+1:nA)
		
		PRINT *, 'A*nullspace'
		error_tol = SUM(ABS(MATMUL(A,nullSpace)))
		error_tol = MAXVAL( (/error_tol,1.0E-7/) )
		PRINT *, 'error_tol', error_tol
		
		
		!!!!! solve for particular solution using X = A^(-1) Y = V SigmaInv UT Y = pinv(A) Y
		A_pinv = pseudoInv(A)
		xParticular = MATMUL(A_pinv,B)
		! PRINT *, 'xParticular'
		! PRINT *, xParticular
		
		! PRINT *, 'particular error'
		! PRINT *, MATMUL(A,xParticular)-B
		
		smallestErr = 10;
			
		IF (nullDim > 0) THEN
			
			cStep = 0.1
			maxNumIt = 200;
			PRINT *, 'cStep, maxNumIt', cStep, maxNumIt
			X = minimizeError(xParticular, nullSpace, cStep, error_tol, maxNumIt)
			
		ELSE
			! matrix is not singular
			X = cAXB(A,B)
		END IF
			
		DEALLOCATE(SVD_V,SVD_VT,SVD_U,SVD_UT,SVD_SigmaInv,A_pinv,WORK,sigmaVec)
		DEALLOCATE(xParticular, nullSpace, A_temp, B_temp, fixed_pivs, free_pivs)
		
		PRINT *, 'caxb null done'
		
		WRITE (*,"(20(F10.3))") X
		
		PRINT *, 'bc error'
		PRINT *, MATMUL(A,X)-B
		
	END FUNCTION cAXB_null
	
	
END MODULE CONSTANTS