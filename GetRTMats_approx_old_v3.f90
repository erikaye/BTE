SUBROUTINE GetRTMats(matlParams, Ls, write_path, &
				R12_d,T12_d,R21_d,T21_d,g1mCoeff_d,g2pCoeff_d, &
				R12_s,T12_s,R21_s,T21_s,g1mCoeff_s,g2pCoeff_s, &
				numG12d,numG21d,numR12d,numR21d,numT12d,numT21d,&
				numG12s,numG21s,numR12s,numR21s,numT12s,numT21s, fBoth, numFBoth, weightsMatch, spec1, spec2)
! these matrices are represented by a list of indices and their values (row, col, val)
! in the case that lattice 1 = lattice 2 at interface, then we can either do specular or diffuse scattering (one slab or 2 slabs)
! in the case that lattice 1 != lattice 2, only do diffuse scattering bc of lattice mismatch
! for now assume that crystal orientation is the same.. st comparing ky and kz will give you this information

! newer version modified from GetRTMats_dbal
! if lattice mismatch, assume all modes are mapped to (though not equally transmitted to)
! if specular, can have k-conserved scattering

USE Constants
USE omp_lib

IMPLICIT NONE

TYPE(material), DIMENSION(2), INTENT(IN) :: matlParams
REAL, DIMENSION(2), INTENT(IN) :: Ls
CHARACTER(longStr), INTENT(IN) :: write_path
TYPE(sIndVal), INTENT(OUT) :: R12_d(:), T12_d(:), R21_d(:), T21_d(:), g1mCoeff_d(:), g2pCoeff_d(:), &
					R12_s(:), T12_s(:), R21_s(:), T21_s(:), g1mCoeff_s(:), g2pCoeff_s(:)
INTEGER, INTENT(OUT) :: numG12d,numG21d,numR12d,numR21d,numT12d,numT21d,numG12s,numG21s,numR12s,numR21s,numT12s,numT21s, numFBoth, fBoth(:,:)
REAL*8, INTENT(OUT) :: spec1(:), spec2(:), weightsMatch(:)
LOGICAL :: filesExist

REAL, ALLOCATABLE :: kx1(:), kx2(:), ky1(:), ky2(:), kz1(:), kz2(:), dkx1(:), dkx2(:), Kn1(:), Kn2(:), &
					freq1(:), freq2(:), weights1(:), weights2(:), uniqueFreqs1(:), uniqueFreqs2(:)
REAL :: df1, df2, rho1, rho2
INTEGER, ALLOCATABLE :: numEachFreq1(:), numEachFreq2(:)

INTEGER :: w1, w2, w1d, w2d, i, j, m1, m2, indNum1, indNum2
! INTEGER, DIMENSION(:,:), ALLOCATABLE :: fBoth
REAL, DIMENSION(:), ALLOCATABLE :: sumRowG12_s, sumRowR12_d,sumRowR12_s,sumRowT12_d,sumRowT12_s, &
									sumRowG21_s, sumRowR21_d,sumRowR21_s,sumRowT21_d,sumRowT21_s, &
									R1d,R1s,T1d,T1s,R2d,R2s,T2d,T2s,G1d,G1s,G2d,G2s
DOUBLE COMPLEX :: Zin
REAL*8, ALLOCATABLE :: Gamma(:), ZSL(:)
REAL*8 :: Z0, ZL, kxpt
INTEGER :: numInt, m
									! rWeights1, rWeights2, cWeights1, cWeights2  ! scaled weighting (s-right, ss-left)
REAL, DIMENSION(:), ALLOCATABLE :: rmu1_d,rmu1_s,tmu1_d,tmu1_s, &
									rmu2_d,rmu2_s,tmu2_d,tmu2_s, &
									gmu1_d,gmu1_s,gmu2_d,gmu2_s
									
REAL :: kylatt1, kylatt2, kzlatt1, kzlatt2, ksig, fsig1, fsig2, fsig
REAL*8 :: weightT
REAL, ALLOCATABLE :: k2List(:)
REAL*8, ALLOCATABLE :: k2weightsT(:), tempWeights1(:), tempWeights2(:)
INTEGER :: numMatchK
LOGICAL :: sameLatt

REAL*8, ALLOCATABLE :: oneVec(:), redG(:,:), spec1d(:,:), spec1Temp(:,:)
CHARACTER(5) :: LString



w1 = matlParams(1)%wTot;
w2 = matlParams(2)%wTot;
w1d = matlParams(1)%numUniqueFreq;
w2d = matlParams(2)%numUniqueFreq;

ALLOCATE(kx1(w1),dkx1(w1),ky1(w1),kz1(w1), freq1(w1),Kn1(w1),weights1(w1), uniqueFreqs1(w1d),numEachFreq1(w1d))
ALLOCATE(kx2(w2),dkx2(w2),ky2(w2),kz2(w2), freq2(w2),Kn2(w2),weights2(w2), uniqueFreqs2(w2d),numEachFreq2(w2d))

kx1 = matlParams(1)%kx(1:w1);
dkx1 = matlParams(1)%dkx(1:w1);
ky1 = matlParams(1)%ky(1:w1);
kz1 = matlParams(1)%kz(1:w1);
freq1 = matlParams(1)%freq(1:w1)
Kn1 = matlParams(1)%MFP(1:w1) / Ls(1)
weights1 = matlParams(1)%weights(1:w1)
uniqueFreqs1 = matlParams(1)%uniqueFreqs(1:w1d)
numEachFreq1 = matlParams(1)%numEachFreq(1:w1d)
df1 = matlParams(1)%df
rho1 = matlParams(1)%rho
spec1 = matlParams(1)%spec

PRINT *, w1, w2, SIZE(kx2)

kx2 = matlParams(2)%kx(1:w2);
dkx2 = matlParams(2)%dkx(1:w2);
ky2 = matlParams(2)%ky(1:w2);
kz2 = matlParams(2)%kz(1:w2);
freq2 = matlParams(2)%freq(1:w2)
Kn2 = matlParams(2)%MFP(1:w2) / Ls(2)
weights2 = matlParams(2)%weights(1:w2)
uniqueFreqs2 = matlParams(2)%uniqueFreqs(1:w2d)
numEachFreq2 = matlParams(2)%numEachFreq(1:w2d)
df2 = matlParams(2)%df
rho2 = matlParams(2)%rho
spec2 = matlParams(2)%spec

PRINT *, 'RT mats'

kylatt1 = MAXVAL(ABS(ky1));
kylatt2 = MAXVAL(ABS(ky2));
kzlatt1 = MAXVAL(ABS(kz1));
kzlatt2 = MAXVAL(ABS(kz2));
sameLatt = ((kylatt1 == kylatt2) .AND. (kzlatt1 == kzlatt2));

PRINT *, 'klatts', kylatt1, kylatt2, kzlatt1, kzlatt2


PRINT *, 'w1', w1, 'w2', w2

! ALLOCATE(rWeights1(w1),rWeights2(w2),cWeights1(w1),cWeights2(w2))
! rWeights1 = 1.0
! rWeights2 = 1.0
! cWeights1 = 1.0
! cWeights2 = 1.0

! arrays with info on sum of non-zero element weights in each row of R12, R21, T12, T21
ALLOCATE(sumRowG12_s(w1), sumRowR12_d(w1),sumRowR12_s(w1),sumRowT12_d(w1),sumRowT12_s(w1), &  
		 sumRowG21_s(w2), sumRowR21_d(w2),sumRowR21_s(w2),sumRowT21_d(w2),sumRowT21_s(w2))

! arrays describing |reflection/transmission| of each mode
ALLOCATE(R1d(w1), R1s(w1), T1d(w1), T1s(w1), R2d(w2), R2s(w2), T2d(w2), T2s(w2), &
		 G1d(w1), G1s(w1), G2d(w2), G2s(w2), ZSL(w1), Gamma(w1))

! diffuse case:  energy is conserved, but tangential momentum is not
! include weighting along column (diag(weights)*mat)  
! (different weights for each row, same weights for each column)

sumRowG12_s = 0; sumRowR12_d = 0; sumRowR12_s = 0; sumRowT12_d = 0; sumRowT12_s = 0;  ! sum over mu2
sumRowG21_s = 0; sumRowR21_d = 0; sumRowR21_s = 0; sumRowT21_d = 0; sumRowT21_s = 0;

numR12d = 0;  numR21d = 0;
numT12d = 0;  numT21d = 0;

PRINT *, 'start'

! reflection matrices for the two materials 
DO m1 = 1,w1d
	indNum1 = SUM(numEachFreq1(1:m1-1))
	DO i = indNum1+1,indNum1+numEachFreq1(m1)  ! row
		IF (i > w1) THEN
			PRINT *, 'R12 ind1 > w1', i
		END IF
		DO j = indNum1+1,indNum1+numEachFreq1(m1)  ! col
			IF (j > w1) THEN
				PRINT *, 'R12 ind1 > w1', i
			END IF
			numR12d = numR12d+1;
			sumRowR12_d(i) = sumRowR12_d(i) + 1.0;  ! getting num elements along row
			! R12_d(numR12d) = sIndVal(i,j,weights1(j))   ! same weights for each row
			R12_d(numR12d) = sIndVal(i,j,1.0)   ! same weights for each row
			! using weights helps us approximately take degeneracies into account for weighting of how much each mode gets transmitted to
			
		END DO
	END DO
END DO

PRINT *, 'loop1'
DO m2 = 1,w2d
	indNum2 = SUM(numEachFreq2(1:m2-1))
	DO i = indNum2+1,indNum2+numEachFreq2(m2)  ! row
		IF (i > w2) THEN
			PRINT *, 'R21 ind1 > w2', i
		END IF
		DO j = indNum2+1,indNum2+numEachFreq2(m2)  ! col
			IF (j > w2) THEN
				PRINT *, 'R21 ind2 > w2', i
			END IF
			numR21d = numR21d+1;
			sumRowR21_d(i) = sumRowR21_d(i) + 1.0;  ! getting num elements along row
			! R21_d(numR21d) = sIndVal(i,j,weights2(j))
			R21_d(numR21d) = sIndVal(i,j,1.0)
		END DO
	END DO
END DO



numFBoth = 0;
fBoth = 0;
weightsMatch = 0;

R1d = 1.0;  R2d = 1.0;  ! initially assume that all modes are mismatched at interface

! T12_d, T21_d
! use gaussian smearing to determine transmission probability from one frequency to another
weightT = 0.0;

IF (NOT(sameLatt)) THEN  
			! different lattice constants at interface, different materials (assume all Lindsay-like data)
			! all materials are fcc though
	DO m1 = 1, w1d  ! row / col
		DO m2 = 1, w2d  ! col / row
			indNum1 = SUM(numEachFreq1(1:m1-1));
			indNum2 = SUM(numEachFreq2(1:m2-1));
			
			fsig = (df1+df2)/2/2
		
			weightT = 1./(fsig*SQRT(2*pi))*EXP(-1.0/2.0*(uniqueFreqs1(m1)-uniqueFreqs2(m2))**2/fsig**2)
		
			IF (weightT > fgs_tol/(fsig*SQRT(2*pi))) THEN     ! if we predict non-negligible transmission
				numFBoth = numFBoth + 1;
				 
				fBoth(numFBoth, 1:2) = (/ m1, m2 /);    ! index pairs (in unique freq)
				weightsMatch(numFBoth) = weightT;
				
				DO i = indNum1+1, indNum1+numEachFreq1(m1)   ! material 1
					IF (i > w1) THEN
						PRINT *, 'T12 ind1 > w1', i
					END IF
					
					DO j = indNum2+1, indNum2+numEachFreq2(m2)   ! material 2
						IF (j > w2) THEN
							PRINT *, 'T12 ind2 > w2', j
						END IF
						numT12d = numT12d + 1;
						sumRowT12_d(j) = sumRowT12_d(j) + 1.0
						! T12_d(numT12d) = sIndVal(j,i,weightT*weights2(i))		! different weights (dmu1)*weightT for each row
						T12_d(numT12d) = sIndVal(j,i,weightT)		! different weights (dmu1)*weightT for each row
						
						numT21d = numT21d + 1;
						sumRowT21_d(i) = sumRowT21_d(i) + 1.0   
						! T21_d(numT21d) = sIndVal(i,j,weightT*weights1(j))  		! different weights (dmu2)*weightT for each row
						T21_d(numT21d) = sIndVal(i,j,weightT)  		! different weights (dmu2)*weightT for each row
						
						R1d(i) = 0.0          ! there can be transmission through these modes at interface
						R2d(j) = 0.0
					END DO
				END DO			
			END IF
		END DO
	END DO
ELSE
	! treat this case separately bc uniquefreqs were defined differently (via sharp cut-off)
	DO m1 = 1, w1d  ! row / col
		indNum1 = SUM(numEachFreq1(1:m1-1));
		
		numFBoth = numFBoth + 1;
			
		fBoth(numFBoth, 1:2) = (/ m1, m1 /);    ! index pairs (in unique freq)
		weightsMatch(numFBoth) = 1.0;
			
		DO i = indNum1+1, indNum1+numEachFreq1(m1)   ! material 1
			IF (i > w1) THEN
				PRINT *, 'T12 ind1 > w1', i
			END IF
			
			DO j = indNum1+1, indNum1+numEachFreq2(m1)   ! material 2
				IF (j > w2) THEN
					PRINT *, 'T12 ind2 > w2', j
				END IF
				numT12d = numT12d + 1;
				! sumRowT12_d(j) = sumRowT12_d(j) + weights2(i)
				! T12_d(numT12d) = sIndVal(j,i,weights2(i))  ! bc mat1 = mat2, numfreq1 = numfreq2
				sumRowT12_d(j) = sumRowT12_d(j) + 1.0
				! T12_d(numT12d) = sIndVal(j,i,weights1(i))  
				T12_d(numT12d) = sIndVal(j,i,1.0)  
				
				
				numT21d = numT21d + 1;
				! sumRowT21_d(i) = sumRowT21_d(i) + weights1(j)
				! T21_d(numT21d) = sIndVal(i,j,weights1(j))
				sumRowT21_d(i) = sumRowT21_d(i) + 1.0
				! T21_d(numT21d) = sIndVal(i,j,weights2(j))
				T21_d(numT21d) = sIndVal(i,j,1.0)
				
				R1d(i) = 0.0          ! there can be transmission through these modes at interface
				R2d(j) = 0.0
			END DO
		END DO			
	END DO
END IF	
	
PRINT *, 'ex fSig', fSig
PRINT *, 'ex dW', uniqueFreqs1(w1d) - uniqueFreqs1(w1d-1), uniqueFreqs2(w2d) - uniqueFreqs2(w2d-1)

PRINT *, 'diffuse mats'
PRINT *, 'numFBoth:', numFBoth
PRINT *, 'maxNumBoth:', w1d, w2d
PRINT *, 'numT12:', numT12d


!!!!!!!!!!!!!!!!!!!!!!!!!!
! specular case: energy is conserved, tangential momentum is also conserved
! no mode conversion
! ignoring crystal momentum for T12, T21 for now


numR12s = 0;  numT12s = 0; 
numR21s = 0;  numT21s = 0;

! only loop through points in which energy (frequency) is conserved
! R12_s
IF (sameLatt) THEN  
	DO m1 = 1, w1d
		indNum1 = SUM(numEachFreq1(1:m1-1))
		! lattice is matched 
		DO i = indNum1 + 1, indNum1 + numEachFreq1(m1)
			DO j = indNum1 + 1, indNum1 + numEachFreq1(m1)
				! ! IF ((ABS(ky1(j) - ky1(i))/ky1(i) <= ktol) .AND. (ABS(kz1(j)-kz1(i))/kz1(i) <= ktol)) THEN
                IF ((ky1(j) == ky1(i)) .AND. (kz1(j) == kz1(i))) THEN
				!!! also matching kx means that no mode conversion is allowed (k_tot 1 = k_tot 2)
				!!! to allow for mode conversion, only match ky, kz
				!!! cross matching ky1=kz2, ky2=kz1 should've been done in orderparams (remove degeneracy)
				! IF ((kx1(j) == kx1(i)) .AND.(ky1(j) == ky1(i)) .AND. (kz1(j) == kz1(i))) THEN
					numR12s = numR12s + 1;

					sumRowG12_s(i) = sumRowG12_s(i) + 1.0;
					! sumRowR12_s(i) = sumRowR12_s(i) + &
						! 0.25*ABS((1.0-kx1(j))*(1.0-kx1(i))/SQRT(kx1(j)**2 + ky1(j)**2 + kz1(j)**2)/SQRT(kx1(i)**2 + ky1(i)**2 + kz1(i)**2))
					! R12_s(numR12s) = sIndVal(i,j,&
						! 0.25*ABS((1.0-kx1(j))*(1.0-kx1(i))/SQRT(kx1(j)**2 + ky1(j)**2 + kz1(j)**2)/SQRT(kx1(i)**2 + ky1(i)**2 + kz1(i)**2)))
					sumRowR12_s(i) = sumRowR12_s(i) + 1.0;
					R12_s(numR12s) = sIndVal(i,j,1.0)
				END IF
			END DO
		END DO
	END DO
ELSE
	! in the case of lattice mismatch, maps to all modes of same frequency, which is the diffuse case
	numR12s = numR12d;
	sumRowR12_s = sumRowR12_d;
	R12_s = R12_d;  
END IF


! R21_s
IF (sameLatt) THEN  
	DO m2 = 1, w2d
		indNum2 = SUM(numEachFreq2(1:m2-1))
		DO i = indNum2 + 1, indNum2 + numEachFreq2(m2)  ! row
			DO j = indNum2 + 1, indNum2 + numEachFreq2(m2)  ! col
				! ! IF ((ABS(ky2(j) - ky2(i))/ky2(i) <= ktol) .AND. (ABS(kz2(j)-kz2(i))/kz2(i) <= ktol)) THEN
				IF ((ky2(j) == ky2(i)) .AND. (kz2(j) == kz2(i))) THEN
				! IF ((kx2(j) == kx2(i)) .AND. (ky2(j) == ky2(i)) .AND. (kz2(j) == kz2(i))) THEN
					numR21s = numR21s + 1;
					sumRowG21_s(i) = sumRowG21_s(i) + 1.0;
					! sumRowR21_s(i) = sumRowR21_s(i) + &
						! 0.25*ABS((1.0-kx2(i))*(1.0-kx2(j))/SQRT(kx2(j)**2 + ky2(j)**2 + kz2(j)**2)/SQRT(kx2(i)**2 + ky2(i)**2 + kz2(i)**2))
					! R21_s(numR21s) = sIndVal(i,j,&
						! 0.25*ABS((1.0-kx2(i))*(1.0-kx2(j))/SQRT(kx2(j)**2 + ky2(j)**2 + kz2(j)**2)/SQRT(kx2(i)**2 + ky2(i)**2 + kz2(i)**2)))
					sumRowR21_s(i) = sumRowR21_s(i) + 1.0;
					R21_s(numR21s) = sIndVal(i,j,1.0)
				END IF
			END DO
		END DO
	END DO
ELSE
	numR21s = numR21d;
	sumRowR21_s = sumRowR21_d;
	R21_s = R21_d;
END IF	

PRINT *, 'spec Rs'

! T12_s, T21_s

! R1s = 1.0;  R2s = 1.0;      ! initially assume no modes can go through interface
numMatchK = 0

IF (sameLatt) THEN
	PRINT *, 'same latt'
	DO m1 = 1, numFBoth
		indNum1 = SUM(numEachFreq1(1:fBoth(m1,1)-1))
		indNum2 = SUM(numEachFreq2(1:fBoth(m1,2)-1))
		
		DO i = (indNum1 + 1), (indNum1 + numEachFreq1(fBoth(m1,1)))  ! material 1
			! call function to match ky2 to ky1, and kz2 to kz1
			! generates list of ky/z2 (indices) for each ky/z1 that can be transmitted to
			
			ALLOCATE(k2List(numEachFreq2(fBoth(m1,2)))) !,k2weightsT(numEachFreq2(fBoth(m1,2))))
			
			CALL matchKs(ky1(i),kz1(i),kylatt1, kylatt2, kzlatt1, kzlatt2, &
						ky2((indNum2+1):(indNum2+numEachFreq2(fBoth(m1,2)))), &
						kz2((indNum2+1):(indNum2+numEachFreq2(fBoth(m1,2)))), &
						indNum2, k2List, numMatchK)
			
			DO j = 1,numMatchK
			
                ! IF (1==1) THEN
				IF (kx1(i) == kx2(k2List(j))) THEN
				
                    numT12s = numT12s + 1;
                    sumRowT12_s(k2List(j)) = sumRowT12_s(k2List(j)) + 1.0
                    T12_s(numT12s) = sIndVal(k2List(j),i,1.0)
                    ! sumRowT12_s(k2List(j)) = sumRowT12_s(k2List(j)) + &
                            ! ABS(kx1(i)*kx2(k2List(j))/SQRT(kx1(i)**2 + ky1(i)**2 + kz1(i)**2)&
                                    ! /SQRT(kx1(k2List(j))**2 + ky1(k2List(j))**2 + kz1(k2List(j))**2))
                    ! T12_s(numT12s) = sIndVal(k2List(j),i,&
                            ! ABS(kx1(i)*kx2(k2List(j))/SQRT(kx1(i)**2 + ky1(i)**2 + kz1(i)**2) &
                                    ! /SQRT(kx1(k2List(j))**2 + ky1(k2List(j))**2 + kz1(k2List(j))**2)))
                    
                    numT21s = numT21s + 1;
                    sumRowT21_s(i) = sumRowT21_s(i) + 1.0  ! sum along col, over rows. size w2
                    T21_s(numT21s) = sIndVal(i,k2List(j),1.0)
                    ! sumRowT21_s(i) = sumRowT21_s(i) + &
                            ! ABS(kx1(i)*kx2(k2List(j))/SQRT(kx1(i)**2 + ky1(i)**2 + kz1(i)**2)&
                                    ! /SQRT(kx1(k2List(j))**2 + ky1(k2List(j))**2 + kz1(k2List(j))**2))
                    ! T21_s(numT21s) = sIndVal(i,k2List(j), &
                            ! ABS(kx1(i)*kx2(k2List(j))/SQRT(kx1(i)**2 + ky1(i)**2 + kz1(i)**2)&
                                    ! /SQRT(kx1(k2List(j))**2 + ky1(k2List(j))**2 + kz1(k2List(j))**2)))
                    
                    R1s(i) = 0.0		! there can be transmission through these modes at interface
                    R2s(j) = 0.0
					
				END IF			
			END DO
			
			DEALLOCATE(k2List)
			
		END DO
		
	END DO
ELSE
	numT12s = numT12d;
	sumRowT12_s = sumRowT12_d;
	T12_s = T12_d;
	
	numT21s = numT21d;
	sumRowT21_s = sumRowT21_d;
	T21_s = T21_d;
	
	R1s = R1d;
	R2s = R2d;
END IF

PRINT*, 'no kx cons: mode conversion allowed (for 0)'
PRINT*, 'kx cons: no mode conversion allowed (for R12,R21,T)'

PRINT *, numT12s, numT21s, numR21s, numR12s, numT12d, numT21d, numR21d, numR12d

! now add values for reflection, transmission
ALLOCATE(rmu1_d(w1),rmu1_s(w1),tmu2_d(w1),tmu2_s(w1), &
		 rmu2_d(w2),rmu2_s(w2),tmu1_d(w2),tmu1_s(w2), &
		 gmu1_d(w1),gmu1_s(w1),gmu2_d(w2),gmu2_s(w2))
		 
rmu1_d = 0; rmu1_s = 0; tmu2_d = 0; tmu2_s = 0;
rmu2_d = 0; rmu2_s = 0; tmu1_d = 0; tmu1_s = 0;
gmu1_d = 0; gmu1_s = 0; gmu2_d = 0; gmu2_s = 0;

! set reflection coefficients for each row
! values along rows are constant (but different for each row (for T mats at least))

! not all modes across two materials are matched... R = 1.
! R1d = 1.0;         R1s = 1.0          ! size(w1) (R1d,R1s defined above in loops)
! can do R1d = MIN(1.0, R1d + frequency dependent shape of transmission) to define all coeffs.

! R1d = MIN(1.0, R1d + 1.0/(1+EXP((Kn1-1)/(roughness+1.0E-12))));		! use fermi dirac fct, short MFP < L1 reflected, large MFP transmitted 
! R1s = MIN(1.0, R1s + 1.0/(1+EXP((Kn1-1)/(roughness+1.0E-12))));

R1d = 0.5;
R1s = 0.0;

PRINT *, 'unity trans (AMM), DMM'

!!! QWP reflectivity

!! original version
! Z0 = 1.0;  ZL = 0.5; Gamma = (Z0-ZL)/(Z0+ZL);

! !! new test (2/17/18)
! Z0 = 1.0; ZL = 0.5;  
! ZSL = (Z0-ZL)/2.0/(1.0+(2.0*pi/kx1/Ls(1)))+Z0
! Gamma = (ZSL-ZL)/(ZSL+ZL)
! ! if wavelength << L, then Z_superlattice looks like Z_0
! ! if mfp and wavelength >> L, then Z_superlattice looks like average of Z_0 and Z_L

! ! ! Zin = Z0*(ZL+DCMPLX(0,Z0*TAN(kx1*Ls(1)/2.0/pi)))/(Z0+DCMPLX(0,ZL*TAN(kx1*Ls(1)/2.0/pi)))
! ! Zin = Z0*(ZL+DCMPLX(0,Z0*TAN(kx1*1.0E-8/2.0/pi)))/(Z0+DCMPLX(0,ZL*TAN(kx1*1.0E-8/2.0/pi)))
! ! R1s = ((REAL(Zin)-Z0)/(REAL(Zin)+Z0))**2
! ! R1s = R1s/MAXVAL(R1s)
! !! this normalization is also potentially wrong..? but we also did a non-normalized version

! !! do numerical integration to find average reflectivity for each mode
! R1s = 0.0;
! DO m1 = 1,w1
	! ! numInt = 10*(INT(dkx1(m1)*Ls(2)/2.0/pi)+1);
	! ! numInt = 10*(INT(dkx1(m1)*0.9E-9/2.0/pi)+1);
	! numInt = 10*(INT(dkx1(m1)*1.0E-8/2.0/pi)+1);
	! DO m = 1,numInt
		! kxpt = kx1(m1) + dkx1(m1)*m/REAL(numInt)
		
		! !!!! first equation (which was wrong)
		! ! Zin(1) = Z0*(ZL+DCMPLX(0,Z0*TAN(kxpt*Ls(2)/2.0/pi)))/(Z0+DCMPLX(0,ZL*TAN(kxpt*Ls(2)/2.0/pi)))
		! ! ! Zin(1) = Z0*(ZL+DCMPLX(0,Z0*TAN(kxpt*1.0E-8/2.0/pi)))/(Z0+DCMPLX(0,ZL*TAN(kxpt*1.0E-8/2.0/pi)))
		! ! R1s(m1) = R1s(m1) + (((REAL(Zin(1))-Z0)/(REAL(Zin(1))+Z0))**2)/REAL(numInt)  ! average reflectivity
		! ! ! by def, Z1 is actually a vector of size w1. but here, we're using it as a scalar
		
		! Zin = ZL*((1+Gamma(m1)*EXP(DCMPLX(0,2*kxpt*Ls(2))))/(1-Gamma(m1)*EXP(DCMPLX(0,2*kxpt*Ls(2)))))
		! ! Zin(1) = ZL*((1+Gamma*EXP(DCMPLX(0,2*kxpt*0.9E-9)))/(1-Gamma*EXP(DCMPLX(0,2*kxpt*0.9E-9))))
		! ! Zin = ZL*((1+Gamma(m1)*EXP(DCMPLX(0,2*kxpt*1.0E-8)))/(1-Gamma(m1)*EXP(DCMPLX(0,2*kxpt*1.0E-8))))
		! R1s(m1) = R1s(m1) + (ABS((Zin-ZSL(m1))/(Zin+ZSL(m1)))**2)/REAL(numInt)  ! average reflectivity
		
	! END DO
! END DO

! PRINT *, 'QWP, ZSL:', Ls(2),Ls(1)

! ! ! R1s = R1s/MAXVAL(R1s)
! ! ! uncommented in QWPint run.  WRONG!

! R1s = 0;
! WHERE(ABS(kx1/2.0/pi/Ls(1) - 1) < 0.25) R1s = 1;


! PRINT *, 'L', 0.9E-9, 'sumR1', SUM(R1s)
! PRINT *, 'L', Ls(2), 'sumR1', SUM(R1s)
! PRINT *, 'L', 1.0E-8, 'sumR1', SUM(R1s)


!! angle dependent reflectivity
! R1s = ABS(kx1/SQRT(kx1**2 + ky1**2 + kz1**2))			! unphysical
! PRINT *, 'unphysical angle dep reflectivity'
! R1s = 1.0-ABS(kx1/SQRT(kx1**2 + ky1**2 + kz1**2))		! physical
! PRINT *, 'physical angle dep reflectivity'

! !! frequency dep reflectivity
! R1s = 0
! WHERE(freq1 < MAXVAL(freq1)*0.6) R1s = 1;
! PRINT *, 'frequency dependent reflectivity (high pass)'

! ! ! WHERE(freq1 < MAXVAL(freq1)*0.6) R1s = 0;
! WHERE(freq1 > MAXVAL(freq1)*0.3) R1s = 1.0;
! PRINT *, 'spec frequency dependent reflectivity (low pass)'
! WHERE(freq1 > MAXVAL(freq1)*0.3) R1d = 1.0;
! PRINT *, 'diff frequency dependent reflectivity (low pass)'

! WHERE(freq1 < MAXVAL(freq1)*0.3) R1s = 0.05*(1.0-ABS(kx1/SQRT(kx1**2 + ky1**2 + kz1**2)))		! physical;
! PRINT *, 'spec angle dependent reflectivity (low pass) angle*0.05'


R1d = MIN(1.0, R1d + 0.0);		! diffuse matrices -> R = T = 0.5 (DMM for Si/Si interface)
R1s = MIN(1.0, R1s + 0.0);


! write R1s to a txt file
WRITE(LString, '(i5)') INT(Ls(2)*1E9)
OPEN(5,FILE = './output_R1s.dat')
WRITE(5,"(40000(ES15.7,','))") R1s
CLOSE(5)

T1d = 1.0 - R1d;   T1s = 1.0 - R1s     ! size(w1)
G1d = 1.0;         G1s = 1.0;
rmu1_d = R1d/(sumRowR12_d);  
rmu1_s = R1s/(sumRowR12_s);  ! size(w1), r_mu1_mu1' for each (non-zero) element in row for mu1
gmu1_d = G1d/(sumRowR12_d);  ! sumRowR12_d = sumRowG12_d
gmu1_s = G1s/(sumRowR12_s);
tmu1_d = T1d/(sumRowT12_d+1.0E-30);  
tmu1_s = T1s/(sumRowT12_s+1.0E-30);  ! size(w2), t_mu1_mu2 for each (non-zero) element in row mu1 (sum over mu2)

WHERE(sumRowT12_d == 0) tmu1_d = 0.0
WHERE(sumRowT12_s == 0) tmu1_s = 0.0
WHERE(sumRowR12_d == 0) rmu1_d = 0.0
WHERE(sumRowR12_s == 0) rmu1_s = 0.0


T2d = SUM(sSPARSE_MATMUL_R(cDIAG(tmu1_d),T12_d(1:numT12d),w2,w1),1)        ! size(w1).  sum along cols of T12 = sum along rows of T21.
T2s = SUM(sSPARSE_MATMUL_R(cDIAG(tmu1_s),T12_s(1:numT12s),w2,w1),1)        ! can't do SUM(t12Mat*T12,2) because the weights applied so far made it non-symmetric
! sum elements in column -> dim = 1
! sum elements in row -> dim = 2


R2d = 1.0 - T2d;   R2s = 1.0 - T2s;    ! size(w2).  should be consistent with what modes were perfectly reflecting, found above.
G2d = 1.0;         G2s = 1.0;
rmu2_d = R2d/(sumRowR21_d);  rmu2_s = R2s/(sumRowR21_s);   ! size(w2)
gmu2_d = G2d/(sumRowR21_d);  gmu2_s = G2s/(sumRowR21_s);
tmu2_d = T2d/(sumRowT21_d+1.0E-30);  tmu2_s = T2s/(sumRowT21_s+1.0E-30);   ! size(w1), sum over mu2

PRINT *, SUM(weights1)


! reflected mode should at least have itself so never sumRowR21 is never 0
WHERE(sumRowT21_d == 0) tmu2_d = 0.0
WHERE(sumRowT21_s == 0) tmu2_s = 0.0
WHERE(sumRowR21_d == 0) rmu2_d = 0.0
WHERE(sumRowR21_s == 0) rmu2_s = 0.0

! write R2s to a txt file
OPEN(5,FILE = './output_R2d.dat')
WRITE(5,"(40000(ES15.7,','))") R2d
CLOSE(5)

OPEN(5,FILE = './output_T2d.dat')
WRITE(5,"(40000(ES15.7,','))") T2d
CLOSE(5)


numG12d = numR12d; numG21d = numR21d;
numG12s = numR12s; numG21s = numR21s;


! definition of 1->2 matrices
g1mCoeff_d(1:numR12d) = sDIAG_SPARSE(gmu1_d,R12_d(1:numR12d));
g1mCoeff_s(1:numR12s) = sDIAG_SPARSE(gmu1_s,R12_s(1:numR12s));
R12_d(1:numR12d) = sDIAG_SPARSE(rmu1_d,R12_d(1:numR12d));  ! list of (row,col,val)
R12_s(1:numR12s) = sDIAG_SPARSE(rmu1_s,R12_s(1:numR12s));
T12_d(1:numT12d) = sDIAG_SPARSE(tmu1_d,T12_d(1:numT12d)); 
T12_s(1:numT12s) = sDIAG_SPARSE(tmu1_s,T12_s(1:numT12s));

! ! ! satisfy T12(mu1,mu2) = T21(mu2,mu1)
! T21_d(1:numT21d) = cTRANSPOSE_SPARSE(T12_d(1:numT12d))   ! numT21s,d = numT12s,d
! T21_s(1:numT21s) = cTRANSPOSE_SPARSE(T12_s(1:numT12s))

! T21_d = T12_d;
! T21_s = T12_s;

! definition of 2->1 matrices
g2pCoeff_d(1:numR21d) = sDIAG_SPARSE(gmu2_d,R21_d(1:numR21d));
g2pCoeff_s(1:numR21s) = sDIAG_SPARSE(gmu2_s,R21_s(1:numR21s));
R21_d(1:numR21d) = sDIAG_SPARSE(rmu2_d,R21_d(1:numR21d)); 
R21_s(1:numR21s) = sDIAG_SPARSE(rmu2_s,R21_s(1:numR21s));
T21_d(1:numT21d) = sDIAG_SPARSE(tmu2_d,T21_d(1:numT21d)); 
T21_s(1:numT21s) = sDIAG_SPARSE(tmu2_s,T21_s(1:numT21s));

! ! weighting along col
! g1mCoeff_d(1:numR12d) = sSPARSE_DIAG(g1mCoeff_d(1:numR12d),weights1);
! g1mCoeff_s(1:numR12s) = sSPARSE_DIAG(g1mCoeff_s(1:numR12s),weights1);
! R12_d(1:numR12d) = sSPARSE_DIAG(R12_d(1:numR12d),weights1);  ! list of (row,col,val)
! R12_s(1:numR12s) = sSPARSE_DIAG(R12_s(1:numR12s),weights1);
! T12_d(1:numT12d) = sSPARSE_DIAG(T12_d(1:numT12d),weights1); 
! T12_s(1:numT12s) = sSPARSE_DIAG(T12_s(1:numT12s),weights1);

! g2pCoeff_d(1:numR21d) = sSPARSE_DIAG(g2pCoeff_d(1:numR21d),weights2);
! g2pCoeff_s(1:numR21s) = sSPARSE_DIAG(g2pCoeff_s(1:numR21s),weights2);
! R21_d(1:numR21d) = sSPARSE_DIAG(R21_d(1:numR21d),weights2); 
! R21_s(1:numR21s) = sSPARSE_DIAG(R21_s(1:numR21s),weights2);
! T21_d(1:numT21d) = sSPARSE_DIAG(T21_d(1:numT21d),weights2);
! T21_s(1:numT21s) = sSPARSE_DIAG(T21_s(1:numT21s),weights2);

! spec1 = Kn2**2/(Kn2**2+1.0)			! Kn >> 1 --> ballistic --> specular; Kn << 1 --> diffusive --> diffuse
! spec2 = Kn2**2/(Kn2**2+1.0)
! PRINT *, 'L2 dependent spec 1, Kn2**2'
!! this is ok bc Kn isn't used anywhere else in this code
! spec1 = 1.0;
! spec2 = 1.0;
! PRINT *, 'spec 1'
PRINT *, 'original spec'

! ALLOCATE(spec1d(w1,1),spec1Temp(w1,1))
! spec1Temp(1:w1,1) = spec1;
! spec1d = scSPARSE_MATMUL_L(g1mCoeff_d(1:numR12d),spec1Temp,w1,1)
! spec1 = spec1d(:,1)
! DEALLOCATE(spec1d,spec1Temp)
! ALLOCATE(spec1d(w2,1),spec1Temp(w2,1))
! spec1Temp(1:w2,1) = spec2;
! spec1d = scSPARSE_MATMUL_L(g2pCoeff_d(1:numR21d),spec1Temp,w2,1)
! spec2 = spec1d(:,1)
! PRINT *, 'avged spec'
! DEALLOCATE(spec1d, spec1Temp)
PRINT *, 'no avg'

PRINT *, SIZE(numEachFreq1)


!!!!!!!!!!!!!!!!!!!!! WRITING R/T FILES !!!!!!!!!!!!!!!!!!!!!!!!
! OPEN(5,FILE = TRIM(write_path) // 'spec1.dat')
! WRITE(5,"(40000(ES16.7E3,','))") spec1
! CLOSE(5)

! OPEN(5,FILE = TRIM(write_path) // 'spec2.dat')
! WRITE(5,"(40000(ES16.7E3,','))") spec2
! CLOSE(5)

! ! CALL OMP_INIT_LOCK(LCK)

INQUIRE(FILE = TRIM(write_path) // 'R12s.dat', EXIST = filesExist)
filesExist = .FALSE.

IF (filesExist == .FALSE.) THEN    ! write files

! !$OMP PARALLEL NUM_THREADS(4) DEFAULT(SHARED) PRIVATE(m1)

! !$OMP SINGLE
! WRITE (*,*) 'Parallel part num threads: ', OMP_GET_NUM_THREADS()
! !$OMP END SINGLE

! !$OMP SECTIONS
	! !$OMP SECTION		
		OPEN(1,FILE = TRIM(write_path) // 'R12s.dat')
		WRITE(1,*) w1, w1, numR12s
		DO m1 = 1,numR12s
			WRITE(1,*) R12_s(m1)%row, R12_s(m1)%col, R12_s(m1)%indVal
		END DO
		CLOSE(1)


	! !$OMP SECTION
		OPEN(2,FILE = TRIM(write_path) // 'T12s.dat')
		WRITE(2,*) w2, w1, numT12s
		DO m1 = 1,numT12s
			WRITE(2,*) T12_s(m1)%row, T12_s(m1)%col, T12_s(m1)%indVal
		END DO
		CLOSE(2)
	
	! !$OMP SECTION
		OPEN(13,FILE = TRIM(write_path) // 'G12s.dat')
		WRITE(13,*) w1, w1, numG12s
		DO m1 = 1,numG12s
			WRITE(13,*) g1mCoeff_s(m1)%row, g1mCoeff_s(m1)%col, g1mCoeff_s(m1)%indVal
		END DO
		CLOSE(13)
	
	! !$OMP SECTION
		OPEN(4,FILE = TRIM(write_path) // 'R21s.dat')
		WRITE(4,*) w2, w2, numR21s
		DO m1 = 1,numR21s
			WRITE(4,*) R21_s(m1)%row, R21_s(m1)%col, R21_s(m1)%indVal
		END DO
		CLOSE(4)

	! !$OMP SECTION
		OPEN(8,FILE = TRIM(write_path) // 'T21s.dat')
		WRITE(8,*) w1, w2, numT21s
		DO m1 = 1,numT21s
			WRITE(8,*) T21_s(m1)%row, T21_s(m1)%col, T21_s(m1)%indVal
		END DO
		CLOSE(8)
		
	! !$OMP SECTION
		OPEN(7,FILE = TRIM(write_path) // 'G21s.dat')
		WRITE(7,*) w2, w2, numG21s
		DO m1 = 1,numG21s
			WRITE(7,*) g2pCoeff_s(m1)%row, g2pCoeff_s(m1)%col, g2pCoeff_s(m1)%indVal
		END DO
		CLOSE(7)
		
	! !$OMP SECTION
		! ! write R12_d to a txt file _OK
		OPEN(11,FILE = TRIM(write_path) // 'R12d.dat')
		
		WRITE(11,*) w1, w1, numR12d
		DO m1 = 1,numR12d
			WRITE(11,*) R12_d(m1)%row, R12_d(m1)%col, R12_d(m1)%indVal
		END DO
		CLOSE(11)

	! !$OMP SECTION
		! ! write T12_d to a txt file _OK
		OPEN(12,FILE = TRIM(write_path) // 'T12d.dat')
		WRITE(12,*) w2, w1, numT12d
		DO m1 = 1,numT12d
			WRITE(12,*) T12_d(m1)%row, T12_d(m1)%col, T12_d(m1)%indVal
		END DO
		CLOSE(12)
	
	! !$OMP SECTION
		! ! write G12_d to a txt file _OK
		OPEN(3,FILE = TRIM(write_path) // 'G12d.dat')
		WRITE(3,*) w1, w1, numG12d
		DO m1 = 1,numG12d
			WRITE(3,*) g1mCoeff_d(m1)%row, g1mCoeff_d(m1)%col, g1mCoeff_d(m1)%indVal
		END DO
		CLOSE(3)
	
	! !$OMP SECTION
		! ! write R21_s to a txt file
		OPEN(14,FILE = TRIM(write_path) // 'R21d.dat')
		WRITE(14,*) w2, w2, numR21d
		DO m1 = 1,numR21d
			WRITE(14,*) R21_d(m1)%row, R21_d(m1)%col, R21_d(m1)%indVal
		END DO
		CLOSE(14)
		
	! !$OMP SECTION
		OPEN(18,FILE = TRIM(write_path) // 'T21d.dat')
		WRITE(18,*) w1, w2, numT21d
		DO m1 = 1,numT21d
			WRITE(18,*) T21_d(m1)%row, T21_d(m1)%col, T21_d(m1)%indVal
		END DO
		CLOSE(18)
		
	! !$OMP SECTION
		! ! write G21_s to a txt file
		OPEN(17,FILE = TRIM(write_path) // 'G21d.dat')
		WRITE(17,*) w2, w2, numG21d
		DO m1 = 1,numG21d
			WRITE(17,*) g2pCoeff_d(m1)%row, g2pCoeff_d(m1)%col, g2pCoeff_d(m1)%indVal
		END DO
		CLOSE(17)
		
	! !$OMP SECTION
		! ! write f-matching data to a txt file 
		OPEN(9,FILE = TRIM(write_path) // 'fMatch.dat')
		WRITE(9,*) numFBoth
		DO m1 = 1,numFBoth
			WRITE(9,*) fBoth(m1,1), fBoth(m1,2), weightsMatch(m1)
		END DO
		CLOSE(9)

		! OPEN(5,FILE = TRIM(write_path) // 'R1s.dat')
		! WRITE(5,"(40000(ES15.7,','))") R1s
		! CLOSE(5)
		
		OPEN(5,FILE = TRIM(write_path) // 'spec1.dat')
		WRITE(5,"(40000(ES15.7,','))") spec1
		CLOSE(5)

		! OPEN(5,FILE = TRIM(write_path) // 'spec1d.dat')
		! WRITE(5,"(40000(ES15.7,','))") spec1d
		! CLOSE(5)
		
		OPEN(10,FILE = TRIM(write_path) // 'spec2.dat')
		WRITE(10,"(40000(ES15.7,','))") spec2
		CLOSE(10)
		
		OPEN(10,FILE = TRIM(write_path) // 'kx1.dat')
		WRITE(10,"(40000(ES15.7,','))") kx1
		CLOSE(10)
		
		OPEN(10,FILE = TRIM(write_path) // 'kx2.dat')
		WRITE(10,"(40000(ES15.7,','))") kx2
		CLOSE(10)
		
	! !$OMP END SECTIONS
! !$OMP END PARALLEL 
	PRINT *, 'wrote R/T files'
ELSE 
	PRINT *, 'did NOT write R/T files'
END IF

DEALLOCATE(R1d,R1s,T1d,T1s,R2d,R2s,T2d,T2s,G1d,G1s,G2d,G2s, ZSL,Gamma, &
			rmu1_d,rmu1_s,tmu1_d,tmu1_s,rmu2_d,rmu2_s,tmu2_d,tmu2_s, &
			gmu1_d,gmu1_s,gmu2_d,gmu2_s)
			
DEALLOCATE(sumRowG12_s,sumRowR12_d,sumRowR12_s,sumRowT12_d,sumRowT12_s, &
			sumRowG21_s,sumRowR21_d,sumRowR21_s,sumRowT21_d,sumRowT21_s)
			
CONTAINS

	SUBROUTINE matchKs(ky1f,kz1f,kylatt1f,kylatt2f, kzlatt1f, kzlatt2f, &
							ky2f, kz2f, indNum2f, kMatchListf, numMatchf)
		REAL, INTENT(IN) :: ky1f, kz1f, kylatt1f, kylatt2f, kzlatt1f, kzlatt2f, ky2f(:), kz2f(:)
		INTEGER, INTENT(IN) :: indNum2f
		REAL, DIMENSION(:), INTENT(OUT) :: kMatchListf
		INTEGER, INTENT(OUT) :: numMatchf
		INTEGER :: nmf, jf, maxNMf
		REAL, DIMENSION(:), ALLOCATABLE :: tempf, temp2f, myListf, mzListf
		
		numMatchf = 0;
		
		maxNMf = 500;
		
		ALLOCATE(tempf(maxNMf*2+1), temp2f(maxNMf*2+1), myListf(maxNMf*2+1), mzListf(maxNMf*2+1))

		IF (kylatt1f /= kylatt2f) THEN
			
			DO jf = 1, SIZE(ky2f)
				! finding "m" that satisfies
				! k2-k1 = 2pi/klatt2*(n-m) + (klatt2-klatt1)*m
				tempf = (/ (kylatt2f*nmf, nmf = -maxNMf,maxNMf) /) - ABS(ky2f(jf)-ky1f)
				temp2f = tempf/(kylatt1f-kylatt2f)
				myListf = ABS(temp2f-NINT(temp2f))
				
				tempf = (/ (kzlatt2f*nmf, nmf = -maxNMf,maxNMf) /) - ABS(kz2f(jf)-kz1f)
				temp2f = tempf/(kzlatt1f-kzlatt2f)
				mzListf = ABS(temp2f-NINT(temp2f))
			END DO
		
			IF (ANY(myListf <= 0.01).AND.ANY(mzListf <= 0.01)) THEN
				! IF (myList(j) <= ktol) THEN
				numMatchf = numMatchf + 1;
				kMatchListf(numMatchf) = jf + indNum2f  ! index of first kpoint in ky2, kz2.
			END IF			
		ELSE
			! kylatt1 = kylatt2 when transmitting to the same material
			! ignores mode conversion for reflection (which needs interface k vector)
			! which would be the other case when kylatt1 = kylatt2
			myListf = ABS(ky2f-ky1f)/ky1f
			mzListf = ABS(kz2f-kz1f)/kz1f
			
			DO jf = 1, SIZE(ky2f)
				IF ((myListf(jf) < ktol).AND.(mzListf(jf) < ktol)) THEN
					numMatchf = numMatchf + 1;
					kMatchListf(numMatchf) = jf + indNum2f
				END IF
			END DO
		END IF
		
		DEALLOCATE(tempf,temp2f,myListf,mzListf)
		
	END SUBROUTINE matchKs
	
END SUBROUTINE GetRTMats