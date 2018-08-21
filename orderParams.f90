SUBROUTINE orderParams(kx,ky,kz,dkx,dky,dkz,vgx,vgy,vgz,tau,freq,C,MFP,weights,gamma, &
						numUniqueFreq, numEachFreq, uniqueFreqs, numPts, wTot, df, dataset)
! orders modes by frequency, adds together degenerate modes


USE Constants
IMPLICIT NONE

REAL, DIMENSION(:) :: kx,ky,kz,dkx,dky,dkz,vgx,vgy,vgz,tau,freq,C,MFP,weights
REAL*8, DIMENSION(:) :: gamma
INTEGER, INTENT(IN) :: numPts
INTEGER :: wTot
CHARACTER(4), INTENT(IN) :: dataset

INTEGER :: i, j, w, m, mm, indO, indNum1
INTEGER, INTENT(OUT) :: numUniqueFreq, numEachFreq(:)
REAL, INTENT(OUT) :: uniqueFreqs(:), df

INTEGER :: numKx, numKy, numKz, ind
REAL, DIMENSION(:), ALLOCATABLE :: kxTemp, kxUnique, kyTemp, kyUnique, kzTemp, kzUnique, dkxU, dkyU, dkzU
INTEGER, DIMENSION(:), ALLOCATABLE :: kxInd, kyInd, kzInd, kxIndAll, kyIndAll, kzIndAll, dkxIndAll, dkyIndAll, dkzIndAll

REAL, DIMENSION(:), ALLOCATABLE :: kxO,kyO,kzO,dkxO,dkyO,dkzO,vgxO,vgyO,vgzO,tauO,freqO,CO,MFPO,weightsO
REAL, DIMENSION(:), ALLOCATABLE :: kxF,kyF,kzF,dkxF,dkyF,dkzF,vgxF,vgyF,vgzF,tauF,freqF,CF,MFPF,weightsF
REAL*8, DIMENSION(:), ALLOCATABLE :: gammaO, gammaF
INTEGER, DIMENSION(:), ALLOCATABLE :: indToRemove, numElDegen, indList

REAL :: freqTemp


PRINT *, 'order params hoare quicksort'

! material parameters from DFT for 2 adjacent blocks
w = SIZE(freq)  ! = wTot

! frequencies for Al, Mingo selected by vgx > 0, so k's can be negative!
! find and sort unique k's to get dk's
! note which index elements in kCol map to

ALLOCATE(kxTemp(numPts),kxUnique(numPts),kyTemp(numPts),kyUnique(numPts),kzTemp(numPts),kzUnique(numPts),&
			kxInd(numPts),kyInd(numPts),kzInd(numPts),dkxU(numPts),dkyU(numPts),dkzU(numPts),&
			kxIndAll(wTot),kyIndAll(wTot),kzIndAll(wTot),dkxIndAll(wTot),dkyIndAll(wTot),dkzIndAll(wTot))
			
kxUnique = 1E30; kyUnique = 1E30; kzUnique = 1E30; kxInd = 0; kyInd = 0; kzInd = 0; dkxU = 0; dkyU = 0; dkzU = 0;

kxTemp = kx(1:numPts);
numKx = 0;
DO WHILE (ANY(kxTemp < 1E30)) 
	ind = MINLOC(kxTemp,1)
	numKx = numKx + 1;
	kxUnique(numKx) = kxTemp(ind);	! kUnique should be sorted in ascending order
	WHERE(kx(1:numPts) == kxUnique(numKx)) kxTemp = 1E30;
	WHERE(kx(1:numPts) == kxUnique(numKx)) kxInd = numKx;
END DO

kyTemp = ky(1:numPts);
numKy = 0;
DO WHILE (ANY(kyTemp < 1E30)) 
	ind = MINLOC(kyTemp,1)
	numKy = numKy + 1;
	kyUnique(numKy) = kyTemp(ind);	! kUnique should be sorted in ascending order
	WHERE(ky(1:numPts) == kyUnique(numKy)) kyTemp = 1E30;
	WHERE(ky(1:numPts) == kyUnique(numKy)) kyInd = numKy;
END DO

kzTemp = kz(1:numPts);
numKz = 0;
DO WHILE (ANY(kzTemp < 1E30)) 
	ind = MINLOC(kzTemp,1)
	numKz = numKz + 1;
	kzUnique(numKz) = kzTemp(ind);	! kUnique should be sorted in ascending order
	WHERE(kz(1:numPts) == kzUnique(numKz)) kzTemp = 1E30;
	WHERE(kz(1:numPts) == kzUnique(numKz)) kzInd = numKz;
END DO

PRINT *, numKx, numKy, numKz

! dk's for unique kx/ky/kz
dkxU(2:numKx) = kxUnique(2:numKx)-kxUnique(1:numKx)
dkyU(2:numKy) = kyUnique(2:numKy)-kyUnique(1:numKy)
dkzU(2:numKz) = kzUnique(2:numKz)-kzUnique(1:numKz)
IF (dataset .EQ. 'Lind') THEN
	dkxU(1) = MINVAL(ABS(kxUnique(1:numKx)));   
	dkyU(1) = MINVAL(ABS(kyUnique(1:numKy)));   
	dkzU(1) = MINVAL(ABS(kzUnique(1:numKz)));   
ELSE
	dkxU(1) = dkxU(2);
	dkyU(1) = dkyU(2);
	dkzU(1) = dkzU(2);
END IF

PRINT *, dkxU(1), dkyU(1), dkzU(1)
	! for Al, Mingo, all dk should be the same = dk closest to 0.  for Lindsay, MINVAL(kxUnique) is the one closest to 0

! index to list of unique k's for all modes (consistent with how paramCol is defined in getMatlParams)
kxIndAll = RESHAPE(kxInd,(/wTot/),PAD=kxInd)
kyIndAll = RESHAPE(kyInd,(/wTot/),PAD=kyInd)
kzIndAll = RESHAPE(kzInd,(/wTot/),PAD=kzInd)

! dk for all modes
dkx = dkxU(kxIndAll); dky = dkyU(kyIndAll); dkz = dkzU(kzIndAll);
df = MINVAL(ABS(vgx*dkx + vgy*dky + vgz*dkz))
PRINT *, 'delta freq', df


! ! ! ! now to sort pts by frequency !(and order pts by increasing frequency)
ALLOCATE(kxO(w),kyO(w),kzO(w),dkxO(w),dkyO(w),dkzO(w),vgxO(w),vgyO(w),vgzO(w),tauO(w),&
			freqO(w),CO(w),MFPO(w),weightsO(w),gammaO(w))
			
ALLOCATE(kxF(w),kyF(w),kzF(w),dkxF(w),dkyF(w),dkzF(w),vgxF(w),vgyF(w),vgzF(w),tauF(w),&
			freqF(w),CF(w),MFPF(w),weightsF(w),gammaF(w))
			
ALLOCATE(indList(w))

			
! get list of indices of modes with increasing frequency
! does frequency need to be in increasing order? or just grouped together...
			
			
numEachFreq = 0;
numUniqueFreq = 0;
indO = 0;

! put data in order of increasing frequency
indList = (/ (m, m = 1,w) /)
CALL QUICKSORT(freq,indList,1,w)
! indList and freq are now sorted by increasing frequency

! check success:

PRINT *, 'ordered freq?', ALL((freq(2:w)-freq(1:w-1)) >= 0)

kxF = kx(MOD(indList-1,numPts)+1);
kyF = ky(MOD(indList-1,numPts)+1);
kzF = kz(MOD(indList-1,numPts)+1);
dkxF = dkx(indList)
dkyF = dky(indList)
dkzF = dkz(indList)
vgxF = vgx(indList)
vgyF = vgy(indList)
vgzF = vgz(indList)
tauF = tau(indList)
CF = C(indList)
MFPF = MFP(indList)
weightsF = weights(indList)
gammaF = gamma(indList)

! PRINT *, 'size input', SIZE(vgx), SIZE(kx), SIZE(gamma), w, wTot

! PRINT *, 'sumks', SUM(kx(1:numPts)), SUM(ky(1:numPts)), SUM(kz(1:numPts))
! PRINT *, 'sumkFs', SUM(kxF)/6, SUM(kyF)/6, SUM(kzF)/6
! PRINT *, 'sumdks', SUM(dkx), SUM(dky), SUM(dkz)
! PRINT *, 'sumdkFs', SUM(dkxF), SUM(dkyF), SUM(dkzF)
! PRINT *, 'vg', SUM(vgx), SUM(vgy), SUM(vgz)
! PRINT *, 'vgF', SUM(vgxF), SUM(vgyF), SUM(vgzF)
! PRINT *, 'others', SUM(tau), SUM(C), SUM(MFP), SUM(weights), SUM(gamma)
! PRINT *, 'others F', SUM(tauF), SUM(CF), SUM(MFPF), SUM(weightsF), SUM(gammaF)

! PRINT *, 'indList', indList(1:10)
! PRINT *, 'indList mod', MOD(indList-1,5)+1

! find unique freqs, group them together (not necessarily in increasing frequency order)
DO i = 1, w
	IF (freq(i) /= 0) THEN
		freqTemp = freq(i)
		numUniqueFreq = numUniqueFreq + 1;
		uniqueFreqs(numUniqueFreq) = freqTemp		
		DO j = i,w
			! because we're searching the same sequence, we can reduce # of points we search
			IF (ABS(freq(j) - freqTemp) < ABS(df)) THEN
				numEachFreq(numUniqueFreq) = numEachFreq(numUniqueFreq) + 1;
				indO = indO + 1;
				freqO(indO) = freq(j)
				! freqO(indO) = freqTemp  ! artificially pin s.t. these points have the same frequency
				freq(j) = 0;
					! so that it's not double counted
					
				! kxO(indO) = kx(MOD(j-1,numPts)+1)
				! kyO(indO) = ky(MOD(j-1,numPts)+1)
				! kzO(indO) = kz(MOD(j-1,numPts)+1)
				kxO(indO) = kxF(j)
				kyO(indO) = kyF(j)
				kzO(indO) = kzF(j)
				! dkxO(indO) = dk(kxInd(MOD(j-1,numPts)+1))
				! dkyO(indO) = dk(kyInd(MOD(j-1,numPts)+1))
				! dkzO(indO) = dk(kzInd(MOD(j-1,numPts)+1))
				! kxO(indO) = kx(j)
				! kyO(indO) = ky(j)
				! kzO(indO) = kz(j)
				dkxO(indO) = dkxF(j)
				dkyO(indO) = dkyF(j)
				dkzO(indO) = dkzF(j)
				vgxO(indO) = vgxF(j);  
				vgyO(indO) = vgyF(j); 
				vgzO(indO) = vgzF(j);
				tauO(indO) = tauF(j);
				CO(indO) = CF(j); 
				! CO(indO) = C(i);  ! bc it's determined by freq
				MFPO(indO) = MFPF(j);
				weightsO(indO) = weightsF(j);
				gammaO(indO) = gammaF(j);
			END IF
		END DO
	END IF
END DO


PRINT *, 'num Unique freq', numUniqueFreq
PRINT *, 'wTot', wTot

! check
IF (w /= indO) THEN
	PRINT *, 'w', w, 'sumEachFreq', indO
ELSE
	PRINT *, 'sorting is all ok'
END IF



PRINT *, 'end ordering, start removing degeneracies'

! only loop through points in which energy (frequency) is conserved

ALLOCATE(indToRemove(w), numElDegen(numUniqueFreq))

numElDegen = 0
indToRemove = 0;

DO m = 1, numUniqueFreq
	indNum1 = SUM(numEachFreq(1:m-1))
	
	DO i = indNum1 + 1, indNum1 + numEachFreq(m)
		IF (weightsO(i) /= 0) THEN		! skip indices that are going to be removed to avoid double counting
			DO j = i+1, indNum1 + numEachFreq(m)
				! if ky(i),kz(i) = ky(j),kz(j) (rotational symmetry--can't distinguish btwn ky,kz) and if kx(i) = kx(i)
				IF ((((ABS(kyO(j) - kyO(i))/kyO(i) <= ktol) .AND. (ABS(kzO(j)-kzO(i))/kzO(i) <= ktol)) .OR. &
						((ABS(kyO(j) - kzO(i))/kyO(i) <= ktol) .AND. (ABS(kzO(j)-kyO(i))/kzO(i) <= ktol))) .AND. &
						(ABS(kxO(j) - kxO(i))/kxO(i) <= ktol)) THEN
				! IF (((ABS(kyO(j) - kyO(i))/kyO(i) <= ktol) .AND. (ABS(kzO(j)-kzO(i))/kzO(i) <= ktol)) .AND. &
						! (ABS(kxO(j) - kxO(i))/kxO(i) <= ktol)) THEN
					numElDegen(m) = numElDegen(m) + 1;
					indToRemove(SUM(numElDegen)) = j;
					weightsO(i) = weightsO(i) + weightsO(j);
					! CO(i) = CO(i) + CO(j);		! heat capacity accounts for degeneracy.  
						!in getCosExpCoeffs, we multiply C by weight, so weighting will be double counted...
					weightsO(j) = 0;					
				END IF
			END DO
		END IF
	END DO
END DO

numEachFreq(1:numUniqueFreq) = numEachFreq(1:numUniqueFreq) - numElDegen;  ! since we're going to remove this kpoint from the list

i = 1;
mm = 1;
DO m = 1,wTot
	! indToRemove is not necessarily in increasing order!
	IF (NOT(ANY(indToRemove==m))) THEN          
		kx(mm) = kxO(m); ky(mm) = kyO(m); kz(mm) = kzO(m); 
		dkx(mm) = dkxO(m); dky(mm) = dkyO(m); dkz(mm) = dkzO(m);
		! should the tau's, vg's  be averages?
		freq(mm) = freqO(m); tau(mm) = tauO(m);
		vgx(mm) = vgxO(m); vgy(mm) = vgyO(m); vgz(mm) = vgzO(m); 
		C(mm) = CO(m); MFP(mm) = MFPO(m); 
		weights(mm) = weightsO(m); gamma(mm) = gammaO(m);
		mm = mm + 1;
	ELSE
		! PRINT *, m, indToRemove(i)
		i = i+1
	END IF
	
END DO


PRINT *, 'numdegen:', SUM(numElDegen)
PRINT *, 'new wTot:', mm-1, 'oldwTot:', w
PRINT *, SUM(C(1:wTot)), SUM(vgx(1:wTot)), SUM(tau(1:wTot)), SUM(weights(1:wTot))

kx(mm:wTot) = 0; ky(mm:wTot) = 0; kz(mm:wTot) = 0;
dkx(mm:wTot) = 0; dky(mm:wTot) = 0; dkz(mm:wTot) = 0;
freq(mm:wTot) = 0; tau(mm:wTot) = 0; kz(mm:wTot) = 0;
vgx(mm:wTot) = 0; vgy(mm:wTot) = 0; vgz(mm:wTot) = 0;
C(mm:wTot) = 0;  MFP(mm:wTot) = 0; weights(mm:wTot) = 0; gamma(mm:wTot) = 0;

wTot = mm-1;

DEALLOCATE(kxO,kyO,kzO,dkxO,dkyO,dkzO,vgxO,vgyO,vgzO,tauO,freqO,CO,MFPO,weightsO,gammaO)
DEALLOCATE(kxF,kyF,kzF,dkxF,dkyF,dkzF,vgxF,vgyF,vgzF,tauF,freqF,CF,MFPF,weightsF,gammaF)
DEALLOCATE(indToRemove, numElDegen, indList)
DEALLOCATE(kxTemp,kxUnique,kyTemp,kyUnique,kzTemp,kzUnique,kxInd,kyInd,kzInd,&
			kxIndAll,kyIndAll,kzIndAll,dkxU,dkyU,dkzU,dkxIndAll,dkyIndAll,dkzIndAll)
			
CONTAINS

	FUNCTION PARTITION (A, B, lo, hi) RESULT (breakpoint)		! unstable bc of repeated elements?
		IMPLICIT NONE
		
		REAL, DIMENSION(:) :: A  !! list of #s that we are comparing
		INTEGER, DIMENSION(:) :: B	!! list of indices corresponding to elements in B
		INTEGER :: lo, hi, breakpoint
		INTEGER :: temp2, i, j
		REAL :: piv, temp
		
		piv = A((lo+hi)/2);		! make to be some random value between lo, hi
		i = lo;
		j = hi;
		
		DO 			!! goes on forever
			DO WHILE (A(i) < piv)   ! pivot is larger element--move to +1 element in A
				i = i + 1
			END DO
			
			DO WHILE (A(j) > piv)	! pivot is smaller element--move to -1 element in A
				j = j - 1
			END DO
			
			IF (i >= j) THEN 		! swap point reached, or all elements are equal
				breakpoint = j
				EXIT
				! what happens when A(i) = A(j) = piv?
			ELSE IF ((A(i) == piv).AND.(A(j) == piv)) THEN
				i = i + 1
			ELSE
				temp = A(i)
				A(i) = A(j)
				A(j) = temp
				
				temp2 = B(i)
				B(i) = B(j)
				B(j) = temp2			
			END IF
		END DO
		
	END FUNCTION
	
	FUNCTION PARTITION_lomuto (A, B, lo, hi) RESULT (breakpoint)
		IMPLICIT NONE
		
		REAL, DIMENSION(:) :: A  !! list of #s that we are comparing
		INTEGER, DIMENSION(:) :: B	!! list of indices corresponding to elements in B
		INTEGER :: lo, hi, breakpoint
		INTEGER :: temp2, i, j
		REAL :: piv, temp
		
		piv = A(hi);		! make to be some random value between lo, hi
		i = lo;
		j = hi;
		
		DO j = lo, hi-1 			!! goes on forever
			IF (A(j) <= piv) THEN
				i = i + 1;
				temp = A(i)
				A(i) = A(j)
				A(j) = temp
				
				temp2 = B(i)
				B(i) = B(j)
				B(j) = temp2
			END IF
		END DO
		
		temp = A(i+1)
		A(i+1) = A(j)
		A(j) = temp
		
		temp2 = B(i+1)
		B(i+1) = B(j)
		B(j) = temp2
		
		breakpoint = i + 1;
		
	END FUNCTION
	
	RECURSIVE SUBROUTINE QUICKSORT(A, B, lo, hi)
		IMPLICIT NONE
		
		REAL, DIMENSION(:) :: A
		INTEGER, DIMENSION(:) :: B
		INTEGER :: lo, hi
		INTEGER :: bkpt
		
		IF (lo < hi) THEN
			bkpt = PARTITION(A, B, lo, hi)
			! PRINT *, bkpt, lo, hi
			CALL QUICKSORT(A(lo:bkpt), B(lo:bkpt), 1, (bkpt-lo))
			CALL QUICKSORT(A(bkpt+1:hi), B(bkpt+1:hi), 1, hi-bkpt)
		ELSE
			! PRINT *, 'end recursion?'
			A = A; B = B;
		END IF
	END SUBROUTINE

END SUBROUTINE orderParams