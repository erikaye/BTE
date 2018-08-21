FUNCTION GetMatlParams(DFT_path) RESULT(allMatls)

USE Constants
IMPLICIT NONE

! REAL, INTENT(IN) :: roughness
INTEGER :: m, ierr, numModes, numPts, wTot, wMax
CHARACTER(longStr), INTENT(IN) :: DFT_path
CHARACTER(longStr) :: fileName, fileLoc
CHARACTER(medStr) :: purity
CHARACTER(shortStr) :: matl, tempstr
TYPE(material), DIMENSION(numMatls) :: allMatls

REAL, DIMENSION(:), ALLOCATABLE :: kx, ky, kz
REAL, DIMENSION(:,:), ALLOCATABLE :: vgx, vgy, vgz, tau, freq, weights

INTEGER :: numUniqueFreq
INTEGER, DIMENSION(maxPts*maxModes) :: numEachFreq
REAL, DIMENSION(maxPts*maxModes) :: uniqueFreqs
REAL :: rho

INTERFACE
	SUBROUTINE orderParams(kx,ky,kz,dkx,dky,dkz,vgx,vgy,vgz,tau,freq,C,MFP,weights,gamma, &
							numUniqueFreq, numEachFreq, uniqueFreqs, numPts, wTot, df, dataset)
		USE Constants
		IMPLICIT NONE

		REAL, DIMENSION(:) :: kx,ky,kz,dkx,dky,dkz,vgx,vgy,vgz,tau,freq,C,MFP,weights
		REAL*8, DIMENSION(:) :: gamma
		
		INTEGER, INTENT(IN) :: numPts
		CHARACTER(4), INTENT(IN) :: dataset
		INTEGER :: wTot
		INTEGER, INTENT(OUT) :: numUniqueFreq, numEachFreq(:)
		REAL, INTENT(OUT) :: uniqueFreqs(:), df
	END SUBROUTINE
END INTERFACE


WRITE(tempstr,'(I4)') T 

ALLOCATE(kx(maxPts),ky(maxPts),kz(maxPts),vgx(maxPts,maxModes),vgy(maxPts,maxModes),&
		 vgz(maxPts,maxModes), tau(maxPts,maxModes),freq(maxPts,maxModes),weights(maxPts,maxModes))

wMax = maxPts*maxModes;
		 
DO m = 1,numMatls
    matl = TRIM(matlist(m))
    WRITE(*,*), matl
    IF ( matl == 'Al') THEN
        numModes = 3
        purity = 'natural'
		rho = 2830 ! kg/m^3
        allMatls(m) = getAlParams()
    ELSE IF ((matl=='Si').OR.(matl=='SiX').OR.(matl=='MgO').OR.(matl=='Ge').OR.(matl=='AlP')) THEN
        numModes = 6
        purity = 'natural'
		IF ((matl=='Si').OR.(matl=='SiX')) THEN 
			rho = 2328 
		ELSE IF (matl=='MgO') THEN 
			rho = 3600 
		ELSE IF (matl=='Ge') THEN 
			rho = 5323 
		ELSE	! AlP
			rho = 2850 
		END IF
        allMatls(m) = getLindsayParams()
    ELSE
        IF (matl == 'SiGe') THEN
            purity = 'natural'
            ! Si can be natural or pure
        ELSE
            purity = 'pure'
        END IF
        numModes = 6
        allMatls(m) = getMingoParams()
		rho = 4000;  ! place holder...
    END IF
    ! Si can be from Lindsay or Mingo
    
END DO 

PRINT *, 'got material params'

DEALLOCATE(kx,ky,kz,vgx,vgy,vgz,tau,freq,weights)

! Lindsay:  DFT data is in the positive kx,ky,kz octant of BZ.  
! Mingo and Al:  DFT data is in the whole BZ.  Crystal momentum vector can be added to negative values to bring 
! find k-pts into the positive octant of BZ.  (k-pts are not evenly distributed in BZ)

CONTAINS 
    FUNCTION getAlParams() RESULT(matData)
        IMPLICIT NONE
        
        INTEGER :: i, i2, notIncl
        ! single line of data from the file
        REAL, DIMENSION(:), ALLOCATABLE :: kxT, kyT, kzT     ! T for temporary
        REAL, DIMENSION(:,:), ALLOCATABLE :: vgxT, vgyT, vgzT, freqT
		REAL, DIMENSION(:), ALLOCATABLE :: kxCol, kyCol, kzCol, dkx, dky, dkz     ! ordered column vectors
        REAL, DIMENSION(:), ALLOCATABLE :: vgxCol, vgyCol, vgzCol, freqCol, tauCol, weightsCol,&
								C_mu, lambda_mu, spec
		REAL*8, DIMENSION(:), ALLOCATABLE :: gamma_mu, E_mu, gBE
		REAL :: df, kBulk
        
        TYPE(material) :: matData
        
        ! matrix inversion to find unit cell volume
        INTEGER :: loc_x, loc_y, loc_z
        REAL, DIMENSION(3,3) :: rLattVecs, LattVecs
        REAL :: vol
        
		
        ! re-initiate
        kx = 0; ky = 0; kz = 0; vgx = 0; vgy = 0; vgz = 0; tau = 0; freq = 0;
		
		ALLOCATE(kxT(maxPts),kyT(maxPts),kzT(maxPts),vgxT(maxPts,numModes),vgyT(maxPts,numModes),&
				 vgzT(maxPts,numModes),freqT(maxPts,numModes))
		
        fileName = TRIM(purity)//'_'//TRIM(ADJUSTL(tempstr))//'K_'//TRIM(matl)//'.txt'
        !fileLoc = './DFT_Data/Al/'//fileName
		! fileLoc = './../DFT_Data/Al/'//fileName
		fileLoc = TRIM(DFT_path)//'DFT_Data/Al/'//fileName
        ! WRITE(*,*) fileLoc

        OPEN (1,FILE=fileLoc)
        
        READ(1,*)
        ! Skip the first line
        
        ReadLoop: DO i = 1,maxLines  !info for modes in different column
            READ (1,*,IOSTAT=ierr) kxT(i), kyT(i), kzT(i), freqT(i,1), freqT(i,2), freqT(i,3), &
                                    vgxT(i,1), vgyT(i,1), vgzT(i,1), &
                                    vgxT(i,2), vgyT(i,2), vgzT(i,2), &
                                    vgxT(i,3), vgyT(i,3), vgzT(i,3)
            
            
            IF ( ierr /= 0 ) THEN
                numPts = i - 1
                IF ( ierr == -1 ) THEN
                    WRITE(*,*) 'end of file'
                    WRITE(*,*) numPts
                    EXIT ReadLoop
                ELSE
                    WRITE(*,'(/ "Error on read: ", I0 )' ) ierr
                    WRITE(*,*) numPts
                    STOP
                END IF
            END IF
        END DO ReadLoop
        
        IF (ierr .EQ. 0) THEN
            numPts = i
            WRITE(*,*) 'Increase maxLines'
        END IF
        
        CLOSE(1)
		
        i2 = 0;
		notIncl = 0;
        DO i=1,numPts
            ! IF ((kxT(i) >= 0) .AND. (kyT(i) >= 0) .AND. (kzT(i) >= 0)) THEN
			! IF ((vgxT(i,1) > 0) .AND. (SQRT(vgxT(i,1)**2+vgyT(i,1)**2+vgzT(i,1)) > 0) & 
					! .AND. (SQRT(vgxT(i,2)**2+vgyT(i,2)**2+vgzT(i,2)) > 0) & 
					! .AND. (SQRT(vgxT(i,3)**2+vgyT(i,3)**2+vgzT(i,3)) > 0)) THEN    
			IF (vgxT(i,1) >= 0) THEN
				! select only forward propagating mode 1
                i2 = i2+1;
        
				kx(i2) = kxT(i); ky(i2) = kyT(i); kz(i2) = kzT(i);
                vgx(i2,1:numModes) = vgxT(i,:); vgy(i2,1:numModes) = vgyT(i,:); 
				! not all velocities for all modes = 0... 
                vgz(i2,1:numModes) = vgzT(i,:); freq(i2,1:numModes) = freqT(i,:); 
			ELSE
				notIncl = notIncl + 1;
            END IF
        END DO
        
        numPts = i2;   ! redefine numPts for BZ octant
        
		wTot = numPts*numModes;
		PRINT *, 'Al num modes', numModes
		PRINT *, numPts, notIncl
		
		ALLOCATE(kxCol(wMax),kyCol(wMax),kzCol(wMax), dkx(wMax),dky(wMax),dkz(wMax), &
				 vgxCol(wMax),vgyCol(wMax),vgzCol(wMax), &
				 freqCol(wMax),tauCol(wMax),weightsCol(wMax),&
				 E_mu(wMax),gBE(wMax),C_mu(wMax),lambda_mu(wMax),gamma_mu(wMax),spec(wMax))
				 
		kxCol(1:numPts) = kx(1:numPts); kyCol(1:numPts) = ky(1:numPts); kzCol(1:numPts) = kz(1:numPts)
		vgxCol(1:wTot) = RESHAPE(vgx(1:numPts,1:numModes),(/wTot/))
		vgyCol(1:wTot) = RESHAPE(vgy(1:numPts,1:numModes),(/wTot/))
		vgzCol(1:wTot) = RESHAPE(vgz(1:numPts,1:numModes),(/wTot/))
		vgxCol = vgxCol + 1.0E-12; vgyCol = vgyCol + 1.0E-12; vgzCol = vgzCol + 1.0E-12;  ! offset to avoid 0s
		freqCol(1:wTot) = RESHAPE(freq(1:numPts,1:numModes),(/wTot/))
		
		tauCol = 70.0E-9 / SQRT(vgxCol(1:wTot)**2 + vgyCol(1:wTot)**2 + vgzCol(1:wTot)**2)
        freqCol = freqCol*1.0E12*2.0*pi
		PRINT *, MINVAL(vgxCol(1:wTot)), MINVAL(SQRT(vgxCol(1:wTot)**2 + vgyCol(1:wTot)**2 + vgzCol(1:wTot)**2))
		! can be negative, can be 0...
        
        E_mu = hbar*freqCol;                             ! energy of each mode.
        gBE = E_mu*(1.0/(EXP(E_mu/kB/T)-1));            ! energy*number density of each mode 
        C_mu = 1.0/kB/(T**2)*gBE*(gBE+hbar*freqCol);     ! heat capacity for each mode
		
        lambda_mu = vgxCol*tauCol                          ! mean free path
        ! gamma_mu = CMPLX(1/tauCol,eta)/vgxCol
		gamma_mu = 1.0/lambda_mu
        
        ! find volume of unit cell
        loc_x = MINLOC(kxT,1); loc_y = MINLOC(kyT,1); loc_z = MINLOC(kzT,1);
        ! reciprocal lattice vectors
        rLattVecs(:,1) = (/ kxT(loc_x), kyT(loc_x), kzT(loc_x) /)
        rLattVecs(:,2) = (/ kxT(loc_y), kyT(loc_y), kzT(loc_y) /)
        rLattVecs(:,3) = (/ kxT(loc_z), kyT(loc_z), kzT(loc_z) /)
        ! real space lattice vectors
        LattVecs = 2*pi*sINV(rLattVecs);
        vol = sDET3(LattVecs);
        PRINT*, 'Unit Vol', vol
        weightsCol = 1.0/numPts/vol;
		weightsCol = weightsCol/((numPts+notIncl)/numPts);
			! in our calc, we separate forward propagating and backward propagating, so only counted ~half # kpts in BZ 
			! (approx... bc this data is weird)  only considered numPts/(numPts+notIncl)~2 of BZ
			! sum(weights) over entire BZ = 1/vol
			
		CALL orderParams(kxCol(1:wTot),kyCol(1:wTot),kzCol(1:wTot),&
						dkx(1:wTot),dky(1:wTot),dkz(1:wTot),&
						vgxCol(1:wTot),vgyCol(1:wTot),vgzCol(1:wTot),&
						tauCol(1:wTot),freqCol(1:wTot),C_mu(1:wTot),&
						lambda_mu(1:wTot),weightsCol(1:wTot),gamma_mu(1:wTot), &
						numUniqueFreq, numEachFreq, uniqueFreqs,numPts, wTot, df, 'Al  ')
		
		! spec(1:wTot) = EXP(-4.0*(roughness**2)*kxCol(1:wTot)**2);
		spec(1:wTot) = EXP(-4.0*(kxCol(1:wTot)/MAXVAL(kxCol(1:wTot)))**2);
		
		! df= MINVAL(vgxCol(1:wTot)*dkx(1:wTot) + vgyCol(1:wTot)*dky(1:wTot) + vgzCol(1:wTot)*dkz(1:wTot))
		PRINT *, 'df', df
		
		kBulk = SUM(C_mu(1:wTot)*(vgxCol(1:wTot)**2)*tauCol(1:wTot) * (((numPts+notIncl)/numPts)*weightsCol(1:wTot)) )
		! kBulk = 1.0/3.0*SUM(C_mu(1:wTot)*(vgxCol(1:wTot)**2 + vgyCol(1:wTot)**2 + vgzCol(1:wTot)**2)*tauCol(1:wTot) * &
						! (((numPts+notIncl)/numPts)*weightsCol(1:wTot)) )
			! isotropic
			! *((numPts+notIncl)/numPts) ~ 2 because need to sum over entire BZ
		PRINT *, 'k_bulk', kBulk
		
        matData = material(matl,fileLoc,numModes,wTot,&
                                kxCol,kyCol,kzCol,dkx,dky,dkz,vgxCol,vgyCol,vgzCol,tauCol,freqCol,&
                                C_mu,lambda_mu,weightsCol,spec,gamma_mu,&
								numUniqueFreq,numEachFreq,uniqueFreqs,df,rho)
	
	
		DEALLOCATE(kxT,kyT,kzT,vgxT,vgyT,vgzT,freqT)
		DEALLOCATE(kxCol,kyCol,kzCol,dkx,dky,dkz,vgxCol,vgyCol,vgzCol,freqCol,tauCol,weightsCol,&
				 E_mu,gBE,C_mu,lambda_mu,gamma_mu,spec)
	
    END FUNCTION getAlParams

    
    FUNCTION getLindsayParams() RESULT (matData)
        IMPLICIT NONE
        
        INTEGER :: i, numWeights=0, modenumber
        INTEGER, DIMENSION(numModes) :: kpt
        REAL :: dummy, kBulk
		
		REAL, DIMENSION(:), ALLOCATABLE :: kxCol, kyCol, kzCol, dkx, dky, dkz   ! ordered column vectors
        REAL, DIMENSION(:), ALLOCATABLE :: vgxCol, vgyCol, vgzCol, freqCol, tauCol, weightsCol, &
								C_mu, lambda_mu, spec
		REAL*8, DIMENSION(:), ALLOCATABLE :: gamma_mu, E_mu, gBE
		REAL :: df
        
        TYPE(material) :: matData
        
        CHARACTER(longStr) :: fileLocW
        CHARACTER(shortStr) :: matlReal
        
        
        ! re-initiate
        kx = 0; ky = 0; kz = 0; vgx = 0; vgy = 0; vgz = 0; tau = 0; freq = 0;
        
        matlReal = matl
        IF (matl=='SiX') THEN
            matlReal = 'Si'
        END IF
        
        fileName = TRIM(purity)//'_'//TRIM(ADJUSTL(tempstr))//'K_'//TRIM(matlReal)//'.txt'
        ! fileLoc = './DFT_Data/Lindsay/'//TRIM(matl)//'/'//fileName
		! fileLoc = './../DFT_Data/Lindsay/'//TRIM(matl)//'/'//fileName
		fileLoc = TRIM(DFT_path)//'DFT_Data/Lindsay/'//TRIM(matlReal)//'/'//fileName
        ! WRITE(*,*) fileLoc
        
        OPEN (1,FILE=fileLoc)
        
        ! data is presented mode by mode, so need to count kpts per mode
        kpt = 1;  !all elements in kpt = 1
        
        ReadLoop: DO i = 1,maxLines  !info for all modes in one column
            READ (1,*,IOSTAT=ierr) modenumber, kx(kpt(modenumber)), ky(kpt(modenumber)), & 
                                    kz(kpt(modenumber)), freq(kpt(modenumber),modenumber), &
                                    vgx(kpt(modenumber),modenumber), vgy(kpt(modenumber),modenumber), &
                                    vgz(kpt(modenumber),modenumber),dummy, tau(kpt(modenumber),modenumber)
            
            IF ( ierr /= 0 ) THEN
                numPts = SUM(kpt-1)/numModes  ! kpt(modenumber) gets extra 1 when it moves to next modenumber
                IF ( ierr == -1 ) THEN   ! read error -> ierr = -1
                    WRITE(*,*) 'end of file'
                    WRITE(*,*) ierr
                    WRITE(*,*) numPts
                    EXIT ReadLoop
                ELSE
                    WRITE(*,'(/ "Error on read: ", I0 )' ) ierr
                    STOP
                END IF
            END IF
    
            kpt(modenumber) = kpt(modenumber) + 1;  ! increase index for kpt for given mode
        END DO ReadLoop
        
        
        IF (ierr .EQ. 0) THEN
            numPts = SUM(kpt)/numModes
            WRITE(*,*) 'Increase maxLines'
        END IF
        
        CLOSE(1)

		wTot = numPts*numModes;
		
        ! read integration weights file
        fileName = TRIM(matl)//'_Integration_Weights.txt'
        ! fileLocW = './DFT_Data/Lindsay/'//TRIM(matl)//'/'//fileName
		! fileLocW = './../DFT_Data/Lindsay/'//TRIM(matl)//'/'//fileName
		fileLocW = TRIM(DFT_path)//'DFT_Data/Lindsay/'//TRIM(matlReal)//'/'//fileName
        ! WRITE(*,*) fileLocW

        OPEN (2,FILE=fileLocW)
                
        ReadLoopW: DO i = 1,numPts+1 !info for diff modes in diff columns
            READ (2,*,IOSTAT=ierr) weights(i,:)
            
            IF ( ierr /= 0 ) THEN
                numWeights = i-1
                IF ( ierr == -1 .OR. ierr > 0 ) THEN   ! read error -> ierr = -1
                    WRITE(*,*) 'end of file'
                    WRITE(*,*) numWeights
                    EXIT ReadLoopW
                ELSE
                    WRITE(*,*) numWeights
                    WRITE(*,'(/ "Error on read: ", I0 )' ) ierr
                    STOP
                END IF
            END IF
        END DO ReadLoopW
        
        IF (ierr .EQ. 0) THEN
            numWeights = i
            WRITE(*,*) 'Increase maxLines'
        END IF
        
        CLOSE(2)
        
        IF (numWeights /= numPts) THEN
            WRITE(*,*) 'error: numweights != numPts'
            WRITE(*,*) numWeights, numPts
        END IF
        

        kx = kx/1.0E-10; ky = ky/1.0E-10; kz = kz/1.0E-10
        freq = freq*1.0E12*2.0*pi
        tau = tau*1.0E-12
		weights = weights*1.0E30;   ! weights are in 1/Angstroms.  sum(weights) = 1/Vol unit cell
        
        if (matl=='SiX') THEN
            tau = tau/100.0
            PRINT *, 'SiX tau: tau/100'
        END IF
        
		ALLOCATE(kxCol(wMax),kyCol(wMax),kzCol(wMax), dkx(wMax),dky(wMax),dkz(wMax), &
				 vgxCol(wMax),vgyCol(wMax),vgzCol(wMax), &
				 freqCol(wMax),tauCol(wMax),weightsCol(wMax),&
				 E_mu(wMax),gBE(wMax),C_mu(wMax),lambda_mu(wMax),gamma_mu(wMax),spec(wMax))
		
		kxCol(1:numPts) = kx(1:numPts); kyCol(1:numPts) = ky(1:numPts); kzCol(1:numPts) = kz(1:numPts)
		vgxCol(1:wTot) = ABS(RESHAPE(vgx(1:numPts,1:numModes),(/wTot/)))  ! if sign is flipped, just say we are consider mode with -kx
		vgyCol(1:wTot) = ABS(RESHAPE(vgy(1:numPts,1:numModes),(/wTot/)))
		vgzCol(1:wTot) = ABS(RESHAPE(vgz(1:numPts,1:numModes),(/wTot/)))
		freqCol(1:wTot) = RESHAPE(freq(1:numPts,1:numModes),(/wTot/))
		tauCol(1:wTot) = RESHAPE(tau(1:numPts,1:numModes),(/wTot/))
		weightsCol(1:wTot) = RESHAPE(weights(1:numPts,1:numModes),(/wTot/))
		weightsCol = weightsCol/2.0;  ! in our calc, we separate forward propagating and backward propagating,
									! so our sums are only over half of the BZ
		
        E_mu = hbar*freqCol;                             ! energy of each mode.
        gBE = E_mu*(1.0/(EXP(E_mu/kB/T)-1));            ! energy*number density of each mode
        C_mu = REAL(1.0/kB/(T**2)*gBE*(gBE + E_mu));     ! heat capacity for each mode

        lambda_mu = vgxCol*tauCol                           ! mean free path        
		gamma_mu = 1.0/lambda_mu;
		
		CALL orderParams(kxCol(1:wTot),kyCol(1:wTot),kzCol(1:wTot),&
						dkx(1:wTot),dky(1:wTot),dkz(1:wTot),&
						vgxCol(1:wTot),vgyCol(1:wTot),vgzCol(1:wTot),&
						tauCol(1:wTot),freqCol(1:wTot),C_mu(1:wTot),&
						lambda_mu(1:wTot),weightsCol(1:wTot),gamma_mu(1:wTot), &
						numUniqueFreq, numEachFreq, uniqueFreqs, numPts, wTot, df, 'Lind')
						
        
		spec(1:wTot) = EXP(-4.0*(roughness**2)*kxCol(1:wTot)**2);
		! ! spec(1:wTot) = EXP(-4.0*(roughness**2)*(kxCol(1:wTot)**2 + kyCol(1:wTot)**2 + kzCol(1:wTot)**2));
		PRINT *, 'ziman spec. roughness:', roughness
		PRINT *, 'min/max spec', MINVAL(spec(1:wTot)), MAXVAL(spec(1:wTot))
		! ! ! spec(1:wTot) = EXP(-4.0*(kxCol(1:wTot)/MAXVAL(kxCol(1:wTot)))**2);
		! ! PRINT *, 'ziman spec. roughness:', 1.0/MAXVAL(kxCol(1:wTot))
		spec(1:wTot) = 0.0
		PRINT *, 'spec 0'
		
		! PRINT *, 'check ordering'
		! DO i = 1,wTot
			! IF (freqCol(i) > freqCol(i+1)) PRINT *, 'out of order at', i-1
			! IF (freqCol(i) == freqCol(i+1)) THEN
				! IF (C_mu(i) /= C_mu(i+1)) PRINT *, 'c_mu error'
			! END IF
		! END DO
		
		kBulk = SUM(C_mu(1:wTot)*(vgxCol(1:wTot)**2)*tauCol(1:wTot) * (2.0*weightsCol(1:wTot)) )
			! *2 bc need to sum over modes in entire BZ
		! kBulk = 1.0/3.0*SUM(C_mu(1:wTot)*(vgxCol(1:wTot)**2 + vgyCol(1:wTot)**2 + vgzCol(1:wTot)**2)*tauCol(1:wTot) * &
					! (2.0*weightsCol(1:wTot)) )
			! isotropic
			
		PRINT *, 'k_bulk', kBulk
		
        matData = material(matl,fileLoc,numModes,wTot,&
                                kxCol,kyCol,kzCol,dkx,dky,dkz,vgxCol,vgyCol,vgzCol,tauCol,freqCol,&
                                C_mu,lambda_mu,weightsCol,spec,gamma_mu,&
								numUniqueFreq,numEachFreq,uniqueFreqs,df,rho)
		
		
		DEALLOCATE(kxCol,kyCol,kzCol,dkx,dky,dkz,vgxCol,vgyCol,vgzCol,freqCol,tauCol,weightsCol,&
				 E_mu,gBE,C_mu,lambda_mu,gamma_mu,spec)
    
		
	
    END FUNCTION getLindsayParams

    
    FUNCTION getMingoParams() RESULT(matData)
        IMPLICIT NONE
        
        INTEGER :: i, i2, numPts, modenumber, kpt, notIncl
        REAL :: dummy

		REAL, DIMENSION(:), ALLOCATABLE :: kxT, kyT, kzT     ! T for temporary
        REAL, DIMENSION(:,:), ALLOCATABLE :: vgxT, vtotT, freqT, tauT
		REAL, DIMENSION(:), ALLOCATABLE :: kxCol, kyCol, kzCol, dkx, dky, dkz     ! ordered column vectors
        REAL, DIMENSION(:), ALLOCATABLE :: vgxCol, vgyCol, vgzCol, freqCol, tauCol, weightsCol, &
								C_mu, lambda_mu, spec
		REAL*8, DIMENSION(:), ALLOCATABLE :: gamma_mu, E_mu, gBE
		REAL :: df, kBulk
        
        REAL, DIMENSION(maxPts,maxModes) :: vtot

        TYPE(material) :: matData
        
        ! matrix inversion to find unit cell volume
        INTEGER :: loc_x, loc_y, loc_z
        REAL, DIMENSION(3,3) :: rLattVecs, LattVecs
        REAL :: vol
        
        
        ! re-initiate
        kx = 0; ky = 0; kz = 0; vgx = 0; vgy = 0; vgz = 0; tau = 0; freq = 0;
        
		ALLOCATE(kxT(maxPts),kyT(maxPts),kzT(maxPts),vgxT(maxPts,numModes),vtotT(maxPts,numModes),&
				 freqT(maxPts,numModes),tauT(maxPts,numModes))
		
		kxT = 0; kyT = 0; kzT = 0; vgxT = 0; vtotT = 0; freqT = 0; tauT = 0;
		
        fileName = TRIM(purity)//'_'//TRIM(ADJUSTL(tempstr))//'K_'//TRIM(matl)//'.txt'
		! fileLoc = './../DFT_Data/Mingo/'//TRIM(matl)//'/'//fileName
		fileLoc = TRIM(DFT_path)//'DFT_Data/Mingo/'//TRIM(matl)//'/'//fileName
        
        OPEN (1,FILE=fileLoc)

        ! skip the first 6 lines
        DO i = 1,6
            READ(1,*)
        END DO
        
        ReadLoop: DO i = 1,maxLines  !info for all modes in one column
            ! data for all modes is given for single kpt at once
            kpt = 1 + (i-1)/numModes
            READ (1,*,IOSTAT=ierr) kxT(kpt), kyT(kpt), kzT(kpt), modenumber, freqT(kpt,modenumber), &
                                    vgxT(kpt,modenumber), dummy, vtotT(kpt,modenumber), tauT(kpt,modenumber)
            
            
            IF ( ierr /= 0 ) THEN
                numPts = kpt - 1  
                IF ( ierr == -1 ) THEN   ! read error -> ierr = -1                   
                    ! WRITE(*,*) 'end of file'
                    ! WRITE(*,*) numPts
                    EXIT ReadLoop
                ELSE
                    WRITE(*,'(/ "Error on read: ", I0 )' ) ierr
                    STOP
                END IF
            END IF
        END DO ReadLoop
        
        IF (ierr .EQ. 0) THEN
            numPts = kpt
            WRITE(*,*) 'Increase maxpoints'
        END IF
        
        CLOSE(1)
        		
        i2 = 0;
		notIncl = 0;
        DO i=1,numPts
            ! IF ((kxT(i) >= 0) .AND. (kyT(i) >= 0) .AND. (kzT(i) >= 0)) THEN
			IF (vgxT(i,1) >= 0) THEN    ! select only forward propagating mode 1
										! by symmetry, all forward propagating modes will have backward propagating counterpart
										! this set of data measured entire BZ
				! need to fix order params if choose data points this way
                i2 = i2+1;
                
                kx(i2) = kxT(i); ky(i2) = kyT(i); kz(i2) = kzT(i);
                vgx(i2,:) = ABS(vgxT(i,:)); vtot(i2,:) = vtotT(i,:);
                tau(i2,:) = tauT(i,:); freq(i2,:) = freqT(i,:); 
                
                vgy(i2,:) = SQRT(vtot(i2,:)**2 - vgx(i2,:)**2) * (ky(i2)/SQRT(ky(i2)**2 + kz(i2)**2))
                vgz(i2,:) = SQRT(vtot(i2,:)**2 - vgx(i2,:)**2) * (kz(i2)/SQRT(ky(i2)**2 + kz(i2)**2))
            ELSE
				notIncl = notIncl + 1;
            END IF
        END DO
        
        numPts = i2;   ! redefine numPts for BZ octant
        wTot = numPts*numModes;
		PRINT *, numPts, notIncl
		
		kx = kx*1.0E9; ky = ky*1.0E9; kz = kz*1.0E9
        vgx = vgx*1000.0; vgy = vgy*1000.0; vgz = vgz*1000.0
        ! freq = freq*0.24*1.0E12*2.0*pi;  ! originally in MeV
		freq = freq*1.602E-22/hbar;  ! originally in MeV -> J -> rad/s
		tau = tau*1.0E-12
		
		ALLOCATE(kxCol(wMax),kyCol(wMax),kzCol(wMax), dkx(wMax),dky(wMax),dkz(wMax), &
				 vgxCol(wMax),vgyCol(wMax),vgzCol(wMax), &
				 freqCol(wMax),tauCol(wMax),weightsCol(wMax),&
				 E_mu(wMax),gBE(wMax),C_mu(wMax),lambda_mu(wMax),gamma_mu(wMax),spec(wMax))
		
		kxCol(1:numPts) = kx(1:numPts); kyCol(1:numPts) = ky(1:numPts); kzCol(1:numPts) = kz(1:numPts)
		vgxCol(1:wTot) = RESHAPE(vgx(1:numPts,1:numModes),(/wTot/))
		vgyCol(1:wTot) = RESHAPE(vgy(1:numPts,1:numModes),(/wTot/))
		vgzCol(1:wTot) = RESHAPE(vgz(1:numPts,1:numModes),(/wTot/))
		vgxCol = vgxCol + 1.0E-12; vgyCol = vgyCol + 1.0E-12; vgzCol = vgzCol + 1.0E-12;
		freqCol(1:wTot) = RESHAPE(freq(1:numPts,1:numModes),(/wTot/))
		tauCol(1:wTot) = RESHAPE(tau(1:numPts,1:numModes),(/wTot/))
		
        E_mu = hbar*freqCol;                             ! energy of each mode.
        gBE = E_mu*(1.0/(EXP(E_mu/kB/T)-1));            ! energy*number density of each mode 
        C_mu = 1.0/kB/(T**2)*gBE*(gBE+E_mu);     		! heat capacity for each mode

        lambda_mu = vgxCol*tauCol                           ! mean free path
		gamma_mu = 1.0/lambda_mu

        ! find volume of unit cell
        loc_x = MINLOC(kxT,1); loc_y = MINLOC(kyT,1); loc_z = MINLOC(kzT,1);
		PRINT *, MINVAL(kxT), MINVAL(ABS(kxT)), MINLOC(kxT), MINLOC(ABS(kxT))
        ! reciprocal lattice vectors
        rLattVecs(:,1) = (/ kxT(loc_x), kyT(loc_x), kzT(loc_x) /) * 1.0E9
        rLattVecs(:,2) = (/ kxT(loc_y), kyT(loc_y), kzT(loc_y) /) * 1.0E9
        rLattVecs(:,3) = (/ kxT(loc_z), kyT(loc_z), kzT(loc_z) /) * 1.0E9
        ! real space lattice vectors
        LattVecs = 2.0*pi*sINV(rLattVecs);
        vol = sDET3(LattVecs);
        PRINT *, 'Unit Vol', vol

        weightsCol(1:wTot) = 1.0/numPts/vol;  !kpts*vol = total vol... kpts = numPts, + 1 for the first kpt we skipped?
		weightsCol = weightsCol/((numPts+notIncl)/numPts)  ! sum(weights) over entire BZ = 1/vol
			! in our calc, we separate forward propagating and backward propagating, so only counted ~half of the BZ
			! kpt = numPts = num kpts for forward propagating phonons
		
		DEALLOCATE(kxT,kyT,kzT,vgxT,vtotT,freqT,tauT)
		
		PRINT *, SUM(C_mu(1:wTot)), SUM(vgxCol(1:wTot)), SUM(tauCol(1:wTot)), SUM(weightsCol(1:wTot))
		
		CALL orderParams(kxCol(1:wTot),kyCol(1:wTot),kzCol(1:wTot),&
						dkx(1:wTot),dky(1:wTot),dkz(1:wTot),&
						vgxCol(1:wTot),vgyCol(1:wTot),vgzCol(1:wTot),&
						tauCol(1:wTot),freqCol(1:wTot),C_mu(1:wTot),&
						lambda_mu(1:wTot),weightsCol(1:wTot),gamma_mu(1:wTot), &
						numUniqueFreq, numEachFreq, uniqueFreqs, numPts, wTot, df, 'Ming')
		
		! spec(1:wTot) = EXP(-4.0*(roughness**2)*kxCol(1:wTot)**2);
		spec(1:wTot) = EXP(-4.0*(kxCol(1:wTot)/MAXVAL(kxCol(1:wTot)))**2);
		
		! df = MINVAL(vgxCol(1:wTot)*dkx(1:wTot) + vgyCol(1:wTot)*dky(1:wTot) + vgzCol(1:wTot)*dkz(1:wTot))
		PRINT *, 'df', df
		
		kBulk = SUM(C_mu(1:wTot)*(vgxCol(1:wTot)**2)*tauCol(1:wTot) * (((numPts+notIncl)/numPts)*weightsCol(1:wTot)) )
			! *((numPts+notIncl)/numPts) ~ 2 because sum over entire BZ, not just forward propagating modes
		! kBulk = 1.0/3.0*SUM(C_mu(1:wTot)*(vgxCol(1:wTot)**2 + vgyCol(1:wTot)**2 + vgzCol(1:wTot)**2)*tauCol(1:wTot) * &
						!((numPts+notIncl)/numPts)*weightsCol(1:wTot)) )
			! isotropic 
			
		PRINT *, 'k_bulk', kBulk		
		
        matData = material(matl,fileLoc,numModes,wTot,&
                                kxCol,kyCol,kzCol,dkx,dky,dkz,vgxCol,vgyCol,vgzCol,tauCol,freqCol,&
                                C_mu,lambda_mu,weightsCol,spec,gamma_mu,&
								numUniqueFreq,numEachFreq,uniqueFreqs,df,rho)
    
		DEALLOCATE(kxCol,kyCol,kzCol,dkx,dky,dkz,vgxCol,vgyCol,vgzCol,freqCol,tauCol,weightsCol,&
				 E_mu,gBE,C_mu,lambda_mu,gamma_mu,spec)
    
    END FUNCTION getMingoParams
    
END FUNCTION
