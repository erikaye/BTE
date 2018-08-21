SUBROUTINE GetInterfaceBM_vSpec(matlParams,Ls,Qs,skD,bc, spec1,spec2, &
							G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
							MatABCD, MatEF, MatGH, VecJK, Amat1, f11Coeffs, f12Coeffs, f13Coeffs, Amat2, f21Coeffs, f22Coeffs, f23Coeffs)
! Get interface block matrices (works, copied from _v6 file)

USE omp_lib
USE mkl_service
USE Constants
IMPLICIT NONE

TYPE(material), DIMENSION(2), INTENT(IN) :: matlParams		
			! material parameters from DFT for 2 adjacent blocks			
REAL, DIMENSION(2), INTENT(IN) :: Qs, Ls
            ! heat source magnitude from laser at surface of slab, length of slabs
			! for the two slabs of interest
CHARACTER(2), INTENT(IN) :: bc		! which layer has specular bc's  (s=specular reflection,b=blackbody)
			! (not sure if derivation is independent of direction, or only specular reflection at left boundary)
REAL, INTENT(IN) :: skD, spec1(:,:), spec2(:,:)
			! specularity for each mode
TYPE(sIndVal), DIMENSION(:), INTENT(IN) :: R12_d, T12_d, T21_d, R21_d, R12_s, T12_s, T21_s, R21_s, G12_d, G12_s, G21_d, G21_s
REAL*8, INTENT(OUT) :: MatABCD(:,:), MatEF(:,:), MatGH(:,:), VecJK(:,:), &
						Amat1(:,:), f11Coeffs(:,:), f12Coeffs(:,:), f13Coeffs(:,:), &
						Amat2(:,:), f21Coeffs(:,:), f22Coeffs(:,:), f23Coeffs(:,:)
! REAL*8, ALLOCATABLE :: Amat2(:,:), f21Coeffs(:,:), f22Coeffs(:,:), f23Coeffs(:,:)


INTEGER(KIND=OMP_LOCK_KIND) LCK

INTEGER :: w1, w2, w1d, w2d, wt, statvar, m, loc1, loc2
REAL :: skD1, skD2, rho1, rho2

REAL, DIMENSION(:), ALLOCATABLE :: kx1,kx2,ky1,ky2,kz1,kz2
REAL, DIMENSION(:), ALLOCATABLE :: vgx1, vgy1, vgz1, vgx2, vgy2, vgz2, freq1, freq2, tau1, tau2, &
									C_mu1, C_mu2, Kn_mu1, Kn_mu2, weights1, weights2
REAL*8, DIMENSION(:), ALLOCATABLE :: gamma_mu1, gamma_mu2

INTEGER, ALLOCATABLE :: numEachFreq1(:),numEachFreq2(:)
REAL, ALLOCATABLE :: uniqueFreqs1(:),uniqueFreqs2(:)

REAL*8, DIMENSION(:,:), ALLOCATABLE :: x11Soln, x12Soln, x13Soln, x21Soln, x22Soln, x23Soln, &
										phi1f, phi1b, phi2f, phi2b, phi1fU, phi1bU, phi2fU, phi2bU

REAL*8, DIMENSION(:,:), ALLOCATABLE :: MatABCD_UL, MatEF_UL, MatGH_UL, VecJK_UL, &
										MatABCD_UR, MatEF_UR, MatGH_UR, VecJK_UR, &
										MatABCD_LL, MatEF_LL, MatGH_LL, VecJK_LL, &
										MatABCD_LR, MatEF_LR, MatGH_LR, VecJK_LR
REAL*8, DIMENSION(:,:), ALLOCATABLE :: blockA, blockB, blockC, blockD, blockE, &
										blockF, blockG, blockH, blockJ, blockK, &
										term3_1f, term3_2b
REAL*8, ALLOCATABLE :: temp1(:,:), temp2(:,:), temp3(:)



INTERFACE
	SUBROUTINE GetCosExpCoeffs(Q_mu, skD, L, vgx, vgy, vgz, tau, freq, C_mu, Kn_mu, gamma_mu, weights, &
				Amat,fx1Coeffs,fx2Coeffs,fx3Coeffs,phixb,phixf)
		USE Constants
		IMPLICIT NONE

		REAL, INTENT(IN) :: Q_mu, skD, L			! heat source magnitude from laser at surface of slab, skin depth, layer thickness
		REAL, DIMENSION(:), INTENT(IN) :: vgx, vgy, vgz, tau, freq, C_mu, Kn_mu, weights  ! single column vector containing info for all modes
		REAL*8, DIMENSION(:), INTENT(IN) :: gamma_mu

		REAL*8, DIMENSION(N+1,N+1), INTENT(OUT) :: Amat
		REAL*8, DIMENSION(:,:), INTENT(OUT) :: fx1Coeffs, fx2Coeffs, phixb, phixf
		REAL*8, DIMENSION(N+1), INTENT(OUT) :: fx3Coeffs
	END SUBROUTINE

	SUBROUTINE GetCosExpCoeffs_specBC(Q_mu, skD, L, vgx, vgy, vgz, tau, freq, C_mu, Kn_mu, gamma_mu, weights, &
				Amat,fx1Coeffs,fx2Coeffs,fx3Coeffs,phixb,phixf)
		USE Constants
		IMPLICIT NONE

		REAL, INTENT(IN) :: Q_mu, skD, L			! heat source magnitude from laser at surface of slab, skin depth, layer thickness
		REAL, DIMENSION(:), INTENT(IN) :: vgx, vgy, vgz, tau, freq, C_mu, Kn_mu, weights  ! single column vector containing info for all modes
		REAL*8, DIMENSION(:), INTENT(IN) :: gamma_mu

		REAL*8, DIMENSION(N+1,N+1), INTENT(OUT) :: Amat
		REAL*8, DIMENSION(:,:), INTENT(OUT) :: fx1Coeffs, fx2Coeffs, phixb, phixf
		REAL*8, DIMENSION(N+1), INTENT(OUT) :: fx3Coeffs
	END SUBROUTINE
	
END INTERFACE

! extract material parameters 
skD1 = skD/Ls(1);
skD2 = skD/Ls(2);

w1 = matlParams(1)%wTot
w2 = matlParams(2)%wTot
w1d = matlParams(1)%numUniqueFreq
w2d = matlParams(2)%numUniqueFreq
wt = w1 + w2;


! ALLOCATE(Amat2(N+1,N+1),f21Coeffs(N+1,w2),f22Coeffs(N+1,w2),f23Coeffs(N+1,1))

ALLOCATE(kx1(w1),ky1(w1),kz1(w1),kx2(w2),ky2(w2),kz2(w2))

ALLOCATE(vgx1(w1), vgy1(w1), vgz1(w1), freq1(w1), tau1(w1), C_mu1(w1), Kn_mu1(w1), gamma_mu1(w1), & 
		 weights1(w1), numEachFreq1(w1d), uniqueFreqs1(w1d), phi1f(N+1,w1), phi1b(N+1,w1), phi1fU(N+1,w1), phi1bU(N+1,w1), &
		 x11Soln(N+1,w1), x12Soln(N+1,w1), x13Soln(N+1,1))

ALLOCATE(vgx2(w2), vgy2(w2), vgz2(w2), freq2(w2), tau2(w2), C_mu2(w2), Kn_mu2(w2), gamma_mu2(w2), &
		 weights2(w2), numEachFreq2(w2d), uniqueFreqs2(w2d), phi2f(N+1,w2), phi2b(N+1,w2), phi2fU(N+1,w2), phi2bU(N+1,w2), &
		 x21Soln(N+1,w2), x22Soln(N+1,w2), x23Soln(N+1,1))
		 
phi1f = 0; phi1b = 0; phi2f = 0; phi2b = 0; 
		 
kx1 = matlParams(1)%kx(1:w1);
ky1 = matlParams(1)%ky(1:w1);
kz1 = matlParams(1)%kz(1:w1);
		 
vgx1 = matlParams(1)%vgx(1:w1)
vgy1 = matlParams(1)%vgy(1:w1)
vgz1 = matlParams(1)%vgz(1:w1)
freq1 = matlParams(1)%freq(1:w1)
tau1 = matlParams(1)%tau(1:w1)
rho1 = matlParams(1)%rho
numEachFreq1 = matlParams(1)%numEachFreq(1:w1d)
uniqueFreqs1 = matlParams(1)%uniqueFreqs(1:w1d)

C_mu1 = matlParams(1)%C(1:w1)
Kn_mu1 = matlParams(1)%MFP(1:w1) / Ls(1)
gamma_mu1 = matlParams(1)%gamma(1:w1) * Ls(1)
weights1 = matlParams(1)%weights(1:w1)

kx2 = matlParams(2)%kx(1:w2);
ky2 = matlParams(2)%ky(1:w2);
kz2 = matlParams(2)%kz(1:w2);

vgx2 = matlParams(2)%vgx(1:w2)
vgy2 = matlParams(2)%vgy(1:w2)
vgz2 = matlParams(2)%vgz(1:w2)
freq2 = matlParams(2)%freq(1:w2)
tau2 = matlParams(2)%tau(1:w2)
rho2 = matlParams(2)%rho
numEachFreq2 = matlParams(2)%numEachFreq(1:w2d)
uniqueFreqs2 = matlParams(2)%uniqueFreqs(1:w2d)

C_mu2 = matlParams(2)%C(1:w2)
Kn_mu2 = matlParams(2)%MFP(1:w2) / Ls(2)
gamma_mu2 = matlParams(2)%gamma(1:w2) * Ls(2)
weights2 = matlParams(2)%weights(1:w2)


! cosine expansion coefficients
IF (bc(1).EQ.'b') THEN
	CALL GetCosExpCoeffs(Qs(1), skD, Ls(1), vgx1, vgy1, vgz1, tau1, &
					freq1, C_mu1, Kn_mu1, gamma_mu1, weights1, &
					Amat1,f11Coeffs,f12Coeffs,f13Coeffs(:,1),phi1b,phi1f)
ELSE 
	CALL GetCosExpCoeffs_specBC(Qs(1), skD, Ls(1), vgx1, vgy1, vgz1, tau1, &
					freq1, C_mu1, Kn_mu1, gamma_mu1, weights1, &
					Amat1,f11Coeffs,f12Coeffs,f13Coeffs(:,1),phi1b,phi1f)
END

IF (bc(2).EQ.'b') THEN
	CALL GetCosExpCoeffs(Qs(2), skD, Ls(2), vgx2, vgy2, vgz2, tau2, &
					freq2, C_mu2, Kn_mu2, gamma_mu2, weights2, &
					Amat2,f21Coeffs,f22Coeffs,f23Coeffs(:,1),phi2b,phi2f)
ELSE
	CALL GetCosExpCoeffs_specBC(Qs(2), skD, Ls(2), vgx2, vgy2, vgz2, tau2, &
					freq2, C_mu2, Kn_mu2, gamma_mu2, weights2, &
					Amat2,f21Coeffs,f22Coeffs,f23Coeffs(:,1),phi2b,phi2f)
END

phi1bU = phi1b
phi2bU = phi2b
phi1fU = phi1f
phi2fU = phi2f
				
CALL sDIAG_MATMUL_R_SUB(phi1b,weights1);		! phi1b = (N+1 x w1)
CALL sDIAG_MATMUL_R_SUB(phi1f,weights1);
CALL sDIAG_MATMUL_R_SUB(phi2b,weights2);
CALL sDIAG_MATMUL_R_SUB(phi2f,weights2);

PRINT *, 'successful cosExp call'




ALLOCATE(term3_1f(w1,1),term3_2b(w2,1))

PRINT *, 'Relevant num elements: w1',w1,'w2',w2,'w1d',w1d,'w2d',w2d

DEALLOCATE(kx1,kx2,ky1,ky2,kz1,kz2,vgy1,vgy2,vgz1,vgz2)


! ################################################################
! build specular block matrices


PRINT *, 'skD1', skD1, EXP(-1/skD1)
term3_1f(:,1) = weights1*skD1/4/pi*Qs(1)* &
		(vgx1*tau1/Kn_mu1*(EXP(-1/skD1)-EXP(-gamma_mu1))/(skD1*gamma_mu1-1));
	
term3_2b(:,1) = weights2*skD2/4/pi*Qs(2)* &
		(vgx2*tau2/Kn_mu2*(1-EXP(-1/skD2-gamma_mu2))/(skD2*gamma_mu2+1));

PRINT *, 'term3s defined'
		
x11Soln = cAXB(Amat1,f11Coeffs);
x12Soln = cAXB(Amat1,f12Coeffs);
x13Soln = cAXB(Amat1,f13Coeffs);

x21Soln = cAXB(Amat2,f21Coeffs);
x22Soln = cAXB(Amat2,f22Coeffs);
x23Soln = cAXB(Amat2,f23Coeffs);


! ########################################################################

ALLOCATE(blockB(w1,w1),blockC(w1,w1d), blockD(w1,w1),blockE(w1,w1d),temp3(w1))

blockB = sSPARSE_MATMUL_L(G12_d,cDIAG(1.0-spec1(:,1)),w1,w1);
blockC = cBLOCKSUM(blockB,(/w1,w1d/),numEachFreq1)

blockD = sSPARSE_MATMUL_R(cDIAG(cVECAVG(1.0-spec1(:,1),numEachFreq1)),G12_d,w1,w1);
		! shouldn't need cVECAVG
blockE = cBLOCKSUM(blockD,(/w1,w1d/),numEachFreq1)

temp3 = cVECAVG(1.0-spec1(:,1),numEachFreq1)

! PRINT *, blockC(1,1), blockC(2,2), blockC(3:5,3)
! PRINT *, 1.0-spec1(1:10,1)
! PRINT *, temp3(1:10)
! PRINT *, numEachFreq1(1:5)

IF (ALL(ABS(blockE-blockC) < 1.0E-20)) THEN
	PRINT *, 'UR equiv', SUM(ABS(blockE-blockC))/w1d/w1
ELSE
	PRINT *, 'UR not equiv', SUM(ABS(blockE-blockC))/w1d/w1
END IF

DEALLOCATE(blockB,blockC,blockD,blockE,temp3)



ALLOCATE(blockB(w1,w1),blockC(w1d,w1d), blockD(w1,w1),blockE(w1d,w1d))


blockB = cDIAG(-1.0*spec1(:,1))
blockC = cBLOCKSUM(cBLOCKSUM(blockB,(/w1,w1d/),numEachFreq1),(/w1d,w1d/),numEachFreq1)

blockD = sSPARSE_MATMUL_R(cDIAG(cVECAVG(-1.0*spec1(:,1),numEachFreq1)),G12_d,w1,w1)
		! shouldn't need cVECAVG
blockE = cBLOCKSUM(cBLOCKSUM(blockD,(/w1,w1d/),numEachFreq1),(/w1d,w1d/),numEachFreq1)

IF (ALL(ABS(blockE-blockC) < 1.0E-20)) THEN
	PRINT *, 'LR equiv', SUM(ABS(blockE-blockC))/w1d**2
ELSE
	PRINT *, 'LR not equiv', SUM(ABS(blockE-blockC))/w1d**2
END IF

PRINT *, blockC(1,1), blockC(2,2), blockC(3:5,3)
PRINT *, blockE(1,1), blockE(2,2), blockE(3:5,3)
PRINT *, spec1(1:10,1)
PRINT *, numEachFreq1(1:5)

DEALLOCATE(blockB,blockC,blockD,blockE)


ALLOCATE(blockB(w2,w2),blockC(w2,w2d), blockD(w2,w2),blockE(w2,w2d))

blockB = sSPARSE_MATMUL_L(G21_d,cDIAG(1.0-spec2(:,1)),w2,w2);
blockC = cBLOCKSUM(blockB,(/w2,w2d/),numEachFreq2)

blockD = sSPARSE_MATMUL_R(cDIAG(cVECAVG(1.0-spec2(:,1),numEachFreq2)),G21_d,w2,w2);
		! shouldn't need cVECAVG
blockE = cBLOCKSUM(blockD,(/w2,w2d/),numEachFreq2)

IF (ALL(ABS(blockE-blockC) < 1.0E-20)) THEN
	PRINT *, 'UR equiv', SUM(ABS(blockE-blockC))/w2d/w2
ELSE
	PRINT *, 'UR not equiv', SUM(ABS(blockE-blockC))/w2d/w2
END IF

DEALLOCATE(blockB,blockC,blockD,blockE)



ALLOCATE(blockB(w2,w2),blockC(w2d,w2d), blockD(w2,w2),blockE(w2d,w2d),temp1(w2d,w2d))


blockB = cDIAG(-1.0*spec2(:,1))
blockC = cBLOCKSUM(cBLOCKSUM(blockB,(/w2,w2d/),numEachFreq2),(/w2d,w2d/),numEachFreq2)

blockD = sSPARSE_MATMUL_R(cDIAG(cVECAVG(-1.0*spec2(:,1),numEachFreq2)),G21_d,w2,w2)
		! shouldn't need cVECAVG
blockE = cBLOCKSUM(cBLOCKSUM(blockD,(/w2,w2d/),numEachFreq2),(/w2d,w2d/),numEachFreq2)

temp1 = ABS(blockE-blockC)

IF (ALL(ABS(blockE-blockC) < 1.0E-20)) THEN
	PRINT *, 'LR equiv', SUM(ABS(blockE-blockC))/w2d**2
ELSE
	PRINT *, 'LR not equiv', SUM(ABS(blockE-blockC))/w2d**2
	loc1 = MAXLOC(MAXVAL(temp1,1),1)
	loc2 = MAXLOC(temp1(loc1,:),1)
	PRINT *, 'max err LR', temp1(loc1,loc2), loc1, loc2
END IF

DEALLOCATE(blockB,blockC,blockD,blockE,temp1)


! ########################################################################


MatABCD = 0; MatEF = 0; MatGH = 0; VecJK = 0;
CALL OMP_INIT_LOCK(LCK)

!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(blockA,blockB,blockC,blockD,blockE,blockF,blockG,blockH,blockJ,blockK)

!$OMP SINGLE
WRITE (*,*) 'Parallel part num threads: ', OMP_GET_NUM_THREADS()
!$OMP END SINGLE

!$OMP SECTIONS
!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL_A', OMP_GET_THREAD_NUM()
PRINT *, 'sec team size level 0', OMP_GET_TEAM_SIZE(0)
PRINT *, 'sec team size level 1', OMP_GET_TEAM_SIZE(1)
PRINT *, 'sec team size level 2', OMP_GET_TEAM_SIZE(2)
PRINT *, 'current level', OMP_GET_ACTIVE_LEVEL()
CALL OMP_UNSET_LOCK(LCK)

ALLOCATE(blockA(w2,w1))
blockA = -scSPARSE_MATMUL_L(T12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w2,w1);
! blockA = -scSPARSE_MATMUL_L(T12_s,MATMUL(TRANSPOSE(phi1fU),x12Soln),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockA);
MatABCD(1:w2,1:w1) = blockA;
DEALLOCATE(blockA)	

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL B', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockB(w2,w2))
! blockB = sSPARSE_MATMUL_L(G21_s,cDIAG(weights2*vgx2),w2,w2) - scSPARSE_MATMUL_L(R21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);  ! singular
blockB = cDIAG(weights2*vgx2) - scSPARSE_MATMUL_L(R21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);
! blockB = sSPARSE_MATMUL_L(G21_s,cDIAG(vgx2),w2,w2) - scSPARSE_MATMUL_L(R21_s,MATMUL(TRANSPOSE(phi2bU),x21Soln),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockB);
MatABCD(1:w2,(w1+1):wt) = blockB;
! PRINT *, 'Mat UL_B done'
DEALLOCATE(blockB)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL C', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockC(w1,w1))
! blockC = sSPARSE_MATMUL_L(G12_s,cDIAG(weights1*vgx1),w1,w1) - scSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);  ! singular
blockC = cDIAG(weights1*vgx1) - scSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);
! blockC = sSPARSE_MATMUL_L(G12_s,cDIAG(vgx1),w1,w1) - scSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1fU),x12Soln),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockC);
MatABCD((w2+1):wt,1:w1) = blockC;
! PRINT *, 'Mat UL_C done'
DEALLOCATE(blockC)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL D', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockD(w1,w2))
blockD = -1.0*scSPARSE_MATMUL_L(T21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w1,w2);
! blockD = -1.0*scSPARSE_MATMUL_L(T21_s,MATMUL(TRANSPOSE(phi2bU),x21Soln),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockD)
MatABCD((w2+1):wt,(w1+1):wt) = blockD;
! PRINT *, 'Mat UL_D done'	
DEALLOCATE(blockD)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL E', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockE(w2,w2))
blockE = scSPARSE_MATMUL_L(R21_s,(cDIAG(weights2*vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln)),w2,w2);
! blockE = scSPARSE_MATMUL_L(R21_s,(cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2bU),x22Soln)),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockE);
MatEF(1:w2,1:w2) = blockE;
! PRINT *, 'Mat UL_E done'
DEALLOCATE(blockE)

!$OMP SECTION
ALLOCATE(blockF(w1,w2))
blockF = scSPARSE_MATMUL_L(T21_s,(cDIAG(weights2*vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln)),w1,w2);
! blockF = scSPARSE_MATMUL_L(T21_s,(cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2bU),x22Soln)),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockF)
MatEF((w2+1):wt,1:w2) = blockF;
! PRINT *, 'Mat UL_F done'
DEALLOCATE(blockF)

!$OMP SECTION
ALLOCATE(blockG(w2,w1))
blockG = scSPARSE_MATMUL_L(T12_s,(MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(weights1*vgx1*EXP(-gamma_mu1))),w2,w1);
! blockG = scSPARSE_MATMUL_L(T12_s,(MATMUL(TRANSPOSE(phi1fU),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1))),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockG)
MatGH(1:w2,1:w1) = blockG;
! PRINT *, 'Mat UL_G done'
DEALLOCATE(blockG)

!$OMP SECTION
ALLOCATE(blockH(w1,w1))
blockH = scSPARSE_MATMUL_L(R12_s,(MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(weights1*vgx1*EXP(-gamma_mu1))),w1,w1);
! blockH = scSPARSE_MATMUL_L(R12_s,(MATMUL(TRANSPOSE(phi1fU),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1))),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockH)
MatGH((w2+1):wt,1:w1) = blockH;
! PRINT *, 'Mat UL_H done'
DEALLOCATE(blockH)



! ! ################################################################
! ! build LOWER LEFT (diffuse/specular) block matrices

!$OMP SECTION
! ALLOCATE(blockA(w2,w1))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LL blockA', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)

MatABCD((wt+1):(wt+w2d),1:w1) = 0.0
! DEALLOCATE(blockA)

!$OMP SECTION
ALLOCATE(blockB(w2,w2))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LL blockB', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)

blockB = sSPARSE_MATMUL_L(G21_d,cDIAG((1.0-spec2(:,1))*weights2),w2,w2);
! blockB = sSPARSE_MATMUL_R(cDIAG(cVECAVG(1.0-spec2(:,1),numEachFreq2)),G21_d,w2,w2);
		! shouldn't need cVECAVG

MatABCD((wt+1):(wt+w2d),(w1+1):wt) = cBLOCKSUM(blockB,(/w2d,w2/),numEachFreq2);
! PRINT *, 'Mat LL_B done'
DEALLOCATE(blockB)

!$OMP SECTION
ALLOCATE(blockC(w1,w1))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LL blockC', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)

blockC = sSPARSE_MATMUL_L(G12_d,cDIAG((1.0-spec1(:,1))*weights1),w1,w1);
! blockC = sSPARSE_MATMUL_R(cDIAG(cVECAVG(1.0-spec1(:,1),numEachFreq1)),G12_d,w1,w1);

MatABCD((wt+w2d+1):(wt+w2d+w1d),1:w1) = cBLOCKSUM(blockC,(/w1d,w1/),numEachFreq1);
! PRINT *, 'Mat LL_C done'
DEALLOCATE(blockC)

!$OMP SECTION
! ALLOCATE(blockD(w1,w2))
MatABCD((wt+w2d+1):(wt+w2d+w1d),(w1+1):wt) = 0.0;
! PRINT *, 'Mat LL_D done'
! DEALLOCATE(blockD)

!$OMP SECTION
! ALLOCATE( blockE(w2,w2))
MatEF((wt+1):(wt+w2d),1:w2) = 0.0;
! PRINT *, 'Mat LL_E done'
! DEALLOCATE(blockE)

!$OMP SECTION
! ALLOCATE( blockF(w1,w2))
MatEF((wt+w2d+1):(wt+w2d+w1d),1:w2) = 0.0;
! PRINT *, 'Mat LL_F done'
! DEALLOCATE(blockF)

!$OMP SECTION
! ALLOCATE( blockG(w2,w1))
MatGH((wt+1):(wt+w2d),1:w1) = 0.0;
! PRINT *, 'Mat LL_G done'
! DEALLOCATE(blockG)

!$OMP SECTION
! ALLOCATE(blockH(w1,w1))
MatGH((wt+w2d+1):(wt+w1d+w2d),1:w1) = 0.0;
! PRINT *, 'Mat LL_H done'
! DEALLOCATE(blockH)



!!!!!!!!!!!! UR section !!!!!!!!!!!!!!!

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, ' UR blockA', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockA(w2,w1))
blockA = -scSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w2,w1);
! blockA = -scSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1fU),x12Soln),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockA)
MatABCD(1:w2,(wt+1):(wt+w1d)) = cBLOCKSUM(blockA,(/w2,w1d/),numEachFreq1)
! PRINT *, 'Mat UR_A done'
DEALLOCATE(blockA)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'UR blockB', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockB(w2,w2))
blockB = sSPARSE_MATMUL_L(G21_d,cDIAG(weights2*vgx2),w2,w2) - scSPARSE_MATMUL_L(R21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);
! blockB = sSPARSE_MATMUL_L(G21_d,cDIAG(vgx2),w2,w2) - scSPARSE_MATMUL_L(R21_d,MATMUL(TRANSPOSE(phi2bU),x21Soln),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockB)
MatABCD(1:w2,(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(blockB,(/w2,w2d/),numEachFreq2);
! PRINT *, 'Mat UR_B done'
DEALLOCATE(blockB)

!$OMP SECTION
ALLOCATE(blockC(w1,w1))
blockC = sSPARSE_MATMUL_L(G12_d,cDIAG(weights1*vgx1),w1,w1) - scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);
! blockC = sSPARSE_MATMUL_L(G12_d,cDIAG(vgx1),w1,w1) - scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1fU),x12Soln),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockC)
MatABCD((w2+1):wt,(wt+1):(wt+w1d)) = cBLOCKSUM(blockC,(/w1,w1d/),numEachFreq1);
! PRINT *, 'Mat UR_C done'
DEALLOCATE(blockC)

!$OMP SECTION
ALLOCATE(blockD(w1,w2))
blockD = -scSPARSE_MATMUL_L(T21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w1,w2);
! blockD = -scSPARSE_MATMUL_L(T21_d,MATMUL(TRANSPOSE(phi2bU),x21Soln),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockD)
MatABCD((w2+1):wt,(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(blockD,(/w1,w2d/),numEachFreq2);
! PRINT *, 'Mat UR_D done'
DEALLOCATE(blockD)

!$OMP SECTION
ALLOCATE(blockE(w2,w2))			
blockE = scSPARSE_MATMUL_L(R21_d,sSPARSE_MATMUL_L(G21_d,cDIAG(weights2*vgx2*EXP(-gamma_mu2)),w2,w2) &
				+ MATMUL(TRANSPOSE(phi2b),x22Soln),w2,w2);
! blockE = scSPARSE_MATMUL_L(R21_d,cDIAG(weights2*vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w2,w2);
! blockE = scSPARSE_MATMUL_L(R21_d,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2bU),x22Soln),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockE)
MatEF(1:w2,(w2+1):(w2+w2d)) = cBLOCKSUM(blockE,(/w2,w2d/),numEachFreq2);
! PRINT *, 'Mat UR_E done'
DEALLOCATE(blockE)

!$OMP SECTION
ALLOCATE(blockF(w1,w2))
blockF = scSPARSE_MATMUL_L(T21_d,sSPARSE_MATMUL_L(G21_d,cDIAG(weights2*vgx2*EXP(-gamma_mu2)),w2,w2) &
				+ MATMUL(TRANSPOSE(phi2b),x22Soln),w1,w2);
! blockF = scSPARSE_MATMUL_L(T21_d,cDIAG(weights2*vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w1,w2);
! blockF = scSPARSE_MATMUL_L(T21_d,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2bU),x22Soln),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockF)
MatEF((w2+1):(w2+w1),(w2+1):(w2+w2d)) = cBLOCKSUM(blockF,(/w1,w2d/),numEachFreq2);
! PRINT *, 'MatUR_F done'
DEALLOCATE(blockF)

!$OMP SECTION
ALLOCATE(blockG(w2,w1))
blockG = scSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + &
				scSPARSE_MATMUL_L(G12_d,cDIAG(weights1*vgx1*EXP(-gamma_mu1)),w1,w1) ,w2,w1);
! blockG = scSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(weights1*vgx1*EXP(-gamma_mu1)) ,w2,w1);
! blockG = scSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1fU),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockG)
MatGH(1:w2,(w1+1):(w1+w1d)) = cBLOCKSUM(blockG,(/w2,w1d/),numEachFreq1);
! PRINT *, 'Mat UR_G done'
DEALLOCATE(blockG)

!$OMP SECTION
ALLOCATE(blockH(w1,w1))
blockH = scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + &
				scSPARSE_MATMUL_L(G12_d,cDIAG(weights1*vgx1*EXP(-gamma_mu1)),w1,w1) ,w1,w1);
! blockH = scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(weights1*vgx1*EXP(-gamma_mu1)),w1,w1);
! blockH = scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1fU),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockH);
MatGH((w2+1):wt,(w1+1):(w1+w1d)) = cBLOCKSUM(blockH,(/w1,w1d/),numEachFreq1);
! PRINT *, 'Mat UR_H done'
DEALLOCATE(blockH)




!!!!!!!! LR SECTION !!!!!!!!!!

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LR blockA', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
! ALLOCATE(blockA(w2d,w1d))
MatABCD((wt+1):(wt+w2d),(wt+1):(wt+w1d)) = 0.0;
! PRINT *, 'Mat LR_A done'
! DEALLOCATE(blockA)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LR blockB', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockB(w2,w2))

blockB = scSPARSE_MATMUL_L(G21_d,cDIAG(-1.0*spec2(:,1)*weights2),w2,w2)
! blockB = sSPARSE_MATMUL_R(cDIAG(cVECAVG(-1.0*spec2(:,1),numEachFreq2)),G21_d,w2,w2)
MatABCD((wt+1):(wt+w2d),(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(cBLOCKSUM(blockB,(/w2,w2d/),numEachFreq2),&
																(/w2d,w2d/),numEachFreq2);
! PRINT *, 'Mat LR_B done'
DEALLOCATE(blockB)

!$OMP SECTION
ALLOCATE(blockC(w1,w1))
blockC = scSPARSE_MATMUL_L(G12_d,cDIAG(-1.0*spec1(:,1)*weights1),w1,w1)
! blockC = sSPARSE_MATMUL_R(cDIAG(cVECAVG(-1.0*spec1(:,1),numEachFreq1)),G12_d,w1,w1)
MatABCD((wt+w2d+1):(wt+w2d+w1d),(wt+1):(wt+w1d)) = cBLOCKSUM(cBLOCKSUM(blockC,(/w1,w1d/),numEachFreq1),&
																(/w1d,w1d/),numEachFreq1);
! PRINT *, 'Mat LR_C done'
DEALLOCATE(blockC)
	
!$OMP SECTION
! ALLOCATE(blockD(w1d,w2d))
MatABCD((wt+w2d+1):(wt+w2d+w1d),(wt+w1d+1):(wt+w1d+w2d)) = 0.0;
! PRINT *, 'Mat LR_D done'
! DEALLOCATE(blockD)

!$OMP SECTION
! ALLOCATE(blockE(w2,w2))
MatEF((wt+1):(wt+w2d),(w2+1):(w2+w2d)) = 0.0;
! PRINT *, 'Mat LR_E done'
! DEALLOCATE(blockE)

!$OMP SECTION
! ALLOCATE(blockF(w1,w2))
MatEF((wt+w2d+1):(wt+w2d+w1d),(w2+1):(w2+w2d)) = 0.0;
! PRINT *, 'Mat LR_F done'
! DEALLOCATE(blockF)

!$OMP SECTION
! ALLOCATE(blockG(w2,w1))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat LR_G', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)

MatGH((wt+1):(wt+w2d),(w1+1):(w1+w1d)) = 0.0;
! PRINT *, 'Mat LR_G done'
! DEALLOCATE(blockG)

!$OMP SECTION
! ALLOCATE(blockH(w1,w1))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat LR_H', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)

MatGH((wt+w2d+1):(wt+w1d+w2d),(w1+1):(w1+w1d)) = 0.0;
! DEALLOCATE(blockH)



!!!!!!!!!!!!!!!!!!!! VecJK !!!!!!!!!!!!!!!!!!!!!!!!

!$OMP SECTION
!!!!!!!!!!!!! UL !!!!!!!!!!
ALLOCATE(blockJ(w2,1))
blockJ = scSPARSE_MATMUL_L(T12_s,spec2*(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w2,1) &
			+ scSPARSE_MATMUL_L(R21_s,spec2*(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b),w2,1);
VecJK(1:w2,:) = blockJ;
DEALLOCATE(blockJ)

ALLOCATE(blockK(w1,1))
blockK = scSPARSE_MATMUL_L(R12_s,spec1*(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w1,1) &
			+ scSPARSE_MATMUL_L(T21_s,spec1*(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b),w1,1);
VecJK((1+w2):wt,:) = blockK;
DEALLOCATE(blockK)


!!!!!!!!!!!!! UR !!!!!!!!!!!!!
ALLOCATE(blockJ(w2,1))
blockJ = scSPARSE_MATMUL_L(T12_d,(1.0-spec2)*(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w2,1) &
			+ scSPARSE_MATMUL_L(R21_d,(1.0-spec2)*(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b),w2,1);
VecJK(1:w2,:) = VecJK(1:w2,:) + blockJ;
DEALLOCATE(blockJ)

ALLOCATE(blockK(w1,1))
blockK = scSPARSE_MATMUL_L(R12_d,(1.0-spec1)*(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w1,1) &
			+ scSPARSE_MATMUL_L(T21_d,(1.0-spec1)*(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b),w1,1);
VecJK((1+w2):wt,:) = VecJK((1+w2):wt,:) + blockK;
DEALLOCATE(blockK)

!!!!!!!!!! LL, LR !!!!!!!!!!!
VecJK((wt+1):(wt+w2d),:) = 0.0;
VecJK((wt+w2d+1):(wt+w2d+w1d),:) = 0.0;


!$OMP END SECTIONS
!$OMP END PARALLEL

CALL OMP_DESTROY_LOCK(LCK)




PRINT *, 'full blockMats done'

DEALLOCATE(vgx1, freq1, tau1, C_mu1, Kn_mu1, gamma_mu1, weights1, phi1f, phi1b, phi1fU, phi1bU, &
		 x11Soln, x12Soln, x13Soln)

DEALLOCATE(vgx2, freq2, tau2, C_mu2, Kn_mu2, gamma_mu2, weights2, phi2f, phi2b, phi2fU, phi2bU, &
		 x21Soln, x22Soln, x23Soln)
		 
DEALLOCATE(term3_1f,term3_2b)
! DEALLOCATE(Amat2,f21Coeffs,f22Coeffs,f23Coeffs)

END SUBROUTINE GetInterfaceBM_vSpec