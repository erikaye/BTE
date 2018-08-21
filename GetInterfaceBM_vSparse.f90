SUBROUTINE GetInterfaceBM(matlParams,Ls,Qs,skD, spec1,spec2, &
							G12s, R12s, T12s, G21s, R21s, T21s, G12d, R12d, T12d, G21d, R21d, T21d, &
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
REAL, INTENT(IN) :: skD, spec1(:,:), spec2(:,:)
			! specularity for each mode
TYPE(sIndVal), DIMENSION(:), INTENT(IN) :: R12d, T12d, T21d, R21d, R12s, T12s, T21s, R21s, G12d, G12s, G21d, G21s
REAL*8, INTENT(OUT) :: MatABCD(:,:), MatEF(:,:), MatGH(:,:), VecJK(:,:), &
						Amat1(:,:), f11Coeffs(:,:), f12Coeffs(:,:), f13Coeffs(:,:), &
						Amat2(:,:), f21Coeffs(:,:), f22Coeffs(:,:), f23Coeffs(:,:)
! REAL*8, ALLOCATABLE :: Amat2(:,:), f21Coeffs(:,:), f22Coeffs(:,:), f23Coeffs(:,:)


INTEGER(KIND=OMP_LOCK_KIND) LCK

INTEGER :: w1, w2, w1d, w2d, wt, statvar, m
REAL :: skD1, skD2, rho1, rho2

REAL, DIMENSION(:), ALLOCATABLE :: kx1,kx2,ky1,ky2,kz1,kz2
REAL, DIMENSION(:), ALLOCATABLE :: vgx1, vgy1, vgz1, vgx2, vgy2, vgz2, freq1, freq2, tau1, tau2, &
									C_mu1, C_mu2, Kn_mu1, Kn_mu2, weights1, weights2
REAL*8, DIMENSION(:), ALLOCATABLE :: gamma_mu1, gamma_mu2

INTEGER, ALLOCATABLE :: numEachFreq1(:),numEachFreq2(:)
REAL, ALLOCATABLE :: uniqueFreqs1(:),uniqueFreqs2(:)

REAL*8, DIMENSION(:,:), ALLOCATABLE :: x11Soln, x12Soln, x13Soln, x21Soln, x22Soln, x23Soln, &
										phi1f, phi1b, phi2f, phi2b
TYPE(sIndVal), DIMENSION(:), ALLOCATABLE :: R12_d, T12_d, T21_d, R21_d, R12_s, T12_s, T21_s, R21_s, G12_d, G12_s, G21_d, G21_s
											! weighted (along col) version of R12d, T12d, ...
										
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
		 weights1(w1), numEachFreq1(w1d), uniqueFreqs1(w1d), phi1f(N+1,w1), phi1b(N+1,w1), &
		 x11Soln(N+1,w1), x12Soln(N+1,w1), x13Soln(N+1,1))

ALLOCATE(vgx2(w2), vgy2(w2), vgz2(w2), freq2(w2), tau2(w2), C_mu2(w2), Kn_mu2(w2), gamma_mu2(w2), &
		 weights2(w2), numEachFreq2(w2d), uniqueFreqs2(w2d), phi2f(N+1,w2), phi2b(N+1,w2), &
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
CALL GetCosExpCoeffs(Qs(1), skD, Ls(1), vgx1, vgy1, vgz1, tau1, &
				freq1, C_mu1, Kn_mu1, gamma_mu1, weights1, &
				Amat1,f11Coeffs,f12Coeffs,f13Coeffs(:,1),phi1b,phi1f)
CALL GetCosExpCoeffs(Qs(2), skD, Ls(2), vgx2, vgy2, vgz2, tau2, &
				freq2, C_mu2, Kn_mu2, gamma_mu2, weights2, &
				Amat2,f21Coeffs,f22Coeffs,f23Coeffs(:,1),phi2b,phi2f)

PRINT *, 'successful cosExp call'




ALLOCATE(term3_1f(w1,1),term3_2b(w2,1))

PRINT *, 'Relevant num elements: w1',w1,'w2',w2,'w1d',w1d,'w2d',w2d

DEALLOCATE(kx1,kx2,ky1,ky2,kz1,kz2,vgy1,vgy2,vgz1,vgz2)


! ################################################################
! build specular block matrices


PRINT *, 'skD1', skD1, EXP(-1/skD1)
term3_1f(:,1) = skD1/4/pi*Qs(1)* &
		(vgx1*tau1/Kn_mu1*(EXP(-1/skD1)-EXP(-gamma_mu1))/(skD1*gamma_mu1-1));
	
term3_2b(:,1) = skD2/4/pi*Qs(2)* &
		(vgx2*tau2/Kn_mu2*(1-EXP(-1/skD2-gamma_mu2))/(skD2*gamma_mu2+1));

PRINT *, 'term3s defined'
		
x11Soln = cAXB(Amat1,f11Coeffs);
x12Soln = cAXB(Amat1,f12Coeffs);
x13Soln = cAXB(Amat1,f13Coeffs);

x21Soln = cAXB(Amat2,f21Coeffs);
x22Soln = cAXB(Amat2,f22Coeffs);
x23Soln = cAXB(Amat2,f23Coeffs);



! ########################################################################
! modified R, T, G matrices
ALLOCATE(R12_d(SIZE(R12d)), T12_d(SIZE(T12d)), T21_d(SIZE(T21d)), R21_d(SIZE(R21d)), &
			R12_s(SIZE(R12s)), T12_s(SIZE(T12s)), T21_s(SIZE(T21s)), R21_s(SIZE(R21s)), &
			G12_d(SIZE(G12d)), G12_s(SIZE(G12s)), G21_d(SIZE(G21d)), G21_s(SIZE(G21s)))



G12_d = sSPARSE_DIAG(G12d,weights1);	G12_s = sSPARSE_DIAG(G12s,weights1)
R12_d = sSPARSE_DIAG(R12d,weights1);	R12_s = sSPARSE_DIAG(R12s,weights1)
T12_d = sSPARSE_DIAG(T12d,weights1);	T12_s = sSPARSE_DIAG(T12s,weights1)

G21_d = sSPARSE_DIAG(G21d,weights2);	G21_s = sSPARSE_DIAG(G21s,weights2)
R21_d = sSPARSE_DIAG(R21d,weights2);	R21_s = sSPARSE_DIAG(R21s,weights2)
T21_d = sSPARSE_DIAG(T21d,weights2);	T21_s = sSPARSE_DIAG(T21s,weights2)

OPEN(1,FILE = './output_R12s_new.dat')
WRITE(1,*) w1, w1, SIZE(R12_s)
DO m = 1,SIZE(R12_s)
	! CALL OMP_SET_LOCK(LCK)
	WRITE(1,*) R12_s(m)%row, R12_s(m)%col, R12_s(m)%indVal
	! CALL OMP_UNSET_LOCK(LCK)
END DO
CLOSE(1)

OPEN(2,FILE = './output_T12s_new.dat')
WRITE(2,*) w2, w1, SIZE(T12_s)
DO m = 1,SIZE(T12_s)
	WRITE(2,*) T12_s(m)%row, T12_s(m)%col, T12_s(m)%indVal
END DO
CLOSE(2)

! G12_d = G12d;	G12_s = G12s
! R12_d = R12d;	R12_s = R12s
! T12_d = T12d;	T12_s = T12s

! G21_d = G21d;	G21_s = G21s
! R21_d = R21d;	R21_s = R21s
! T21_d = T21d;	T21_s = T21s


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
! CALL sDIAG_MATMUL_L_SUB(weights2,blockA);
MatABCD(1:w2,1:w1) = blockA;
DEALLOCATE(blockA)	

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL B', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockB(w2,w2))
blockB = sSPARSE_MATMUL_L(G21_s,cDIAG(vgx2),w2,w2) - sSPARSE_MATMUL_L(R21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockB);
MatABCD(1:w2,(w1+1):(w1+w2)) = blockB;
DEALLOCATE(blockB)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL C', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockC(w1,w1))
blockC = sSPARSE_MATMUL_L(G12_s,cDIAG(vgx1),w1,w1) - sSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockC);
MatABCD((w2+1):(w2+w1),1:w1) = blockC;
DEALLOCATE(blockC)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL D', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockD(w1,w2))
blockD = -1.0*sSPARSE_MATMUL_L(T21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockD)
MatABCD((w2+1):(w2+w1),(w1+1):(w1+w2)) = blockD;
DEALLOCATE(blockD)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL E', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockE(w2,w2))
blockE = sSPARSE_MATMUL_L(R21_s,(cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln)),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockE);
MatEF(1:w2,1:w2) = blockE;
DEALLOCATE(blockE)

!$OMP SECTION
ALLOCATE(blockF(w1,w2))
blockF = sSPARSE_MATMUL_L(T21_s,(cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln)),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockF)
MatEF((1+w2):(w1+w2),1:w2) = blockF;
DEALLOCATE(blockF)

!$OMP SECTION
ALLOCATE(blockG(w2,w1))
blockG = sSPARSE_MATMUL_L(T12_s,(MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1))),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockG)
MatGH(1:w2,1:w1) = blockG;
DEALLOCATE(blockG)

!$OMP SECTION
ALLOCATE(blockH(w1,w1))
blockH = scSPARSE_MATMUL_L(R12_s,(MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1))),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockH)
MatGH((1+w2):(w1+w2),1:w1) = blockH;
DEALLOCATE(blockH)

!$OMP SECTION
ALLOCATE(blockJ(w2,1))
blockJ = sSPARSE_MATMUL_L(T12_s,(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w2,1) &
			+ sSPARSE_MATMUL_L(R21_s,(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b),w2,1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockJ);
VecJK(1:w2,:) = spec2*blockJ;
DEALLOCATE(blockJ)

ALLOCATE(blockK(w1,1))
blockK = sSPARSE_MATMUL_L(R12_s,(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w1,1) &
			+ sSPARSE_MATMUL_L(T21_s,(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b),w1,1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockK)
VecJK((1+w2):(w1+w2),:) = spec1*blockK;
DEALLOCATE(blockK)


! ! ! ################################################################
! ! ! build LOWER LEFT (diffuse/specular) block matrices

!$OMP SECTION
ALLOCATE(blockA(w2,w1))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LL blockA', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
blockA = -sSPARSE_MATMUL_L(T12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockA)
MatABCD((wt+1):(wt+w2d),1:w1) = cBLOCKSUM(blockA,(/w2d,w1/),numEachFreq2)
DEALLOCATE(blockA)

!$OMP SECTION
ALLOCATE(blockB(w2,w2))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LL blockB', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
blockB = sSPARSE_MATMUL_L(G21_s,cDIAG(vgx2),w2,w2) - sSPARSE_MATMUL_L(R21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2) + &
			sSPARSE_MATMUL_R(cDIAG((1-spec2(:,1))*cVECDIFF(vgx2,weights2,numEachFreq2)),G21d,w2,w2);
																						! use nonweighted version here
			! need cVECDIFF here bc need to integrate over those values when multiplying coeffs
			
! blockB = sSPARSE_MATMUL_L(G21_s,cDIAG(vgx2),w2,w2) - scSPARSE_MATMUL_L(R21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2) + &
			! sSPARSE_MATMUL_L(G21_d,cDIAG((1-spec2(:,1))*cVECDIFF(vgx2,weights2,numEachFreq2)),w2,w2);
! ! blockB = sSPARSE_MATMUL_L(g2pCoeff_s,sDIAG(vgx2),w2,w2) - scSPARSE_MATMUL_L(R21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2) + &
			! ! sSPARSE_MATMUL_L(g2pCoeff_d,sDIAG((1-spec2(:,1))*vgx2),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockB);
MatABCD((wt+1):(wt+w2d),(w1+1):(w1+w2)) = cBLOCKSUM(blockB,(/w2d,w2/),numEachFreq2);
DEALLOCATE(blockB)

!$OMP SECTION
ALLOCATE(blockC(w1,w1))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LL blockC', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
blockC = sSPARSE_MATMUL_L(G12_s,cDIAG(vgx1),w1,w1) - sSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1) + &
			sSPARSE_MATMUL_R(cDIAG((1-spec1(:,1))*cVECDIFF(vgx1,weights1,numEachFreq1)),G12d,w1,w1);
			! need cVECDIFF here bc need to integrate over those values when multiplying coeffs
			
! blockC = sSPARSE_MATMUL_L(G12_s,cDIAG(vgx1),w1,w1) - scSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1) + &
			! sSPARSE_MATMUL_L(G12_d,cDIAG((1-spec1(:,1))*cVECDIFF(vgx1,weights1,numEachFreq1)),w1,w1);
! ! blockC = sSPARSE_MATMUL_L(g1mCoeff_s,sDIAG(vgx1),w1,w1) - scSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1) + &
			! ! sSPARSE_MATMUL_L(g1mCoeff_d,sDIAG((1-spec1(:,1))*vgx1),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockC)
MatABCD((wt+w2d+1):(wt+w2d+w1d),1:w1) = cBLOCKSUM(blockC,(/w1d,w1/),numEachFreq1);
DEALLOCATE(blockC)

!$OMP SECTION
ALLOCATE(blockD(w1,w2))
blockD = -sSPARSE_MATMUL_L(T21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockD)
MatABCD((wt+w2d+1):(wt+w2d+w1d),(w1+1):(w1+w2)) = cBLOCKSUM(blockD,(/w1d,w2/),numEachFreq1);
DEALLOCATE(blockD)

!$OMP SECTION
ALLOCATE( blockE(w2,w2))
blockE = sSPARSE_MATMUL_L(R21_s,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w2,w2) + &
			sSPARSE_MATMUL_L(R21_d, &
					sSPARSE_MATMUL_R(cDIAG((1-spec2(:,1))*vgx2*EXP(-gamma_mu2)),G21d,w2,w2) &
				,w2,w2);
				! don't need cVECDIFF(vgx2...) here bc R21d takes care of integration over those values
				
! blockE = scSPARSE_MATMUL_L(R21_s,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w2,w2) + &
			! scSPARSE_MATMUL_L(R21_d,cDIAG((1-spec2(:,1))*cVECDIFF(vgx2*EXP(-gamma_mu2),weights2,numEachFreq2)),w2,w2);
! ! blockE = scSPARSE_MATMUL_L(R21_s,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w2,w2) + &
			! ! scSPARSE_MATMUL_L(R21_d,cDIAG((1-spec2(:,1))*vgx2*EXP(-gamma_mu2)),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockE)
MatEF((wt+1):(wt+w2d),1:w2) = cBLOCKSUM(blockE,(/w2d,w2/),numEachFreq2);
DEALLOCATE(blockE)

!$OMP SECTION
ALLOCATE( blockF(w1,w2))
blockF = sSPARSE_MATMUL_L(T21_s,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w1,w2) + &
			sSPARSE_MATMUL_L(T21_d, &
					sSPARSE_MATMUL_R(cDIAG((1-spec2(:,1))*vgx2*EXP(-gamma_mu2)),G21d,w2,w2) &
				,w1,w2);
				! don't need cVECDIFF(vgx2...) here bc R21d takes care of integration over those values
				
! blockF = scSPARSE_MATMUL_L(T21_s,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w1,w2) + &
			! scSPARSE_MATMUL_L(T21_d,cDIAG((1-spec2(:,1))*cVECDIFF(vgx2*EXP(-gamma_mu2),weights2,numEachFreq2)),w1,w2);
! ! blockF = scSPARSE_MATMUL_L(T21_s,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w1,w2) + &
			! ! scSPARSE_MATMUL_L(T21_d,cDIAG((1-spec2(:,1))*vgx2*EXP(-gamma_mu2)),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockF)
MatEF((wt+w2d+1):(wt+w2d+w1d),1:w2) = cBLOCKSUM(blockF,(/w1d,w2/),numEachFreq1);
DEALLOCATE(blockF)

!$OMP SECTION
ALLOCATE( blockG(w2,w1))
blockG = sSPARSE_MATMUL_L(T12_s,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w2,w1) + &
			sSPARSE_MATMUL_L(T12_d, &
					sSPARSE_MATMUL_R(cDIAG((1-spec1(:,1))*vgx1*EXP(-gamma_mu1)),G12d,w1,w1) & 
				,w2,w1);
				! don't need cVECDIFF(vgx2...) here bc R21d takes care of integration over those values
				
! blockG = scSPARSE_MATMUL_L(T12_s,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w2,w1) + &
			! scSPARSE_MATMUL_L(T12_d,cDIAG((1-spec1(:,1))*cVECDIFF(vgx1*EXP(-gamma_mu1),weights1,numEachFreq1)),w2,w1);
! ! blockG = scSPARSE_MATMUL_L(T12_s,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w2,w1) + &
			! ! scSPARSE_MATMUL_L(T12_d,cDIAG((1-spec1(:,1))*vgx1*EXP(-gamma_mu1)),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockG)
MatGH((wt+1):(wt+w2d),1:w1) = cBLOCKSUM(blockG,(/w2d,w1/),numEachFreq2);
DEALLOCATE(blockG)

!$OMP SECTION
ALLOCATE(blockH(w1,w1))
blockH = sSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w1,w1) + &
			sSPARSE_MATMUL_L(R12_d, &
					sSPARSE_MATMUL_R(cDIAG((1-spec1(:,1))*vgx1*EXP(-gamma_mu1)),G12d,w1,w1) &
				,w1,w1);
				! don't need cVECDIFF(vgx2...) here bc R21d takes care of integration over those values
				
! blockH = scSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w1,w1) + &
			! scSPARSE_MATMUL_L(R12_d,cDIAG((1-spec1(:,1))*cVECDIFF(vgx1*EXP(-gamma_mu1),weights1,numEachFreq1)),w1,w1);
! ! blockH = scSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w1,w1) + &
			! ! scSPARSE_MATMUL_L(R12_d,cDIAG((1-spec1(:,1))*vgx1*EXP(-gamma_mu1)),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockH)
MatGH((wt+w2d+1):(wt+w1d+w2d),1:w1) = cBLOCKSUM(blockH,(/w1d,w1/),numEachFreq1);
DEALLOCATE(blockH)

!$OMP SECTION
ALLOCATE(blockJ(w2,1))
blockJ = sSPARSE_MATMUL_L(T12_s,MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f,w2,1) &
			+ sSPARSE_MATMUL_L(R21_s,MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b,w2,1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockJ)
VecJK((wt+1):(wt+w2d),:) = cBLOCKSUM(spec2*blockJ,(/w2d,1/),numEachFreq2);
DEALLOCATE(blockJ)

ALLOCATE( blockK(w1,1))
blockK = sSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f,w1,1) &
			+ sSPARSE_MATMUL_L(T21_s,MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b,w1,1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockK)
VecJK((wt+w2d+1):(wt+w2d+w1d),:) = cBLOCKSUM(spec1*blockK,(/w1d,1/),numEachFreq1);
DEALLOCATE(blockK)


!!!!!!!!!!!! UR section !!!!!!!!!!!!!!!

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, ' UR blockA', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockA(w2,w1))
blockA = -sSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockA)
MatABCD(1:w2,(wt+1):(wt+w1d)) = cBLOCKSUM(blockA,(/w2,w1d/),numEachFreq1)
DEALLOCATE(blockA)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'UR blockB', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockB(w2,w2))
blockB = sSPARSE_MATMUL_L(G21_d,cDIAG(vgx2),w2,w2) - sSPARSE_MATMUL_L(R21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockB)
MatABCD(1:w2,(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(blockB,(/w2,w2d/),numEachFreq2);
DEALLOCATE(blockB)

!$OMP SECTION
ALLOCATE(blockC(w1,w1))
blockC = sSPARSE_MATMUL_L(G12_d,cDIAG(vgx1),w1,w1) - sSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockC)
MatABCD((w2+1):(w2+w1),(wt+1):(wt+w1d)) = cBLOCKSUM(blockC,(/w1,w1d/),numEachFreq1);
DEALLOCATE(blockC)

!$OMP SECTION
ALLOCATE(blockD(w1,w2))
blockD = -sSPARSE_MATMUL_L(T21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockD)
MatABCD((w2+1):(w2+w1),(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(blockD,(/w1,w2d/),numEachFreq2);
DEALLOCATE(blockD)

!$OMP SECTION
ALLOCATE(blockE(w2,w2))			
blockE = sSPARSE_MATMUL_L(R21_d,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockE)
MatEF(1:w2,(w2+1):(w2+w2d)) = cBLOCKSUM(blockE,(/w2,w2d/),numEachFreq2);
DEALLOCATE(blockE)

!$OMP SECTION
ALLOCATE(blockF(w1,w2))
blockF = sSPARSE_MATMUL_L(T21_d,cDIAG(vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockF)
MatEF((w2+1):(w2+w1),(w2+1):(w2+w2d)) = cBLOCKSUM(blockF,(/w1,w2d/),numEachFreq2);
DEALLOCATE(blockF)

!$OMP SECTION
ALLOCATE(blockG(w2,w1))
blockG = sSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockG)
MatGH(1:w2,(w1+1):(w1+w1d)) = cBLOCKSUM(blockG,(/w2,w1d/),numEachFreq1);
DEALLOCATE(blockG)

!$OMP SECTION
ALLOCATE(blockH(w1,w1))
blockH = sSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(vgx1*EXP(-gamma_mu1)),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockH);
MatGH((w2+1):(w1+w2),(w1+1):(w1+w1d)) = cBLOCKSUM(blockH,(/w1,w1d/),numEachFreq1);
DEALLOCATE(blockH)

!$OMP SECTION
ALLOCATE(blockJ(w2,1))
blockJ = sSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f,w2,1) &
			+ sSPARSE_MATMUL_L(R21_d,MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b,w2,1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockJ)
VecJK(1:w2,:) = VecJK(1:w2,:) + (1-spec2)*blockJ;
DEALLOCATE(blockJ)

ALLOCATE(blockK(w1,1))
blockK = sSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f,w1,1) &
			+ sSPARSE_MATMUL_L(T21_d,MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b,w1,1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockK)
VecJK((w2+1):(w2+w1),:) = VecJK((w2+1):(w2+w1),:) + (1-spec1)*blockK;
DEALLOCATE(blockK)



!!!!!!!! LR SECTION !!!!!!!!!!

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LR blockA', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockA(w2,w1))
blockA = -sSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockA)
MatABCD((wt+1):(wt+w2d),(wt+1):(wt+w1d)) = cBLOCKSUM(cBLOCKSUM(blockA,(/w2,w1d/),numEachFreq1),&
														(/w2d,w1d/),numEachFreq2)
DEALLOCATE(blockA)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LR blockB', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockB(w2,w2))
blockB = sSPARSE_MATMUL_L(G21_d,cDIAG((1-spec2(:,1))*vgx2),w2,w2) - sSPARSE_MATMUL_L(R21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);
! blockB = sSPARSE_MATMUL_L(g2pCoeff_d,sDIAG((1-spec2(:,1))*cVECDIFF(vgx2,weights2,numEachFreq2)),w2,w2) &
				! - scSPARSE_MATMUL_L(R21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockB)
MatABCD((wt+1):(wt+w2d),(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(cBLOCKSUM(blockB,(/w2,w2d/),numEachFreq2),&
																(/w2d,w2d/),numEachFreq2);
DEALLOCATE(blockB)

!$OMP SECTION
ALLOCATE(blockC(w1,w1))
blockC = sSPARSE_MATMUL_L(G12_d,cDIAG((1-spec1(:,1))*vgx1),w1,w1) - sSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);
! blockC = sSPARSE_MATMUL_L(g1mCoeff_d,sDIAG((1-spec1(:,1))*cVECDIFF(vgx1,weights1,numEachFreq1)),w1,w1) & 
				! - scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockC)
MatABCD((wt+w2d+1):(wt+w2d+w1d),(wt+1):(wt+w1d)) = cBLOCKSUM(cBLOCKSUM(blockC,(/w1,w1d/),numEachFreq1),&
																		(/w1d,w1d/),numEachFreq1);
DEALLOCATE(blockC)
	
!$OMP SECTION
ALLOCATE(blockD(w1,w2))
blockD = -sSPARSE_MATMUL_L(T21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockD)
MatABCD((wt+w2d+1):(wt+w2d+w1d),(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(cBLOCKSUM(blockD,(/w1,w2d/),numEachFreq2),&
																				(/w1d,w2d/),numEachFreq1);
DEALLOCATE(blockD)

!$OMP SECTION
ALLOCATE(blockE(w2,w2))
blockE = sSPARSE_MATMUL_L(R21_d,cDIAG((1-spec2(:,1))*vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w2,w2);
! blockE = scSPARSE_MATMUL_L(R21_d,cDIAG((1-spec2(:,1))*cVECDIFF(vgx2*EXP(-gamma_mu2),weights2,numEachFreq2)) &
									! + MATMUL(TRANSPOSE(phi2b),x22Soln),w2,w2);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockE)
MatEF((wt+1):(wt+w2d),(w2+1):(w2+w2d)) = cBLOCKSUM(cBLOCKSUM(blockE,(/w2,w2d/),numEachFreq2),&
																(/w2d,w2d/),numEachFreq2);
DEALLOCATE(blockE)

!$OMP SECTION
ALLOCATE(blockF(w1,w2))
blockF = sSPARSE_MATMUL_L(T21_d,cDIAG((1-spec2(:,1))*vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln),w1,w2);
! blockF = scSPARSE_MATMUL_L(T21_d,cDIAG((1-spec2(:,1))*cVECDIFF(vgx2*EXP(-gamma_mu2),weights2,numEachFreq2)) &
									! + MATMUL(TRANSPOSE(phi2b),x22Soln),w1,w2);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockF)
MatEF((wt+w2d+1):(wt+w2d+w1d),(w2+1):(w2+w2d)) = cBLOCKSUM(cBLOCKSUM(blockF,(/w1,w2d/),numEachFreq2),&
																(/w1d,w2d/),numEachFreq1);
DEALLOCATE(blockF)

!$OMP SECTION
ALLOCATE(blockG(w2,w1))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat LR_G', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
blockG = sSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG((1-spec1(:,1))*vgx1*EXP(-gamma_mu1)),w2,w1);
! blockG = scSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) &
								! + cDIAG((1-spec1(:,1))*cVECDIFF(vgx1*EXP(-gamma_mu1),weights1,numEachFreq1)),w2,w1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockG)
MatGH((wt+1):(wt+w2d),(w1+1):(w1+w1d)) = cBLOCKSUM(cBLOCKSUM(blockG,(/w2,w1d/),numEachFreq1),&
																(/w2d,w1d/),numEachFreq2);
PRINT *, 'Mat LR_G done'
DEALLOCATE(blockG)

!$OMP SECTION
ALLOCATE(blockH(w1,w1))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat LR_H', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
blockH = sSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG((1-spec1(:,1))*vgx1*EXP(-gamma_mu1)),w1,w1);
! blockH = scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) &
								! + cDIAG((1-spec1(:,1))*cVECDIFF(vgx1*EXP(-gamma_mu1),weights1,numEachFreq1)),w1,w1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockH)
MatGH((wt+w2d+1):(wt+w1d+w2d),(w1+1):(w1+w1d)) = cBLOCKSUM(cBLOCKSUM(blockH,(/w1,w1d/),numEachFreq1),&
																(/w1d,w1d/),numEachFreq1);
DEALLOCATE(blockH)

!$OMP SECTION
ALLOCATE(blockJ(w2,1))
! PRINT *, 'Mat LR_J', OMP_GET_THREAD_NUM()
blockJ = sSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f,w2,1) &
			+ sSPARSE_MATMUL_L(R21_d,MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b,w2,1);
! CALL sDIAG_MATMUL_L_SUB(weights2,blockJ)
VecJK((wt+1):(wt+w2d),:) = VecJK((wt+1):(wt+w2d),:) + cBLOCKSUM((1-spec2)*blockJ,(/w2d,1/),numEachFreq2);
DEALLOCATE(blockJ)

ALLOCATE(blockK(w1,1))
! PRINT *, 'Mat LR_K', OMP_GET_THREAD_NUM()
blockK = sSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f,w1,1) &
			+ sSPARSE_MATMUL_L(T21_d,MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b,w1,1);
! CALL sDIAG_MATMUL_L_SUB(weights1,blockK)
VecJK((wt+w2d+1):(wt+w2d+w1d),:) = VecJK((wt+w2d+1):(wt+w2d+w1d),:) + cBLOCKSUM((1-spec1)*blockK,(/w1d,1/),numEachFreq1);
DEALLOCATE(blockK)

!$OMP END SECTIONS
!$OMP END PARALLEL

CALL OMP_DESTROY_LOCK(LCK)




PRINT *, 'full blockMats done'

DEALLOCATE(vgx1, freq1, tau1, C_mu1, Kn_mu1, gamma_mu1, weights1, phi1f, phi1b, &
		 x11Soln, x12Soln, x13Soln)

DEALLOCATE(vgx2, freq2, tau2, C_mu2, Kn_mu2, gamma_mu2, weights2, phi2f, phi2b, &
		 x21Soln, x22Soln, x23Soln)
		 
DEALLOCATE(R12_d, T12_d, T21_d, R21_d, R12_s, T12_s, T21_s, R21_s, G12_d, G12_s, G21_d, G21_s)
		 
DEALLOCATE(term3_1f,term3_2b)
! DEALLOCATE(Amat2,f21Coeffs,f22Coeffs,f23Coeffs)

END SUBROUTINE GetInterfaceBM