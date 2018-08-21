SUBROUTINE GetInterfaceBM_PBC(matlParams,Ls,Qs,skD, spec1,spec2, dT21, &
							G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, &
							MatABCD, MatEFGH, VecJK, Amat1, f11Coeffs, f12Coeffs, f13Coeffs, Amat2, f21Coeffs, f22Coeffs, f23Coeffs)
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
REAL, INTENT(IN) :: skD, spec1(:,:), spec2(:,:), dT21
			! specularity for each mode
TYPE(sIndVal), DIMENSION(:), INTENT(IN) :: R12_d, T12_d, T21_d, R21_d, R12_s, T12_s, T21_s, R21_s, G12_d, G12_s, G21_d, G21_s
REAL*8, INTENT(OUT) :: MatABCD(:,:), MatEFGH(:,:), VecJK(:,:), &
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

REAL*8, DIMENSION(:,:), ALLOCATABLE :: blockA, blockB, blockC, blockD, blockE, &
										blockF, blockG, blockH, blockJ, blockK, &
										term3_1f, term3_2b, bcFlux2
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
CALL GetCosExpCoeffs(Qs(1), skD, Ls(1), vgx1, vgy1, vgz1, tau1, &
				freq1, C_mu1, Kn_mu1, gamma_mu1, weights1, &
				Amat1,f11Coeffs,f12Coeffs,f13Coeffs(:,1),phi1b,phi1f)
CALL GetCosExpCoeffs(Qs(2), skD, Ls(2), vgx2, vgy2, vgz2, tau2, &
				freq2, C_mu2, Kn_mu2, gamma_mu2, weights2, &
				Amat2,f21Coeffs,f22Coeffs,f23Coeffs(:,1),phi2b,phi2f)


phi1bU = phi1b
phi2bU = phi2b
phi1fU = phi1f
phi2fU = phi2f
				
CALL sDIAG_MATMUL_R_SUB(phi1b,weights1);		! phi1b = (N+1 x w1)
CALL sDIAG_MATMUL_R_SUB(phi1f,weights1);
CALL sDIAG_MATMUL_R_SUB(phi2b,weights2);
CALL sDIAG_MATMUL_R_SUB(phi2f,weights2);

PRINT *, 'successful cosExp call'

PRINT *, 'Relevant num elements: w1',w1,'w2',w2,'w1d',w1d,'w2d',w2d

DEALLOCATE(kx1,kx2,ky1,ky2,kz1,kz2,vgy1,vgy2,vgz1,vgz2)


! ################################################################
! build specular block matrices


ALLOCATE(term3_1f(w1,1),term3_2b(w2,1), bcFlux2(w2,1))

PRINT *, 'skD1', skD1, EXP(-1/skD1)
term3_1f(:,1) = weights1*skD1/4/pi*Qs(1)* &
		(vgx1*tau1/Kn_mu1*(EXP(-1/skD1)-EXP(-gamma_mu1))/(skD1*gamma_mu1-1));
	
term3_2b(:,1) = weights2*skD2/4/pi*Qs(2)* &
		(vgx2*tau2/Kn_mu2*(1-EXP(-1/skD2-gamma_mu2))/(skD2*gamma_mu2+1));
		
bcFlux2(:,1) = C_mu2*dT21*vgx2*weights2;

PRINT *, 'term3s defined, bcfluxes defined'
		
x11Soln = cAXB(Amat1,f11Coeffs);
x12Soln = cAXB(Amat1,f12Coeffs);
x13Soln = cAXB(Amat1,f13Coeffs);

x21Soln = cAXB(Amat2,f21Coeffs);
x22Soln = cAXB(Amat2,f22Coeffs);
x23Soln = cAXB(Amat2,f23Coeffs);



! ########################################################################


MatABCD = 0; MatEFGH = 0; VecJK = 0;
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


!!! include specularity into these boundary conditions??  (UL, UR)
ALLOCATE(blockA(w2,w1))
blockA = -scSPARSE_MATMUL_L(T12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w2,w1);
MatABCD(1:w2,1:w1) = blockA;
DEALLOCATE(blockA)	

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL B', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockB(w2,w2))
blockB = cDIAG(weights2*vgx2) - scSPARSE_MATMUL_L(R21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);
MatABCD(1:w2,(w1+1):wt) = blockB;
DEALLOCATE(blockB)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL C', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockC(w1,w1))
blockC = cDIAG(weights1*vgx1) - scSPARSE_MATMUL_L(R12_s,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);
MatABCD((w2+1):wt,1:w1) = blockC;
DEALLOCATE(blockC)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL D', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockD(w1,w2))
blockD = -1.0*scSPARSE_MATMUL_L(T21_s,MATMUL(TRANSPOSE(phi2b),x21Soln),w1,w2);
MatABCD((w2+1):wt,(w1+1):wt) = blockD;
DEALLOCATE(blockD)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'Mat UL E', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockE(w2,w2))
blockE = scSPARSE_MATMUL_L(R21_s,(cDIAG(weights2*vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln)),w2,w2);
MatEFGH(1:w2,1:w2) = blockE;
DEALLOCATE(blockE)

!$OMP SECTION
ALLOCATE(blockF(w1,w2))
blockF = scSPARSE_MATMUL_L(T21_s,(cDIAG(weights2*vgx2*EXP(-gamma_mu2)) + MATMUL(TRANSPOSE(phi2b),x22Soln)),w1,w2);
MatEFGH((w2+1):wt,1:w2) = blockF;
DEALLOCATE(blockF)

!$OMP SECTION
ALLOCATE(blockG(w2,w1))
blockG = scSPARSE_MATMUL_L(T12_s,(MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(weights1*vgx1*EXP(-gamma_mu1))),w2,w1);
MatEFGH(1:w2,(w2+1):wt) = blockG;
DEALLOCATE(blockG)

!$OMP SECTION
ALLOCATE(blockH(w1,w1))
blockH = scSPARSE_MATMUL_L(R12_s,(MATMUL(TRANSPOSE(phi1f),x11Soln) + cDIAG(weights1*vgx1*EXP(-gamma_mu1))),w1,w1);
MatEFGH((w2+1):wt,(w2+1):wt) = blockH;
DEALLOCATE(blockH)			




! ! ################################################################
! ! build LOWER LEFT (diffuse/specular) block matrices

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LL blockA', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
MatABCD((wt+1):(wt+w2d),1:w1) = 0.0

ALLOCATE(blockB(w2,w2))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LL blockB', OMP_GET_THREAD_NUM()
PRINT *, 'spec defined by P,B'
CALL OMP_UNSET_LOCK(LCK)

! blockB = sSPARSE_MATMUL_L(G21_d,cDIAG((1.0-spec2(:,1))*weights2),w2,w2);
blockB = sSPARSE_MATMUL_L(G21_s,cDIAG((1.0-spec2(:,1))*weights2),w2,w2);          ! 5/15/18
! blockB = sSPARSE_MATMUL_L(G21_s,cDIAG((1.0-spec2(:,1))*vgx2*weights2),w2,w2);   ! 12/26/17, spec defined by q
MatABCD((wt+1):(wt+w2d),(w1+1):wt) = cBLOCKSUM(blockB,(/w2d,w2/),numEachFreq2);
DEALLOCATE(blockB)

ALLOCATE(blockC(w1,w1))
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LL blockC', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)

! blockC = sSPARSE_MATMUL_L(G12_d,cDIAG((1.0-spec1(:,1))*weights1),w1,w1);
blockC = sSPARSE_MATMUL_L(G12_s,cDIAG((1.0-spec1(:,1))*weights1),w1,w1);         ! 5/15/18
! blockC = sSPARSE_MATMUL_L(G12_s,cDIAG((1.0-spec1(:,1))*vgx1*weights1),w1,w1);  ! 12/26/17
MatABCD((wt+w2d+1):(wt+w2d+w1d),1:w1) = cBLOCKSUM(blockC,(/w1d,w1/),numEachFreq1);
DEALLOCATE(blockC)

MatABCD((wt+w2d+1):(wt+w2d+w1d),(w1+1):wt) = 0.0;

MatEFGH((wt+1):(wt+w2d+w1d),1:w2+w1) = 0.0;


!!!!!!!!!!!! UR section !!!!!!!!!!!!!!!

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, ' UR blockA', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockA(w2,w1))
blockA = -scSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w2,w1);
MatABCD(1:w2,(wt+1):(wt+w1d)) = cBLOCKSUM(blockA,(/w2,w1d/),numEachFreq1)
DEALLOCATE(blockA)

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'UR blockB', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockB(w2,w2))
! blockB = sSPARSE_MATMUL_L(G21_d,cDIAG(weights2*vgx2),w2,w2) - scSPARSE_MATMUL_L(R21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);
blockB = sSPARSE_MATMUL_L(G21_s,cDIAG(weights2*vgx2),w2,w2) - scSPARSE_MATMUL_L(R21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w2,w2);  ! 5/15/18
MatABCD(1:w2,(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(blockB,(/w2,w2d/),numEachFreq2);
DEALLOCATE(blockB)

!$OMP SECTION
ALLOCATE(blockC(w1,w1))
! blockC = sSPARSE_MATMUL_L(G12_d,cDIAG(weights1*vgx1),w1,w1) - scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);
blockC = sSPARSE_MATMUL_L(G12_s,cDIAG(weights1*vgx1),w1,w1) - scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x12Soln),w1,w1);  ! 5/15/18
MatABCD((w2+1):wt,(wt+1):(wt+w1d)) = cBLOCKSUM(blockC,(/w1,w1d/),numEachFreq1);
DEALLOCATE(blockC)

!$OMP SECTION
ALLOCATE(blockD(w1,w2))
blockD = -scSPARSE_MATMUL_L(T21_d,MATMUL(TRANSPOSE(phi2b),x21Soln),w1,w2);
MatABCD((w2+1):wt,(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(blockD,(/w1,w2d/),numEachFreq2);
DEALLOCATE(blockD)

!$OMP SECTION
ALLOCATE(blockE(w2,w2))			
blockE = scSPARSE_MATMUL_L(R21_d,sSPARSE_MATMUL_L(G21_d,cDIAG(weights2*vgx2*EXP(-gamma_mu2)),w2,w2) &
				+ MATMUL(TRANSPOSE(phi2b),x22Soln),w2,w2);
MatEFGH(1:w2,(wt+1):(wt+w2d)) = cBLOCKSUM(blockE,(/w2,w2d/),numEachFreq2);
DEALLOCATE(blockE)

!$OMP SECTION
ALLOCATE(blockF(w1,w2))
blockF = scSPARSE_MATMUL_L(T21_d,sSPARSE_MATMUL_L(G21_d,cDIAG(weights2*vgx2*EXP(-gamma_mu2)),w2,w2) &
				+ MATMUL(TRANSPOSE(phi2b),x22Soln),w1,w2);
MatEFGH((w2+1):(w2+w1),(wt+1):(wt+w2d)) = cBLOCKSUM(blockF,(/w1,w2d/),numEachFreq2);
DEALLOCATE(blockF)

!$OMP SECTION
ALLOCATE(blockG(w2,w1))
blockG = scSPARSE_MATMUL_L(T12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + &
				scSPARSE_MATMUL_L(G12_d,cDIAG(weights1*vgx1*EXP(-gamma_mu1)),w1,w1) ,w2,w1);
MatEFGH(1:w2,(wt+w2d+1):(wt+w2d+w1d)) = cBLOCKSUM(blockG,(/w2,w1d/),numEachFreq1);
DEALLOCATE(blockG)

!$OMP SECTION
ALLOCATE(blockH(w1,w1))
blockH = scSPARSE_MATMUL_L(R12_d,MATMUL(TRANSPOSE(phi1f),x11Soln) + &
				scSPARSE_MATMUL_L(G12_d,cDIAG(weights1*vgx1*EXP(-gamma_mu1)),w1,w1) ,w1,w1);
MatEFGH((w2+1):wt,(wt+w2d+1):(wt+w2d+w1d)) = cBLOCKSUM(blockH,(/w1,w1d/),numEachFreq1);
DEALLOCATE(blockH)



!!!!!!!! LR SECTION !!!!!!!!!!

!$OMP SECTION
CALL OMP_SET_LOCK(LCK)
PRINT *, 'LR blockB', OMP_GET_THREAD_NUM()
CALL OMP_UNSET_LOCK(LCK)
ALLOCATE(blockB(w2,w2))

! blockB = scSPARSE_MATMUL_L(G21_d,cDIAG(-1.0*spec2(:,1)*weights2),w2,w2)
blockB = scSPARSE_MATMUL_L(G21_s,cDIAG(-1.0*spec2(:,1)*weights2),w2,w2)             ! 5/15/18

! 12/26/17 version:  spec defined with heat flux instead of energy density
! blockB = scSPARSE_MATMUL_L(G21_d,cDIAG(-1.0*spec2(:,1)*vgx2*weights2),w2,w2)

! blockB = sSPARSE_MATMUL_R(cDIAG(cVECAVG(-1.0*spec2(:,1),numEachFreq2)),G21_d,w2,w2)
MatABCD((wt+1):(wt+w2d),(wt+w1d+1):(wt+w1d+w2d)) = cBLOCKSUM(cBLOCKSUM(blockB,(/w2,w2d/),numEachFreq2),&
																(/w2d,w2d/),numEachFreq2);
																
! PRINT *, 'Mat LR_B done'
DEALLOCATE(blockB)

ALLOCATE(blockC(w1,w1))
! blockC = scSPARSE_MATMUL_L(G12_d,cDIAG(-1.0*spec1(:,1)*weights1),w1,w1)
blockC = scSPARSE_MATMUL_L(G12_s,cDIAG(-1.0*spec1(:,1)*weights1),w1,w1)             ! 5/15/18

! 12/26/17
! blockC = scSPARSE_MATMUL_L(G12_d,cDIAG(-1.0*spec1(:,1)*vgx1*weights1),w1,w1)

! blockC = sSPARSE_MATMUL_R(cDIAG(cVECAVG(-1.0*spec1(:,1),numEachFreq1)),G12_d,w1,w1)
MatABCD((wt+w2d+1):(wt+w2d+w1d),(wt+1):(wt+w1d)) = cBLOCKSUM(cBLOCKSUM(blockC,(/w1,w1d/),numEachFreq1),&
																(/w1d,w1d/),numEachFreq1);
! PRINT *, 'Mat LR_C done'
DEALLOCATE(blockC)


MatABCD((wt+1):(wt+w2d),(wt+1):(wt+w1d)) = 0.0;
MatABCD((wt+w2d+1):(wt+w2d+w1d),(wt+w1d+1):(wt+w1d+w2d)) = 0.0;
MatEFGH((wt+1):(wt+w2d+w1d),(wt+1):(wt+w1d+w2d)) = 0.0;




!!!!!!!!!!!!!!!!! VecJK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$OMP SECTION
!!!!!!!!!!!!!!!!! UL !!!!!!!!!!!!!!!!!!
ALLOCATE(blockJ(w2,1))
blockJ = scSPARSE_MATMUL_L(T12_s,spec1*(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w2,1) &
			+ scSPARSE_MATMUL_L(R21_s,spec2*(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b + bcFlux2),w2,1) - spec2*bcFlux2;
VecJK(1:w2,:) = blockJ;
DEALLOCATE(blockJ)

ALLOCATE(blockK(w1,1))
blockK = scSPARSE_MATMUL_L(R12_s,spec1*(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w1,1) &
			+ scSPARSE_MATMUL_L(T21_s,spec2*(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b + bcFlux2),w1,1);
VecJK((1+w2):wt,:) = blockK;
DEALLOCATE(blockK)

!!!!!!!!!!!!!!!! UR !!!!!!!!!!!!!!!!!!!
ALLOCATE(blockJ(w2,1))
blockJ = scSPARSE_MATMUL_L(T12_d,(1.0-spec1)*(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w2,1) &
			+ scSPARSE_MATMUL_L(R21_d,(1.0-spec2)*(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b + bcFlux2),w2,1) &
			- scSPARSE_MATMUL_L(G21_d,(1.0-spec2)*bcFlux2,w2,1); 
! VecJK((wt+1):(wt+w2d),:) = cBLOCKSUM(blockJ,(/w2d,1/),numEachFreq2);
VecJK(1:w2,:) = VecJK(1:w2,:) + blockJ;
DEALLOCATE(blockJ)

ALLOCATE(blockK(w1,1))
blockK = scSPARSE_MATMUL_L(R12_d,(1.0-spec1)*(MATMUL(TRANSPOSE(phi1f),x13Soln) + term3_1f),w1,1) &
			+ scSPARSE_MATMUL_L(T21_d,(1.0-spec2)*(MATMUL(TRANSPOSE(phi2b),x23Soln) + term3_2b + bcFlux2),w1,1);
! VecJK((wt+w2d+1):(wt+w2d+w1d),:) = cBLOCKSUM(blockK,(/w1d,1/),numEachFreq1);
VecJK((1+w2):wt,:) = VecJK((1+w2):wt,:) + blockK;
DEALLOCATE(blockK)

!!!!!!!!!!!!!! LR, LL !!!!!!!!!!!!!!!!!
VecJK((wt+1):(wt+w2d+w1d),:) = 0.0;


!$OMP END SECTIONS
!$OMP END PARALLEL

CALL OMP_DESTROY_LOCK(LCK)



! ! try to make matrix less singular by changing last constraint to just define the first element in B1P2
! ! but only for matrix with pbc
! IF (dT21 > 0.0) THEN
	! MatABCD(wt+w1d+w2d,:) = 0;
	! MatABCD(wt+w1d+w2d,1) = 1;
	! VecJK(wt+w1d+w2d,1) = 1e-9;
! END IF

PRINT *, 'full blockMats done'

DEALLOCATE(vgx1, freq1, tau1, C_mu1, Kn_mu1, gamma_mu1, weights1, phi1f, phi1b, phi1fU, phi1bU, &
		 x11Soln, x12Soln, x13Soln)

DEALLOCATE(vgx2, freq2, tau2, C_mu2, Kn_mu2, gamma_mu2, weights2, phi2f, phi2b, phi2fU, phi2bU, &
		 x21Soln, x22Soln, x23Soln)
		 
DEALLOCATE(term3_1f,term3_2b,bcFlux2)

END SUBROUTINE GetInterfaceBM_PBC