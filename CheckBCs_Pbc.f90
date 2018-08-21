SUBROUTINE CheckBCs_Pbc(fileNum, write_path, matlParams, dT21, Qs, skD, Ls, totL, numPer, &
						fX1Coeffs, fX2Coeffs, B1wvecs_s, B1wvecs_d, P2wvecs_s, P2wvecs_d, &
						G12_s, R12_s, T12_s, G21_s, R21_s, T21_s, G12_d, R12_d, T12_d, G21_d, R21_d, T21_d, spec1, spec2, &
						fBoth, numFBoth, weightsMatch, &
						dTCoeffs_s, dTCoeffs_d, kCond)
! with gaussian smearing for matching frequencies

USE Constants
USE omp_lib
IMPLICIT NONE
						 
TYPE(material), DIMENSION(2), INTENT(IN) :: matlParams
CHARACTER(1), INTENT(IN) :: fileNum
CHARACTER(longStr), INTENT(IN) :: write_path
REAL, INTENT(IN) :: Ls(:), Qs(:), skD, totL, dT21
REAL*8, DIMENSION(:,:,:), INTENT(IN) :: fX1Coeffs, fX2Coeffs
REAL*8, DIMENSION(:,:,:), INTENT(IN) :: B1wvecs_s, B1wvecs_d, P2wvecs_s, P2wvecs_d
TYPE(sIndVal), DIMENSION(:), INTENT(IN) :: R12_d, T12_d, T21_d, R21_d, R12_s, T12_s, T21_s, R21_s, G12_d, G12_s, G21_d, G21_s
INTEGER, INTENT(IN) :: fBoth(:,:), numFBoth, numPer
REAL*8, DIMENSION(:), INTENT(IN) :: weightsMatch
REAL*8, DIMENSION(:,:), INTENT(IN) :: dTCoeffs_d, dTCoeffs_s, spec1, spec2

REAL*8, DIMENSION(:,:), ALLOCATABLE :: R1d,R2d,T1d,T2d,R1s,R2s,T1s,T2s

REAL*8, DIMENSION(:,:,:,:), ALLOCATABLE :: intGs, intQs
REAL*8, DIMENSION(:), ALLOCATABLE :: term3_1f, term3_2b
REAL*8, DIMENSION(:,:), ALLOCATABLE :: temp1_s, temp2_s,temp1_d, temp2_d, qInttemp1, qInttemp2, bcFlux2
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: bcErr_s, bcErr_d, bcErr, specErr
REAL*8, DIMENSION(:), INTENT(OUT) :: kCond

! material params
INTEGER :: w1,w2,w1d,w2d, m, nL, mat1,mat2
REAL :: skD1,skD2, norm1,norm2, rho1, rho2
REAL, DIMENSION(:,:), ALLOCATABLE :: kx1,kx2,ky1,ky2,kz1,kz2, &
									vgx1, vgy1, vgz1, vgx2, vgy2, vgz2, freq1, freq2, tau1, tau2, &
									C_mu1, C_mu2, Kn_mu1, Kn_mu2, gamma_mu1, gamma_mu2
REAL, DIMENSION(:,:), ALLOCATABLE :: weights1, weights2
REAL*8, DIMENSION(:,:), ALLOCATABLE :: phi2b, phi1f, phi2b_dT, phi1f_dT

INTEGER, DIMENSION(:), ALLOCATABLE :: numEachFreq1,numEachFreq2
REAL, DIMENSION(:), ALLOCATABLE :: uniqueFreqs1,uniqueFreqs2
! REAL, DIMENSION(:,:), ALLOCATABLE :: temp, dfreq


! heat fluxes
REAL*8, DIMENSION(:,:,:), ALLOCATABLE :: g1pLwTerm_s,g1pLwTerm_d, g2p0wTerm_s,g2p0wTerm_d, &
										g1mLwTerm_s,g1mLwTerm_d, g2m0wTerm_s,g2m0wTerm_d
INTEGER :: m1, m2, indNum1, indNum2
REAL :: fsig1, fsig2, fsig
REAL*8 :: weightT





PRINT *, 'in checkBCs'
kCond = 0;
! nL = SIZE(Ls)*numPer
nL = 2;		! consider only 2 layers at once--fct run multiple times in MultilayerT_BBc_fct
PRINT *, nL

PRINT *, Ls(1), Ls(2), skD
PRINT *, skD/Ls(1)

w1 = matlParams(1)%wTot;
w2 = matlParams(2)%wTot;
w1d = matlParams(1)%numUniqueFreq;
w2d = matlParams(2)%numUniqueFreq;


ALLOCATE(g1pLwTerm_s(MAX(w1,w2),1,nL),g1pLwTerm_d(MAX(w1,w2),1,nL),g2m0wTerm_s(MAX(w1,w2),1,nL),&
			g2m0wTerm_d(MAX(w1,w2),1,nL),g2p0wTerm_s(MAX(w1,w2),1,nL),g2p0wTerm_d(MAX(w1,w2),1,nL),&
			g1mLwTerm_s(MAX(w1,w2),1,nL),g1mLwTerm_d(MAX(w1,w2),1,nL))
ALLOCATE(bcErr_s(w2+w1,1,nL-1), bcErr_d(w2+w1,1,nL-1), bcErr(w2+w1,1,nL-1), specErr(w1d+w2d,1,nL-1), &
			intGs(MAX(w1,w2),1,nL-1,8), intQs(MAX(w1d,w2d),1,nL-1,4))


DO m = 1,nL-1

	! get material parameters
	! either 1 or 2, depending on interface

	mat1 = MODULO(m,nL/numPer)
	mat2 = MODULO(m+1,nL/numPer)
	IF (mat1 == 0) mat1 = nL/numPer;
	IF (mat2 == 0) mat2 = nL/numPer;
	
	mat1 = m;
	mat2 = m + 1;
	
	skD1 = skD/Ls(mat1);
	skD2 = skD/Ls(mat2);
	
	w1 = matlParams(mat1)%wTot;
	w2 = matlParams(mat2)%wTot;
	w1d = matlParams(mat1)%numUniqueFreq;
	w2d = matlParams(mat2)%numUniqueFreq;

	PRINT *, 'w1d',w1d
	PRINT *, 'w2d',w2d

	ALLOCATE(kx1(w1,1),kx2(w2,1),ky1(w1,1),ky2(w2,1),kz1(w1,1),kz2(w2,1), freq1(w1,1), freq2(w2,1), vgx1(w1,1),vgx2(w2,1), tau1(w1,1), tau2(w2,1), &
				weights1(w1,1), weights2(w2,1), C_mu1(w1,1), C_mu2(w2,1), Kn_mu1(w1,1),Kn_mu2(w2,1), gamma_mu1(w1,1),gamma_mu2(w2,1),  &
				numEachFreq1(w1d), numEachFreq2(w2d),uniqueFreqs1(w1d),uniqueFreqs2(w2d), &
				term3_1f(w1),term3_2b(w2),bcFlux2(w2,1),phi2b(N+1,w2),phi1f(N+1,w1),phi2b_dT(w2,N+1),phi1f_dT(w1,N+1))

	PRINT *, 'allocated matl params'
				
	kx1(:,1) = matlParams(mat1)%kx(1:w1);
	kx2(:,1) = matlParams(mat2)%kx(1:w2);
	ky1(:,1) = matlParams(mat1)%ky(1:w1);
	ky2(:,1) = matlParams(mat2)%ky(1:w2);
	kz1(:,1) = matlParams(mat1)%kz(1:w1);
	kz2(:,1) = matlParams(mat2)%kz(1:w2);
	freq1(:,1) = matlParams(mat1)%freq(1:w1);
	freq2(:,1) = matlParams(mat2)%freq(1:w2);
	vgx1(:,1) = matlParams(mat1)%vgx(1:w1);
	vgx2(:,1) = matlParams(mat2)%vgx(1:w2);
	tau1(:,1) = matlParams(mat1)%tau(1:w1);
	tau2(:,1) = matlParams(mat2)%tau(1:w2);
	weights1(:,1) = matlParams(mat1)%weights(1:w1);
	weights2(:,1) = matlParams(mat2)%weights(1:w2);
	C_mu1(:,1) = matlParams(mat1)%C(1:w1);
	C_mu2(:,1) = matlParams(mat2)%C(1:w2);
	Kn_mu1(:,1) = matlParams(mat1)%MFP(1:w1) / Ls(mat1);
	Kn_mu2(:,1) = matlParams(mat2)%MFP(1:w2) / Ls(mat2);
	gamma_mu1(:,1) = matlParams(mat1)%gamma(1:w1) * Ls(mat1);
	gamma_mu2(:,1) = matlParams(mat2)%gamma(1:w2) * Ls(mat2);
	numEachFreq1 = matlParams(mat1)%numEachFreq(1:w1d);
	numEachFreq2 = matlParams(mat2)%numEachFreq(1:w2d);
	uniqueFreqs1 = matlParams(mat1)%uniqueFreqs(1:w1d);
	uniqueFreqs2 = matlParams(mat2)%uniqueFreqs(1:w2d);
	rho1 = matlParams(mat1)%rho;
	rho2 = matlParams(mat2)%rho;


	! calculate needed parameters...

	PRINT *, 'calculating term3s'
	PRINT *, 'skD1', skD1, EXP(-1/skD1)
	term3_1f = skD1/4/pi*Qs(m)* &
			(vgx1(:,1)*tau1(:,1)/Kn_mu1(:,1)*(EXP(-1.0/skD1)-EXP(-gamma_mu1(:,1)))/(skD1*gamma_mu1(:,1)-1.0));
	PRINT *, '1 done'
	term3_2b = skD2/4.0/pi*Qs(m+1)* &
			(vgx2(:,1)*tau2(:,1)/Kn_mu2(:,1)*(1.0-EXP(-1.0/skD2-gamma_mu2(:,1)))/(skD2*gamma_mu2(:,1)+1.0));

	bcFlux2 = C_mu2*dT21*vgx2*weights2;

	norm1 = 2.0*SUM(C_mu1(:,1)/tau1(:,1)*weights1(1:w1,1))
	norm2 = 2.0*SUM(C_mu2(:,1)/tau2(:,1)*weights2(1:w2,1))
	PRINT *, norm1, norm2

	PRINT *, 'calculating phis'
	PRINT *, SHAPE(fX1Coeffs), SHAPE(phi2b)
	phi2b(:,:) = sDIAG_MATMUL_R(fX1Coeffs(:,1:w2,m+1),(C_mu2(:,1)*Ls(mat2)*(norm2/weights2(:,1)/2.0)));
	phi1f(:,:) = sDIAG_MATMUL_R(fX2Coeffs(:,1:w1,m),(C_mu1(:,1)*Ls(mat1)*(norm1/weights1(:,1)/2.0)));
	phi2b(1,:) = 2.0*phi2b(1,:);
	phi1f(1,:) = 2.0*phi1f(1,:);

	phi2b_dT = sSPARSE_MATMUL_L(G21_d,TRANSPOSE(phi2b),w2,N+1)		! R21d*phi2b or T21d*phi2b but without R,T weighting
	phi1f_dT = sSPARSE_MATMUL_L(G12_d,TRANSPOSE(phi1f),w1,N+1)		! R12d*phi1f or T12d*phi1f but without R,T weighting
	
	! calculation of heat flux at interface
	PRINT *, 'heat flux calcs'
	ALLOCATE(temp1_s(w1,1),temp2_s(w2,1),temp1_d(w1,1),temp2_d(w2,1))

	temp1_s(:,1) = MATMUL(TRANSPOSE(phi1f),dTCoeffs_s(:,m)) + term3_1f;
	temp2_s(:,1) = MATMUL(TRANSPOSE(phi2b),dTCoeffs_s(:,m+1)) + term3_2b;
	temp1_d(:,1) = MATMUL(TRANSPOSE(phi1f),dTCoeffs_d(:,m));
	temp2_d(:,1) = MATMUL(TRANSPOSE(phi2b),dTCoeffs_d(:,m+1));
	! temp1_d(:,1) = MATMUL(phi1f_dT,dTCoeffs_d(:,1));
	! temp2_d(:,1) = MATMUL(phi2b_dT,dTCoeffs_d(:,2));

	PRINT *, 'temp vars defined'

	
	! column of P,B vecs correspond to layer or interest
	! third index is from 1 to nL+1 instead of 0 to nL
	
	! g1pLwTerm_s(1:w1,:,m) = scSPARSE_MATMUL_L(G12_s,(P2wvecs_s(1:w1,:,m)*vgx1*EXP(-gamma_mu1) + temp1_s)*weights1,w1,1);
	! g1pLwTerm_d(1:w1,:,m) = P2wvecs_d(1:w1,:,m)*scSPARSE_MATMUL_L(G12_d,vgx1*EXP(-gamma_mu1)*weights1,w1,1) &
								! + scSPARSE_MATMUL_L(G12_d,temp1_d*weights1,w1,1);
	g1pLwTerm_d(1:w1,:,m) = scSPARSE_MATMUL_L(G12_d,(P2wvecs_d(1:w1,:,m)*vgx1*EXP(-gamma_mu1))*weights1,w1,1) + temp1_d*weights1;
	g1pLwTerm_s(1:w1,:,m) = (P2wvecs_s(1:w1,:,m)*vgx1*EXP(-gamma_mu1) + temp1_s)*weights1;
	! g1pLwTerm_d(1:w1,1,m) = (P2wvecs_d(1:w1,1,m)*cVECDIFF(vgx1(:,1)*EXP(-gamma_mu1(:,1)),weights1(:,1),numEachFreq1) + temp1_d(:,1))*weights1(:,1);
	PRINT *, 'g1pLwTerm'

	! g2m0wTerm_s(1:w2,:,m) = scSPARSE_MATMUL_L(G21_s,(B1wvecs_s(1:w2,:,m+2)*vgx2*EXP(-gamma_mu2) + temp2_s)*weights2,w2,1);
	! g2m0wTerm_d(1:w2,:,m) = B1wvecs_s(1:w2,:,m+2)*scSPARSE_MATMUL_L(G21_d,vgx2*EXP(-gamma_mu2)*weights2,w2,1) &
								! + scSPARSE_MATMUL_L(G21_d,temp2_d*weights2,w2,1);
	g2m0wTerm_d(1:w2,:,m) = scSPARSE_MATMUL_L(G21_d,(B1wvecs_d(1:w2,:,m+2)*vgx2*EXP(-gamma_mu2))*weights2  + (1.0-spec2)*bcFlux2,w2,1) + temp2_d*weights2;
	g2m0wTerm_s(1:w2,:,m) = (B1wvecs_s(1:w2,:,m+2)*vgx2*EXP(-gamma_mu2) + temp2_s)*weights2 + spec2*bcFlux2;
	! g2m0wTerm_d(1:w2,1,m) = (B1wvecs_d(1:w2,1,m+2)*cVECDIFF(vgx2(:,1)*EXP(-gamma_mu2(:,1)),weights2(:,1),numEachFreq2) + temp2_d(:,1))*weights2(:,1);
	PRINT *, 'g2m0wTerm'


	! g2p0wTerm_s(1:w2,:,m) = scSPARSE_MATMUL_L(G21_s,P2wvecs_s(1:w2,:,m+1)*vgx2*weights2,w2,1);
	g2p0wTerm_d(1:w2,:,m) = scSPARSE_MATMUL_L(G21_d,P2wvecs_d(1:w2,:,m+1)*vgx2*weights2 + (1.0-spec2)*bcFlux2,w2,1);
	g2p0wTerm_s(1:w2,:,m) = P2wvecs_s(1:w2,:,m+1)*vgx2*weights2 + spec2*bcFlux2;
	! g2p0wTerm_d(1:w2,1,m) = P2wvecs_d(1:w2,1,m+1)*cVECDIFF(vgx2(:,1),weights2(:,1),numEachFreq2)*weights2(:,1);
	PRINT *, 'g2p0wTerm'

	! g1mLwTerm_s(1:w1,:,m) = scSPARSE_MATMUL_L(G12_s,B1wvecs_s(1:w1,:,m+1)*vgx1*weights1,w1,1);
	g1mLwTerm_d(1:w1,:,m) = scSPARSE_MATMUL_L(G12_d,B1wvecs_d(1:w1,:,m+1)*vgx1*weights1,w1,1);
	g1mLwTerm_s(1:w1,:,m) = B1wvecs_s(1:w1,:,m+1)*vgx1*weights1;
	! g1mLwTerm_d(1:w1,1,m) = B1wvecs_d(1:w1,1,m+1)*cVECDIFF(vgx1(:,1),weights1(:,1),numEachFreq1)*weights1(:,1);
	PRINT *, 'g1mLwTerm'
	
	
	
	ALLOCATE(R1s(w1,1), R2s(w2,1), T1s(w1,1), T2s(w2,1), R1d(w1,1), R2d(w2,1), T1d(w1,1), T2d(w2,1))
	R1s(:,1) = SUM_SPARSE(R12_s, (/ w1, w1 /), 2)  ! sum across cols
	R2s(:,1) = SUM_SPARSE(R21_s, (/ w2, w2 /), 2)
	T1s(:,1) = SUM_SPARSE(T21_s, (/ w1, w2 /), 2)
	T2s(:,1) = SUM_SPARSE(T12_s, (/ w2, w1 /), 2)
	R1d(:,1) = SUM_SPARSE(R12_d, (/ w1, w1 /), 2)  ! sum across cols
	R2d(:,1) = SUM_SPARSE(R21_d, (/ w2, w2 /), 2)
	T1d(:,1) = SUM_SPARSE(T21_d, (/ w1, w2 /), 2)
	T2d(:,1) = SUM_SPARSE(T12_d, (/ w2, w1 /), 2)	
	
	
	! ! ! calculation of error
	! bcErr_s(1:w2,:,m) = g2p0wTerm_s(1:w2,:,m) - T2s*g1pLwTerm_s(1:w2,:,m) - R2s*g2m0wTerm_s(1:w2,:,m);
	! bcErr_s((w2+1):(w2+w1),:,m) = g1mLwTerm_s(1:w1,:,m) - R1s*g1pLwTerm_s(1:w1,:,m) - T1s*g2m0wTerm_s(1:w1,:,m)
							
	! ! bcErr_d(1:w2,:,m) = g2p0wTerm_d(1:w2,:,m) - T2d*g1pLwTerm_d(1:w2,:,m) - R2d*g2m0wTerm_d(1:w2,:,m);
	! ! bcErr_d((w2+1):(w2+w1),:,m) = g1mLwTerm_d(1:w1,:,m) - R1d*g1pLwTerm_d(1:w1,:,m) - T1d*g2m0wTerm_d(1:w1,:,m);
							 
	! ! bcErr_s(1:w2,:,m) = scSPARSE_MATMUL_L(G21_s,g2p0wTerm_s(1:w2,:,m),w2,1) - scSPARSE_MATMUL_L(T12_s,g1pLwTerm_s(1:w1,:,m),w2,1) &
							! ! - scSPARSE_MATMUL_L(R21_s,g2m0wTerm_s(1:w2,:,m),w2,1)! + &
			  ! ! ! scSPARSE_MATMUL_L(g2pCoeff_d,g2p0wTerm_d,w2,1) - scSPARSE_MATMUL_L(T12_d,g1pLwTerm_d,w2,1) &
							! ! ! - scSPARSE_MATMUL_L(R21_d,g2m0wTerm_d,w2,1);
	! ! bcErr_s((w2+1):(w2+w1),:,m) = scSPARSE_MATMUL_L(G12_s,g1mLwTerm_s(1:w1,:,m),w1,1) - scSPARSE_MATMUL_L(R12_s,g1pLwTerm_s(1:w1,:,m),w1,1) &
							! ! - scSPARSE_MATMUL_L(T21_s,g2m0wTerm_s(1:w2,:,m),w1,1)! + &
			  ! ! ! scSPARSE_MATMUL_L(g1mCoeff_d,g1mLwTerm_d,w1,1) - scSPARSE_MATMUL_L(R12_d,g1pLwTerm_d,w1,1) &
							! ! ! - scSPARSE_MATMUL_L(T21_d,g2m0wTerm_d,w1,1);

	bcErr_s(1:w2,:,m) = g2p0wTerm_s(1:w2,:,m) - scSPARSE_MATMUL_L(T12_s,g1pLwTerm_s(1:w1,:,m),w2,1) &
							- scSPARSE_MATMUL_L(R21_s,g2m0wTerm_s(1:w2,:,m),w2,1)
	bcErr_s((w2+1):(w2+w1),:,m) = g1mLwTerm_s(1:w1,:,m) - scSPARSE_MATMUL_L(R12_s,g1pLwTerm_s(1:w1,:,m),w1,1) &
							- scSPARSE_MATMUL_L(T21_s,g2m0wTerm_s(1:w2,:,m),w1,1)
							
	bcErr_d(1:w2,:,m) = scSPARSE_MATMUL_L(G21_d,g2p0wTerm_d(1:w2,:,m),w2,1) - scSPARSE_MATMUL_L(T12_d,g1pLwTerm_d(1:w1,:,m),w2,1) &
							 - scSPARSE_MATMUL_L(R21_d,g2m0wTerm_d(1:w2,:,m),w2,1);
	bcErr_d((w2+1):(w2+w1),:,m) = scSPARSE_MATMUL_L(G12_d,g1mLwTerm_d(1:w1,:,m),w1,1) - scSPARSE_MATMUL_L(R12_d,g1pLwTerm_d(1:w1,:,m),w1,1) &
							 - scSPARSE_MATMUL_L(T21_d,g2m0wTerm_d(1:w2,:,m),w1,1);
	! bcErr_d(1:w2,:,m) = g2p0wTerm_d(1:w2,:,m) - scSPARSE_MATMUL_L(T12_d,g1pLwTerm_d(1:w1,:,m),w2,1) &
							 ! - scSPARSE_MATMUL_L(R21_d,g2m0wTerm_d(1:w2,:,m),w2,1);
	! bcErr_d((w2+1):(w2+w1),:,m) = g1mLwTerm_d(1:w1,:,m) - scSPARSE_MATMUL_L(R12_d,g1pLwTerm_d(1:w1,:,m),w1,1) &
							 ! - scSPARSE_MATMUL_L(T21_d,g2m0wTerm_d(1:w2,:,m),w1,1);
							 							
							
	! bcErr(1:w2d,:,m) = cBLOCKSUM(scSPARSE_MATMUL_L(G21_d,bcErr_s(1:w2,:,m)+bcErr_d(1:w2,:,m),w2,1) &
												! ,(/w2d, 1/), numEachFreq2);
	! bcErr((w2d+1):(w2d+w1d),:,m) = cBLOCKSUM(scSPARSE_MATMUL_L(G12_d,bcErr_s((w2+1):(w2+w1),:,m)+bcErr_d((w2+1):(w2+w1),:,m),w1,1) &
												! ,(/w1d,1/), numEachFreq1);
	! bcErr(1:w2d,:,m) = cBLOCKSUM(bcErr_s(1:w2,:,m)+bcErr_d(1:w2,:,m),(/w2d, 1/), numEachFreq2);
	! bcErr((w2d+1):(w2d+w1d),:,m) = cBLOCKSUM(bcErr_s((w2+1):(w2+w1),:,m)+bcErr_d((w2+1):(w2+w1),:,m),(/w1d,1/), numEachFreq1);

	! use cBlOCKSUM not cBLOCKSUM_avg bc weights are included along each row

	bcErr(1:(w2+w1),:,m) = bcErr_d(1:(w2+w1),:,m) + bcErr_s(1:(w2+w1),:,m);
	
	PRINT *, 'bcErrs defined'

	
	
	ALLOCATE(qInttemp1(w2d,1),qInttemp2(w1d,1))

	! ! ! actually calculate heat flux (vs freq)
	! intQs(1:w1d,:,m,1) = cBLOCKSUM( scSPARSE_MATMUL_L(G12_d,(g1pLwTerm_d(:,:,m) - g1mLwTerm_d(:,:,m)),w1,1) &
								! + scSPARSE_MATMUL_L(G12_d,(g1pLwTerm_s(:,:,m) - g1mLwTerm_s(:,:,m)),w1,1) &
							! ,(/w1d,1/),numEachFreq1) 
			! ! 8 from counting all octants of BZ
			! ! heat at interface in material 1
	! intQs(1:w2d,:,m,2) = cBLOCKSUM( scSPARSE_MATMUL_L(G21_d,(g2p0wTerm_d(:,:,m) - g2m0wTerm_d(:,:,m)),w2,1) &
								! + scSPARSE_MATMUL_L(G21_d,(g2p0wTerm_s(:,:,m) - g2m0wTerm_s(:,:,m)),w2,1) &
							! ,(/w2d,1/),numEachFreq2) !* 8
			! ! 8 from counting all octants of BZ
			! ! heat at interface in material 2
	! qInttemp1 = cBLOCKSUM( scSPARSE_MATMUL_L(T12_d,g1pLwTerm_d(:,:,m),w2,1) &
							! + scSPARSE_MATMUL_L(T12_s,g1pLwTerm_s(:,:,m),w2,1), (/w2d,1/), numEachFreq2 );
	! qInttemp2 = cBLOCKSUM( scSPARSE_MATMUL_L(T21_d,g2m0wTerm_d(:,:,m),w1,1) &
							! + scSPARSE_MATMUL_L(T21_s,g2m0wTerm_s(:,:,m),w1,1), (/w1d,1/), numEachFreq1 );
	
	! ! ! actually calculate heat flux (vs freq)
	! 12/10/17
	intQs(1:w1d,:,m,1) = cBLOCKSUM( scSPARSE_MATMUL_L(G12_d,(g1pLwTerm_d(1:w1,:,m) - g1mLwTerm_d(1:w1,:,m)),w1,1) & 
									+ scSPARSE_MATMUL_L(G12_s,(g1pLwTerm_s(1:w1,:,m) - g1mLwTerm_s(1:w1,:,m)),w1,1) &
							,(/w1d,1/),numEachFreq1) 
			! 8 from counting all octants of BZ
			! heat at interface in material 1
	intQs(1:w2d,:,m,2) = cBLOCKSUM( scSPARSE_MATMUL_L(G21_d,(g2p0wTerm_d(1:w2,:,m) - g2m0wTerm_d(1:w2,:,m)),w2,1) &
									+ scSPARSE_MATMUL_L(G21_s,(g2p0wTerm_s(1:w2,:,m) - g2m0wTerm_s(1:w2,:,m)),w2,1) &
							,(/w2d,1/),numEachFreq2)
			! 8 from counting all octants of BZ
			! heat at interface in material 2

	
	! previous version
	! intQs(1:w1d,:,m,1) = cBLOCKSUM((g1pLwTerm_d(1:w1,:,m) - g1mLwTerm_d(1:w1,:,m)) + (g1pLwTerm_s(1:w1,:,m) - g1mLwTerm_s(1:w1,:,m)) &
							! ,(/w1d,1/),numEachFreq1) 
			! ! 8 from counting all octants of BZ
			! ! heat at interface in material 1
	! intQs(1:w2d,:,m,2) = cBLOCKSUM((g2p0wTerm_d(1:w2,:,m) - g2m0wTerm_d(1:w2,:,m)) + (g2p0wTerm_s(1:w2,:,m) - g2m0wTerm_s(1:w2,:,m)) &
							! ,(/w2d,1/),numEachFreq2)
			! ! 8 from counting all octants of BZ
			! ! heat at interface in material 2
			
	! qInttemp1 = cBLOCKSUM( T2d*g1pLwTerm_d(1:w2,:,m) + T2s*g1pLwTerm_s(1:w2,:,m), (/w2d,1/), numEachFreq2 );
	! qInttemp2 = cBLOCKSUM( T1d*g2m0wTerm_d(1:w1,:,m) + T1s*g2m0wTerm_s(1:w1,:,m), (/w1d,1/), numEachFreq1 );
							
	qInttemp1 = cBLOCKSUM( scSPARSE_MATMUL_L(T12_d,g1pLwTerm_d(1:w1,:,m),w2,1) &
							+ scSPARSE_MATMUL_L(T12_s,g1pLwTerm_s(1:w1,:,m),w2,1), (/w2d,1/), numEachFreq2 );
	qInttemp2 = cBLOCKSUM( scSPARSE_MATMUL_L(T21_d,g2m0wTerm_d(1:w2,:,m),w1,1) &
							+ scSPARSE_MATMUL_L(T21_s,g2m0wTerm_s(1:w2,:,m),w1,1), (/w1d,1/), numEachFreq1 );
							
	intQs(1:w2d,1,m,3) = (qInttemp1(:,1) - f1tof2(qInttemp2(:,1),w2d)) !*8; 
	intQs(1:w1d,1,m,4) = (f2tof1(qInttemp1(:,1),w1d) - qInttemp2(:,1)) !*8;
			! heat crossing interface in terms of frequencies in material 1/material 2
			! intQs 3/4 != intQs 1/2 means transmission matrices are not properly defined.
	! all heat fluxes should be equal.

	DEALLOCATE(qInttemp1,qInttemp2)
	PRINT *, 'intQs defined'
	
	
		intQs(1:w1d,:,m,1) = cBLOCKSUM( scSPARSE_MATMUL_L(G12_d,(g1pLwTerm_d(1:w1,:,m) - g1mLwTerm_d(1:w1,:,m)),w1,1) & 
									+ scSPARSE_MATMUL_L(G12_s,(g1pLwTerm_s(1:w1,:,m) - g1mLwTerm_s(1:w1,:,m)),w1,1) &
							,(/w1d,1/),numEachFreq1) 
	
	! specularity error
	specErr(1:w2d,:,m) = cBLOCKSUM( g2p0wTerm_s(1:w2,:,m), (/w2d,1/), numEachFreq2)/ &
									cBLOCKSUM( g2p0wTerm_s(1:w2,:,m) + g2p0wTerm_d(1:w2,:,m), (/w2d,1/), numEachFreq2) &
									- cBLOCKSUM_avg( spec2,(/w2d,1/),numEachFreq2);
	specErr(w2d:w2d+w1d,:,m) = cBLOCKSUM( g1mLwTerm_s(1:w1,:,m), (/w1d,1/), numEachFreq1)/ &
									cBLOCKSUM( g1mLwTerm_s(1:w1,:,m) + g1mLwTerm_d(1:w1,:,m), (/w1d,1/), numEachFreq1) &
									- cBLOCKSUM_avg( spec1,(/w1d,1/),numEachFreq1);
									
	! specErr(w2d:w1d,:,m) = cBLOCKSUM(g1mLwTerm_s(1:w1,:,m)/(g1mLwTerm_s(1:w1,:,m) + g1mLwTerm_d(1:w1,:,m)) - spec1,(/w1d,1/),numEachFreq1);
	
	
	DEALLOCATE(temp1_s,temp2_s,temp1_d,temp2_d)
	DEALLOCATE(kx1,kx2,ky1,ky2,kz1,kz2, freq1,freq2, vgx1,vgx2, tau1,tau2, &
			weights1, weights2, C_mu1, C_mu2, Kn_mu1, Kn_mu2, gamma_mu1, gamma_mu2,  &
			numEachFreq1, numEachFreq2, uniqueFreqs1, uniqueFreqs2, &
			term3_1f, term3_2b, bcFlux2, phi2b, phi1f, phi2b_dT, phi1f_dT)
	DEALLOCATE(R1d,R2d,T1d,T2d,R1s,R2s,T1s,T2s)
	
END DO


intGs(1:w1,:,:,1) = g1pLwTerm_s;  ! includes dk weighting to make consistent with interface matrices
intGs(1:w1,:,:,2) = g1pLwTerm_d;
intGs(1:w1,:,:,3) = g1mLwTerm_s;
intGs(1:w1,:,:,4) = g1mLwTerm_d;
intGs(1:w2,:,:,5) = g2p0wTerm_s;
intGs(1:w2,:,:,6) = g2p0wTerm_d;
intGs(1:w2,:,:,7) = g2m0wTerm_s;
intGs(1:w2,:,:,8) = g2m0wTerm_d;



PRINT *, 'size bcErr', SIZE(bcErr(:,1,1))

DO m = 1,nL-1
	WRITE(2,"(40000(ES16.7E3,','))") bcErr_s(:,1,m)
	WRITE(3,"(40000(ES16.7E3,','))") bcErr_d(:,1,m)
	WRITE(4,"(40000(ES16.7E3,','))") bcErr(:,1,m)
	WRITE(5,"(40000(ES16.7E3,','))") specErr(:,1,m)
END DO

! ! calculate heat flux




! OPEN(3,FILE = TRIM(write_path) // 'output_intQs_' // fileNum // '.dat')
! DO m = 1,SIZE(intQs,3)
	! WRITE(3,"(40000(ES15.7,','))") intQs(:,1,m)
! END DO
! CLOSE(3)

! OPEN(3,FILE = TRIM(write_path) // 'output_intGs_' // fileNum // '.dat')
! DO m = 1,SIZE(intGs,3)
	! WRITE(3,"(40000(ES15.7,','))") intGs(:,1,m)
! END DO
! CLOSE(3)


kCond(:) = -1*SUM(intQs(1:w1d,1,:,1),1)/((dT2-dT1)/totL);
PRINT *, 'kCond no weight 1', kCond
kCond(:) = -1*SUM(intQs(1:w2d,1,:,2),1)/((dT2-dT1)/totL);
PRINT *, 'kCond no weight 2', kCond
kCond(:) = -1*SUM(intQs(1:w2d,1,:,3),1)/((dT2-dT1)/totL);
PRINT *, 'kCond no weight 3', kCond
kCond(:) = -1*SUM(intQs(1:w1d,1,:,4),1)/((dT2-dT1)/totL);
PRINT *, 'kCond no weight 4', kCond

kCond(:) = -1*SUM(intQs(1:w1d,1,:,1),1)/((dT2-dT1)/totL);
PRINT *, 'kCond defined by layer 1', kCond


DEALLOCATE(g1pLwTerm_s,g1pLwTerm_d,g2m0wTerm_s,g2m0wTerm_d)
DEALLOCATE(g2p0wTerm_s,g2p0wTerm_d,g1mLwTerm_s,g1mLwTerm_d)
			
DEALLOCATE(bcErr,bcErr_d,bcErr_s,intGs,intQs, specErr)


CONTAINS
	FUNCTION f2tof1(vec,w1d) RESULT(outVec)
		REAL*8, INTENT(IN) :: vec(:)
		INTEGER, INTENT(IN) :: w1d
		REAL*8, DIMENSION(w1d) :: outVec
		INTEGER :: m
		
		outVec = 0;
		
		DO m = 1,numFBoth
			outVec(fBoth(m,1)) = outVec(fBoth(m,1)) + vec(fBoth(m,2))!*weightsMatch(m);
			!fBoth1(m) = integer between 1,numUniqueFreq1=w1d
			!fBoth2(m) = integer between 1,numUniqueFreq2=w2d
		END DO	
	END FUNCTION
	
	FUNCTION f1tof2(vec,w2d) RESULT(outVec)
		REAL*8, INTENT(IN) :: vec(:)
		INTEGER, INTENT(IN) :: w2d
		REAL*8, DIMENSION(w2d) :: outVec
		INTEGER :: m
		
		outVec = 0;
		
		DO m = 1,numFBoth
			outVec(fBoth(m,2)) = outVec(fBoth(m,2)) + vec(fBoth(m,1))!*weightsMatch(m);
			!fBoth1(m) = integer between 1,numUniqueFreq1=w1d
			!fBoth2(m) = integer between 1,numUniqueFreq2=w2d
		END DO	
	END FUNCTION

END SUBROUTINE