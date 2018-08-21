SUBROUTINE GetCosExpCoeffs_specBC(Q_mu, skD, L, vgx, vgy, vgz, tau, freq, C_mu, Kn_mu, gamma_mu, weights, &
				Amat,fx1Coeffs,fx2Coeffs,fx3Coeffs,phixb,phixf)

USE Constants
IMPLICIT NONE

REAL, INTENT(IN) :: Q_mu, skD, L			! heat source magnitude from laser at surface of slab, skin depth, length of layer
REAL, DIMENSION(:), INTENT(IN) :: vgx, vgy, vgz, tau, freq, C_mu, Kn_mu, weights  ! single column vector containing info for all modes
REAL*8, DIMENSION(:), INTENT(IN) :: gamma_mu

REAL*8, DIMENSION(N+1,N+1), INTENT(OUT) :: Amat
REAL*8, DIMENSION(:,:), INTENT(OUT) :: fx1Coeffs, fx2Coeffs, phixb, phixf
REAL*8, DIMENSION(N+1), INTENT(OUT) :: fx3Coeffs

INTEGER :: totPts

INTEGER :: i,ip,im, m, p
INTEGER :: statvar

REAL :: skDn;

REAL*8, DIMENSION(:), ALLOCATABLE :: denom_m, denom_p
REAL*8, DIMENSION(:,:), ALLOCATABLE :: kxCoeffs
REAL*8 :: term1
REAL :: norm


skDn = skD/L
totPts = SIZE(vgx)

ALLOCATE(kxCoeffs(N+1,N+1), denom_m(totPts), denom_p(totPts), STAT = statvar)

kxCoeffs = 0; denom_m = 0; denom_p = 0;

! ###########################################################################
! find cosine expansion coefficients


norm = 2.0*SUM(C_mu/tau*weights)

DO m = 0,N    ! inhomogenous kernel f coefficients 
	denom_m = gamma_mu*(1.0 + DBLE(m)**2 * pi**2 / (gamma_mu**2))
	
	fx1Coeffs(m+1,:) = 0.0
	fx2Coeffs(m+1,:) = 2.0/norm/tau/denom_m*((EXP(-gamma_mu)-(-1.0)**m*EXP(-2.0*gamma_mu)) + ((-1.0)**m)-EXP(-gamma_mu))*weights
	fx2Coeffs(m+1,:) = 2.0/norm/tau/denom_m*()
	fx3Coeffs(m+1) = SUM(Q_mu/norm/Kn_mu*weights*(skDn/((skDn**2)*(gamma_mu**2)- 1.0) * &
						((EXP(-1.0/skDn-gamma_mu)*(skDn*gamma_mu-1) + &
						  EXP(-1.0/skDn)*((-1.0)**m)*(1.0-skDn*gamma_mu) + &
						  EXP(-gamma_mu)*((-1.0)**m)*(1.0+skDn*gamma_mu) - (1.0+skDn*gamma_mu))/denom_m + &
						2.0*gamma_mu*(skDn**2)*(1.0-(-1.0)**m*EXP(-1.0/skDn))/(1.0+(skDn*m*pi)**2)  )))	
END DO



DO p = 0,N
	denom_p = gamma_mu**2 + DBLE(p)**2*pi**2
	
	DO m = 0,N
		denom_m = gamma_mu**2 + DBLE(m)**2*pi**2
		
		IF (m==p) THEN        ! diagonal terms
			IF (p==0) THEN
				term1 = SUM(C_mu/Kn_mu/tau/denom_p/denom_m*weights*(2.0*(EXP(-gamma_mu) +gamma_mu -1.0)*gamma_mu**2))
			ELSE
				term1 = SUM(C_mu/Kn_mu/tau/denom_p/denom_m*weights* & 
								(EXP(-gamma_mu)*(2.0*((-1.0)**m)*gamma_mu**2) + &
								  gamma_mu*((gamma_mu-2.0)*gamma_mu + m**2.0*pi**2)))
			END IF
		ELSE
			term1 = SUM(C_mu/Kn_mu/tau/denom_p/denom_m*(gamma_mu**2)*weights* & 
								(EXP(-gamma_mu)*((-1.0)**m + (-1.0)**p) - (1.0+(-1.0)**(m+p))))
		END IF
		
		term2 = SUM(C_mu/Kn_mu/tau/denom_p/denom_m*(gamma_mu**2)*weights* &
								(1.0-((-1.0)**m)*EXP(-gamma_mu))*(1.0-((-1.0)**p)*EXP(-gamma_mu)))
		
		kxCoeffs(m+1,p+1) = 4.0/norm*(term1+term2);
		
	END DO
END DO


! phixb = fx1Coeffs*RESHAPE(C_mu*L*(norm/weights),(/N+1, totPts /),&
							! PAD = C_mu*L*(norm/weights), ORDER = (/ 2,1 /) );
! phixf = fx2Coeffs*RESHAPE(C_mu*L*(norm/weights),(/N+1, totPts /),&
							! PAD = C_mu*L*(norm/weights), ORDER = (/ 2,1 /) );
							
phixb = sDIAG_MATMUL_R(fx1Coeffs,C_mu*L*(norm/weights/2.0));
phixf = sDIAG_MATMUL_R(fx2Coeffs,C_mu*L*(norm/weights/2.0));
! phixb = cDIAG_MATMUL_R(fx1Coeffs,C_mu*L*(norm/2.0));
! phixf = cDIAG_MATMUL_R(fx2Coeffs,C_mu*L*(norm/2.0));
! includes vgx term


fx1Coeffs(1,:) = 1.0/2.0*fx1Coeffs(1,:);
fx2Coeffs(1,:) = 1.0/2.0*fx2Coeffs(1,:);
fx3Coeffs(1) = 1.0/2.0*fx3Coeffs(1);

! k1Coeffs(:,1) = 1/2*k1Coeffs(:,1)*2;
kxCoeffs(1,:) = 1.0/2.0*kxCoeffs(1,:);
! PRINT *, 'kxcoeffs', kxCoeffs(1:5,1:5)

! Amat = eye(N+1) - 1/2*k1Coeffs
Amat = -1.0/2.0*kxCoeffs;
DO i = 0,N
	Amat(i+1,i+1) = Amat(i+1,i+1) + 1;
END DO 

DEALLOCATE(kxCoeffs, denom_m, denom_p,STAT = statvar)

END SUBROUTINE