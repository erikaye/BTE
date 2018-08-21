PROGRAM runCalcs_multTfct

USE omp_lib
USE Constants
IMPLICIT NONE

INTERFACE
	SUBROUTINE MultilayerT_BBC_fct_noMPI(fileNum,DFT_path,write_path, matlsInds,Ls,Q0,skD,numPer,dT,kCond)
		USE Constants
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: matlsInds(:), numPer
		CHARACTER(1), INTENT(IN) :: fileNum
		CHARACTER(longStr), INTENT(IN) :: DFT_path, write_path
		REAL, INTENT(IN) :: Q0, skD, Ls(:)
		REAL*8, DIMENSION(:,:), INTENT(OUT) :: dT
		REAL*8, DIMENSION(:) :: kCond
	END SUBROUTINE
	
	SUBROUTINE MultilayerT_PBC_fct(fileNum,DFT_path,write_path, matlsInds,Ls,Q0,skD,numPer,dT,kCond)
		USE Constants
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: matlsInds(:), numPer
		CHARACTER(1), INTENT(IN) :: fileNum
		CHARACTER(longStr), INTENT(IN) :: DFT_path, write_path
		REAL, INTENT(IN) :: Q0, skD, Ls(:)
		REAL*8, DIMENSION(:,:), INTENT(OUT) :: dT
		REAL*8, DIMENSION(:) :: kCond
	END SUBROUTINE
	
	SUBROUTINE MultilayerT_PBC_fct_constL(fileNum,DFT_path,write_path, matlsInds,Ls,Q0,skD,numPer,dT,kCond)
		USE Constants
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: matlsInds(:), numPer
		CHARACTER(1), INTENT(IN) :: fileNum
		CHARACTER(longStr), INTENT(IN) :: DFT_path, write_path
		REAL, INTENT(IN) :: Q0, skD, Ls(:)
		REAL*8, DIMENSION(:,:), INTENT(OUT) :: dT
		REAL*8, DIMENSION(:) :: kCond
	END SUBROUTINE
	
	SUBROUTINE MultilayerT_BBC_fct_noQ(fileNum,DFT_path,write_path, matlsInds,Ls,Q0,skD,numPer,dT,kCond)
		USE Constants
		IMPLICIT NONE
		
		INTEGER, INTENT(IN) :: matlsInds(:), numPer
		CHARACTER(1), INTENT(IN) :: fileNum
		CHARACTER(longStr), INTENT(IN) :: DFT_path, write_path
		REAL, INTENT(IN) :: Q0, skD, Ls(:)
		REAL*8, DIMENSION(:,:), INTENT(OUT) :: dT
		REAL*8, DIMENSION(:) :: kCond
	END SUBROUTINE
	
END INTERFACE


REAL :: skD, Q0!, roughness
INTEGER :: numPer, nL
INTEGER, DIMENSION(:), ALLOCATABLE :: matlsInds
REAL, DIMENSION(:), ALLOCATABLE :: Ls
REAL*8, DIMENSION(:,:), ALLOCATABLE :: dT!, intErr
REAL*8, DIMENSION(:), ALLOCATABLE :: kCond
CHARACTER(64) :: inFile
CHARACTER(1) :: fileNum
CHARACTER(longStr) :: DFT_path, write_path


WRITE(*,*) , 'num procs', OMP_GET_NUM_PROCS()
! CALL OMP_SET_DYNAMIC(.FALSE.)
CALL OMP_SET_NUM_THREADS(OMP_GET_NUM_PROCS())
PRINT *, 'stacksize', KMP_GET_STACKSIZE_S()
PRINT *, 'N', N

CALL GETARG(1,inFile)

Q0 = 0.0;
skD = 1.0e-10;



OPEN(UNIT=1,FILE=TRIM(inFile)//'.txt')
READ(1,*) nL, numPer
READ(1,*) 	! read next line

ALLOCATE(matlsInds(nL),Ls(nL))

READ(1,*) fileNum, Ls, matlsInds
READ(1,*)	! read next line (which contains descriptions)
READ(1,"(A)") DFT_path
READ(1,*)	! read next line (which contains descriptions)
READ(1,"(A)") write_path


PRINT *, 'T', T
PRINT *, 'Ls', Ls
PRINT *, 'matlsInds', matlsInds
PRINT *, 'numPer', numPer
PRINT *, 'fileNum ', fileNum
PRINT *, 'DFT path ', TRIM(DFT_path)
PRINT *, 'write path', TRIM(write_path)


IF (numPer > 10) THEN
	numPer = 1;
	PRINT *, 'periodic bc'
	
	ALLOCATE(dT(Nx,SIZE(Ls)), kCond(SIZE(Ls)))
	PRINT *, SHAPE(dT), SHAPE(kCond)
	PRINT *, Ls
	
	IF (Ls(1) /= Ls(2)) THEN
		CALL MultilayerT_PBC_fct_constL(fileNum,DFT_path,write_path,matlsInds,Ls,Q0,skD,numPer,dT,kCond)
	ELSE
		CALL MultilayerT_PBC_fct(fileNum,DFT_path,write_path,matlsInds,Ls,Q0,skD,numPer,dT,kCond)
	END IF
	
	
ELSE 
	ALLOCATE(dT(Nx,SIZE(Ls)*numPer), kCond(SIZE(Ls)*numPer-1))
	IF (numPer == 1) THEN
		CALL MultilayerT_BBC_fct_noMPI(fileNum,DFT_path,write_path,matlsInds,Ls,Q0,skD,numPer,dT,kCond)
	ELSE 
		CALL MultilayerT_BBC_fct_noQ(fileNum,DFT_path,write_path,matlsInds,Ls,Q0,skD,numPer,dT,kCond)
	END IF
END IF



PRINT*, 'back to main prog'

! DEALLOCATE(matlsInds,Ls,dT,kCond)
DEALLOCATE(matlsInds,Ls)

PRINT*, 'done!'


					
END PROGRAM runCalcs_multTfct