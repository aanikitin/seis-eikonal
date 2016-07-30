program test
USE, INTRINSIC :: ISO_C_BINDING
USE OPENST
use omp_lib
    REAL (C_DOUBLE), TARGET, ALLOCATABLE :: U(:,:,:), V(:,:,:)
    INTEGER (C_SIZE_T) :: NI , NJ, NK
    REAL (C_DOUBLE) :: HI, HJ, HK, SRCI, SRCJ, SRCK, RCVI, RCVJ, RCVK
    TYPE (C_PTR) :: U_PTR, V_PTR
    REAL (C_DOUBLE) :: di, dj, dk, dist
    INTEGER (C_INT) :: errcode
    INTEGER (C_SIZE_T):: i,j,k
    double precision :: start, finish
    
    REAL (C_DOUBLE) :: TSTEP
    REAL (C_DOUBLE), POINTER, DIMENSION(:,:) :: RAY
    TYPE (C_PTR), TARGET :: RAY_PTR, RAY_PTR_PTR
    INTEGER (C_SIZE_T), TARGET :: RAY_NI, RAY_NJ
    TYPE (C_PTR) :: RAY_NI_PTR, RAY_NJ_PTR
    
    NI = 100
    NJ = 100
    NK = 100
    
    RCVI = 0.1
    RCVJ = 0.1
    RCVK = 0.1
    SRCI = 0.5
    SRCJ = 0.5
    SRCK = 0.5
    HI = 1.0 / (NI - 1)
    HJ = 1.0 / (NJ - 1)
    HK = 1.0 / (NK - 1)
    TSTEP = MIN(HI,HJ,HK) / 1.0
    
    print '("TSTEP: ",E15.6)', TSTEP
    
    ALLOCATE ( U (NI, NJ, NK) )
    ALLOCATE ( V (NI, NJ, NK) )
    
    U_PTR = C_LOC(U)
    V_PTR = C_LOC(V)
    
    do k = 1,NK
    do j = 1,NJ
    do i = 1,NI
    V(i,j,k) = 1.0
    di = SRCI - DBLE(i - 1) * HI
    dj = SRCJ - DBLE(j - 1) * HJ
    dk = SRCK - DBLE(k - 1) * HK
    dist = SQRT(di * di + dj * dj + dk * dk)
    U(i,j,k) = dist
    end do
    end do
    end do
    
    RAY_PTR_PTR = C_LOC(RAY_PTR)
    RAY_NI_PTR = C_LOC(RAY_NI)
    RAY_NJ_PTR = C_LOC(RAY_NJ)
    
    start = omp_get_wtime()
    errcode = OpenST_BRT3D_Trace(U_PTR,V_PTR,NI,NJ,NK,HI,HJ,HK,TSTEP,RCVI,RCVJ,RCVK, &
    SRCI,SRCJ,SRCK,RAY_PTR_PTR,RAY_NI_PTR,RAY_NJ_PTR)
    finish = omp_get_wtime()
    print '("BRT3D: ",E15.6," seconds.")',finish-start
    
    if (errcode.NE.0) then
    print '("OPENST_ERROR: code ",i2)',errcode
    STOP 1
    end if
    
    CALL C_F_POINTER(RAY_PTR, RAY, [RAY_NI, RAY_NJ])
    RAY = RESHAPE( RAY, [RAY_NI,RAY_NJ], order=[2,1])
    
    do i = 1,RAY_NI
        print '(I3, E15.6, E15.6, E15.6)', i, RAY(i,2), RAY(i,2), RAY(i,3)
    end do
    
end program
