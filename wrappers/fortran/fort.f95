program test
USE, INTRINSIC :: ISO_C_BINDING
USE OPENST
use omp_lib
    REAL (C_DOUBLE), TARGET, ALLOCATABLE :: U(:,:,:), V(:,:,:)
    CHARACTER (C_CHAR), TARGET, ALLOCATABLE :: LSM_UNLOCKED(:,:,:)
    INTEGER (C_SIZE_T) :: NI , NJ, NK
    REAL (C_DOUBLE) :: HI, HJ, HK, SRCI, SRCJ, SRCK
    TYPE (C_PTR) :: U_PTR, V_PTR, LSM_UNLOCKED_PTR, converged_ptr, it_ptr
    REAL (C_DOUBLE) :: EPS, CHECK, di, dj, dk, dist
    REAL (C_DOUBLE) :: L1, L2, LINF, Umin, Umean, Umax
    INTEGER (C_INT) :: max_iter, errcode
    INTEGER (C_INT), TARGET :: it, converged 
    INTEGER (C_SIZE_T):: i,j,k
    double precision :: start, finish
    
    NI = 320
    NJ = 320
    NK = 320
    
    SRCI = 0.5
    SRCJ = 0.5
    SRCK = 0.5
    HI = 1.0 / (NI - 1)
    HJ = 1.0 / (NJ - 1)
    HK = 1.0 / (NK - 1)
    max_iter = 10
    EPS = 0.0
    
    ALLOCATE ( U (NI, NJ, NK) ) 
    ALLOCATE ( V (NI, NJ, NK) ) 
    ALLOCATE ( LSM_UNLOCKED (NI, NJ, NK) ) 
    
    U_PTR = C_LOC(U)
    V_PTR = C_LOC(V)
    LSM_UNLOCKED_PTR = C_LOC(LSM_UNLOCKED)
    it_ptr = C_LOC(it)
    converged_ptr = C_LOC(converged)
    
    do k = 1,NK
    do j = 1,NJ
    do i = 1,NI
    V(i,j,k) = 1.0
    end do
    end do
    end do
    
    start = omp_get_wtime()
    V = RESHAPE( V, [NI,NJ,NK], order=[3,2,1])
    finish = omp_get_wtime()
    print '("RESHAPE C INPUT: ",f12.6," seconds.")',finish-start
    
    start = omp_get_wtime()
    errcode = OpenST_LSM3D(U_PTR,LSM_UNLOCKED_PTR,V_PTR,NI,NJ,NK,HI,HJ,HK,SRCI,SRCJ,SRCK,EPS,max_iter,it_ptr,converged_ptr)
    finish = omp_get_wtime()
    print '("LSM3D: ",f12.6," seconds, ",i3," iterations, converged ",i1)',finish-start,it,converged
    
    if (errcode.NE.0) then
    print '("OPENST_ERROR: code ",i2)',errcode
    STOP 1
    end if
    
    start = omp_get_wtime()
    U = RESHAPE( U, [NI,NJ,NK], order=[1,2,3])
    finish = omp_get_wtime()
    print '("RESHAPE C OUTPUT: ",f12.6," seconds.")',finish-start
    call flush()
    
    it = 0
    L1 = 0.0
    L2 = 0.0
    Umean = 0.0
    do k = 1,NK
    do j = 1,NJ
    do i = 1,NI
    di = SRCI - DBLE(i - 1) * HI
    dj = SRCJ - DBLE(j - 1) * HJ
    dk = SRCK - DBLE(k - 1) * HK
    dist = SQRT(di * di + dj * dj + dk * dk)
    CHECK = ABS(U(i,j,k) - dist)
    L1 = L1 + CHECK
    L2 = L2 + CHECK * CHECK
    Umean = Umean + U(i,j,k)
    if (it > 0) then
        LINF = MAX(CHECK,LINF)
        Umin = MIN(U(i,j,k),Umin)
        Umax = MAX(U(i,j,k),Umax)
    else
        LINF = CHECK
        Umin = U(i,j,k)
        Umax = U(i,j,k)
    end if
    it = it + 1
    end do
    end do
    end do
    L1 = L1 / DBLE(NI * NJ * NK)
    L2 = L2 / DBLE(NI * NJ * NK)
    Umean = Umean / DBLE(NI * NJ * NK)
    
    print '(A, T10, E15.6)', "L1:", L1, "L2:", L2, "LINF:", LINF
    print '(A, T10, E15.6)', "Umin:", Umin, "Umean:", Umean, "Umax:", Umax

end program
