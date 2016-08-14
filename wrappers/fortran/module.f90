    MODULE OPENST

    ! When using Intel Fortran use /heap-arrays option (Optimization/Heap Arrays)
    ! if encountering stack overflow in tests
    
    INTERFACE
    INTEGER(C_INT) FUNCTION OpenST_LSM3D(U,LSM_UNLOCKED,V, &
        NI,NJ,NK,HI,HJ,HK, &
        SRCI,SRCJ,SRCK, &
        EPS, max_iter, converged, it) &
        BIND(C, NAME='OpenST_LSM3D')
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
    TYPE (C_PTR), VALUE :: U
    TYPE (C_PTR), VALUE :: LSM_UNLOCKED
    TYPE (C_PTR), VALUE :: V
    INTEGER (C_SIZE_T), VALUE :: NI
    INTEGER (C_SIZE_T), VALUE :: NJ
    INTEGER (C_SIZE_T), VALUE :: NK
    REAL (C_DOUBLE), VALUE :: HI
    REAL (C_DOUBLE), VALUE :: HJ
    REAL (C_DOUBLE), VALUE :: HK
    REAL (C_DOUBLE), VALUE :: SRCI
    REAL (C_DOUBLE), VALUE :: SRCJ
    REAL (C_DOUBLE), VALUE :: SRCK
    REAL (C_DOUBLE), VALUE :: EPS
    INTEGER (C_INT), VALUE :: max_iter
    TYPE (C_PTR), VALUE :: it
    TYPE (C_PTR), VALUE :: converged
    END FUNCTION OpenST_LSM3D
    END INTERFACE

    INTERFACE
    INTEGER(C_INT) FUNCTION OpenST_BRT3D_Trace(T,V, &
        NI,NJ,NK,HI,HJ,HK,TSTEP, &
        RCVI,RCVJ,RCVK,SRCI,SRCJ,SRCK,MAX_SEG, &
        RAY,RAY_NI,RAY_NJ) &
        BIND(C, NAME='OpenST_BRT3D_Trace')
    use, intrinsic :: iso_c_binding
    IMPLICIT NONE
    TYPE (C_PTR), VALUE :: T
    TYPE (C_PTR), VALUE :: V
    INTEGER (C_SIZE_T), VALUE :: NI
    INTEGER (C_SIZE_T), VALUE :: NJ
    INTEGER (C_SIZE_T), VALUE :: NK
    REAL (C_DOUBLE), VALUE :: HI
    REAL (C_DOUBLE), VALUE :: HJ
    REAL (C_DOUBLE), VALUE :: HK
    REAL (C_DOUBLE), VALUE :: TSTEP
    REAL (C_DOUBLE), VALUE :: RCVI
    REAL (C_DOUBLE), VALUE :: RCVJ
    REAL (C_DOUBLE), VALUE :: RCVK
    REAL (C_DOUBLE), VALUE :: SRCI
    REAL (C_DOUBLE), VALUE :: SRCJ
    REAL (C_DOUBLE), VALUE :: SRCK
    INTEGER (C_SIZE_T), VALUE :: MAX_SEG
    TYPE (C_PTR), VALUE :: RAY
    TYPE (C_PTR), VALUE :: RAY_NI
    TYPE (C_PTR), VALUE :: RAY_NJ
    END FUNCTION OpenST_BRT3D_Trace
    END INTERFACE

    CONTAINS

    SUBROUTINE OpenST_Fortran_PrintArr(MAT)
    IMPLICIT NONE
    double precision, allocatable, intent(in) :: MAT(:,:,:)
    INTEGER i,j,k
    do i = lbound(MAT, 1), ubound(MAT, 1)
        write(*,*) 'i = ', i
        do j = lbound(MAT, 2), ubound(MAT, 2)
            do k = lbound(MAT, 3), ubound(MAT, 3)
                write(*, '(TL10, E10.3)', ADVANCE='NO') MAT(i,j,k)
            end do
            write(*,*)
        end do
        write(*,*)
    end do
    END SUBROUTINE

    END MODULE OPENST
