MODULE CppWrappers

USE DerivedTypes
USE ISO_C_BINDING

CONTAINS

 SUBROUTINE allocate_ObjectA(ObjectA_cptr)
    TYPE (C_PTR) :: ObjectA_cptr

    TYPE (ObjectA), POINTER :: ObjectA_fptr

    ALLOCATE( ObjectA_fptr )
    ObjectA_cptr = C_LOC(ObjectA_fptr)
  END SUBROUTINE allocate_ObjectA

 SUBROUTINE deallocate_ObjectA(ObjectA_cptr)
    TYPE (C_PTR), VALUE :: ObjectA_cptr

    TYPE (ObjectA), POINTER :: ObjectA_fptr

    CALL C_F_POINTER(ObjectA_cptr, ObjectA_fptr)
    DEALLOCATE( ObjectA_fptr )
  END SUBROUTINE deallocate_ObjectA

 SUBROUTINE allocate_ObjectB(ObjectB_cptr)
    TYPE (C_PTR) :: ObjectB_cptr

    TYPE (ObjectB), POINTER :: ObjectB_fptr

    ALLOCATE( ObjectB_fptr )
    ObjectB_cptr = C_LOC(ObjectB_fptr)
  END SUBROUTINE allocate_ObjectB

 SUBROUTINE deallocate_ObjectB(ObjectB_cptr)
    TYPE (C_PTR), VALUE :: ObjectB_cptr

    TYPE (ObjectB), POINTER :: ObjectB_fptr

    CALL C_F_POINTER(ObjectB_cptr, ObjectB_fptr)
    DEALLOCATE( ObjectB_fptr )
  END SUBROUTINE deallocate_ObjectB

END MODULE CppWrappers
