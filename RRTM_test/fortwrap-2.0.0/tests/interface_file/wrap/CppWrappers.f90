MODULE CppWrappers

USE source
USE ISO_C_BINDING

CONTAINS

 SUBROUTINE allocate_Object(Object_cptr)
    TYPE (C_PTR) :: Object_cptr

    TYPE (Object), POINTER :: Object_fptr

    ALLOCATE( Object_fptr )
    Object_cptr = C_LOC(Object_fptr)
  END SUBROUTINE allocate_Object

 SUBROUTINE deallocate_Object(Object_cptr)
    TYPE (C_PTR), VALUE :: Object_cptr

    TYPE (Object), POINTER :: Object_fptr

    CALL C_F_POINTER(Object_cptr, Object_fptr)
    DEALLOCATE( Object_fptr )
  END SUBROUTINE deallocate_Object

 SUBROUTINE allocate_Object_to_rename(Object_to_rename_cptr)
    TYPE (C_PTR) :: Object_to_rename_cptr

    TYPE (Object_to_rename), POINTER :: Object_to_rename_fptr

    ALLOCATE( Object_to_rename_fptr )
    Object_to_rename_cptr = C_LOC(Object_to_rename_fptr)
  END SUBROUTINE allocate_Object_to_rename

 SUBROUTINE deallocate_Object_to_rename(Object_to_rename_cptr)
    TYPE (C_PTR), VALUE :: Object_to_rename_cptr

    TYPE (Object_to_rename), POINTER :: Object_to_rename_fptr

    CALL C_F_POINTER(Object_to_rename_cptr, Object_to_rename_fptr)
    DEALLOCATE( Object_to_rename_fptr )
  END SUBROUTINE deallocate_Object_to_rename

END MODULE CppWrappers
