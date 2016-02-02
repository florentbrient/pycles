MODULE CppWrappers

USE func_pointers
USE ISO_C_BINDING

CONTAINS

  SUBROUTINE convert_c_funcpointer(cpointer,fpointer)
    USE ISO_C_BINDING
    TYPE(C_FUNPTR), VALUE :: cpointer
    PROCEDURE(), POINTER :: fpointer
    CALL C_F_PROCPOINTER(cpointer,fpointer)
  END SUBROUTINE convert_c_funcpointer

 SUBROUTINE allocate_Container(Container_cptr)
    TYPE (C_PTR) :: Container_cptr

    TYPE (Container), POINTER :: Container_fptr

    ALLOCATE( Container_fptr )
    Container_cptr = C_LOC(Container_fptr)
  END SUBROUTINE allocate_Container

 SUBROUTINE deallocate_Container(Container_cptr)
    TYPE (C_PTR), VALUE :: Container_cptr

    TYPE (Container), POINTER :: Container_fptr

    CALL C_F_POINTER(Container_cptr, Container_fptr)
    DEALLOCATE( Container_fptr )
  END SUBROUTINE deallocate_Container

END MODULE CppWrappers
