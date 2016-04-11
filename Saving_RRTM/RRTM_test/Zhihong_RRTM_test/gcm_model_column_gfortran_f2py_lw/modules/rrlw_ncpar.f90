module rrlw_ncpar
	

	implicit none
        save
	
        real(kind=8), parameter :: cpdair = 1003.5  ! Specific heat capacity of dry air
                                        		 ! at constant pressure at 273 K
                                        		 ! (J kg-1 K-1)

	
	integer(kind=4), parameter :: maxAbsorberNameLength = 5, &
                          Absorber              = 12
    character(len = maxAbsorberNameLength), dimension(Absorber), parameter :: &
    AbsorberNames = (/        &
     				'N2   ',  &
     				'CCL4 ',  &
     				'CFC11',  &
     				'CFC12',  &
     				'CFC22',  &
     				'H2O  ',  &
     				'CO2  ',  &
     				'O3   ',  &
     				'N2O  ',  & 
     				'CO   ',  &
     				'CH4  ',  &
     				'O2   '  /)
	
	integer(kind=4), dimension(40) :: status
	integer(kind=4) :: i
	integer(kind=4), parameter :: keylower  = 9,   &
						  keyupper  = 5,   &
						  Tdiff     = 5,   &
						  ps        = 59,  &
						  plower    = 13,  &
						  pupper    = 47,  &
						  Tself     = 10,  &
						  Tforeign  = 4,   &
						  pforeign  = 4,   &
						  T         = 19,  &
						  Tplanck   = 181, &
						  band      = 16,  &
						  GPoint    = 16,  &
						  GPointSet = 2
						  
	contains 
	
	subroutine getAbsorberIndex(AbsorberName,AbsorberIndex)
		character(len = *), intent(in) :: AbsorberName
		integer(kind=4), intent(out)           :: AbsorberIndex
		
		integer(kind=4) :: m
	
		AbsorberIndex = -1
		do m = 1, Absorber
			if (trim(AbsorberNames(m)) == trim(AbsorberName)) then
				AbsorberIndex = m
			end if
		end do
		
		if (AbsorberIndex == -1) then
			print*, "Absorber name index lookup failed."
		end if
	end subroutine getAbsorberIndex

end module rrlw_ncpar
