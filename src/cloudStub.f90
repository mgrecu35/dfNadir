module cldclass
  integer:: nlayer, nc
  integer:: nRhEofs, nCldwEofs
  real,allocatable :: temp(:), & ! temperature profile used in the radiative transfer module
       rlhm(:), &     ! mean relative humidity profile, currently, not used
       height(:), &   ! height
       press(:), &    ! pressure
       rlhmrc(:,:), & ! relative humidity profiles in precipitation regions, 
       rlhmcs(:,:), & ! relative humidity profiles in clear sky regions
       cc(:,:), &      ! cloud profiles
       rheofs(:,:), cldweofs(:,:), cldwm(:), stdrhPC(:), stdcldwPC(:)
  
 
  real, allocatable :: tpw(:), tpwR(:)
  real, allocatable :: atm_extKa(:,:), atm_extKag(:,:), cld_extKa(:,:), cld_extKag(:,:),   atm_extm(:), &
       cld_extm(:), atm_exts(:), cld_exts(:) 
                      ! atm_extKa(:,:), cld_extKa(:,:) water vapor and cloud extinction at Ka-band
                      ! atm_extm(:), cld_extm(:), atm_exts(:), cld_exts(:) 
                      ! were intended to hold statistics, i.e. mean, 
                      ! standard deviation, but ended up not used
  real,allocatable :: atm_extMw(:,:,:),  atm_extMwRg(:,:,:) 
! precomputed extinction at various heights and frequencies
  real,allocatable :: atm_extMwr(:,:,:) 
! precomputed extinction at various heights and frequencies
  real :: drrte       ! this is the layer thickness in the RT calculations, technically it is not needed
                      ! but used to expedite the calculation of the inclusion of radar derived 
                      ! electromagnetic properties into calculations

end module cldclass
