module mod_cu_kinds

#ifdef MPAS
   use mpas_kind_types, only: kind_phys     => RKIND
   use mpas_kind_types, only: kind_phys_dbl => R8KIND
#elif defined(CCPP)
   use machine,         only: kind_phys
   use machine,         only: kind_dbl_prec
#else
   use iso_fortran_env, only: kind_phys     => real32
   use iso_fortran_env, only: kind_phys_dbl => real64
#endif

   implicit none
   
 end module mod_cu_kinds
 
