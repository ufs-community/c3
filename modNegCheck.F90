module modNegCheck

   use machine, only: kind_phys

   implicit none
   private

   public :: neg_check

contains

   !> Checks for excessive heating/moistening tendencies and prevents
   !! negative water vapor by scaling tendencies and precipitation together.
   !!
   !! Notes:
   !! - q is expected to be dry-air water vapor mixing ratio.
   !! - outq/outt/outu/outv/outqc are tendencies.
   !! - pret is scaled consistently with the tendencies.
   subroutine neg_check(name, j, dt, q, outq, outt, outu, outv, outqc, pret, &
                        its, ite, kts, kte, itf, ktf, ktop)

      implicit none

      character(len=*), intent(in) :: name

      integer, intent(in) :: j
      integer, intent(in) :: its, ite, kts, kte, itf, ktf
      integer, intent(in), dimension(its:) :: ktop

      real(kind=kind_phys), intent(in) :: dt

      real(kind=kind_phys), intent(inout), dimension(its:,kts:) :: q
      real(kind=kind_phys), intent(inout), dimension(its:,kts:) :: outq
      real(kind=kind_phys), intent(inout), dimension(its:,kts:) :: outt
      real(kind=kind_phys), intent(inout), dimension(its:,kts:) :: outu
      real(kind=kind_phys), intent(inout), dimension(its:,kts:) :: outv
      real(kind=kind_phys), intent(inout), dimension(its:,kts:) :: outqc
      real(kind=kind_phys), intent(inout), dimension(its:)      :: pret

      real(kind=kind_phys) :: names
      real(kind=kind_phys) :: scalef
      real(kind=kind_phys) :: thresh
      real(kind=kind_phys) :: qmem
      real(kind=kind_phys) :: qmem1
      real(kind=kind_phys) :: qmem2
      real(kind=kind_phys) :: qmemf
      real(kind=kind_phys) :: qtest

      integer :: i, k
      integer :: icheck

      !-------------------------------------------------------------------
      ! 1. Limit excessive heating/cooling tendencies.
      !-------------------------------------------------------------------

      thresh = 300.01_kind_phys
      names  = 1.0_kind_phys

      if (trim(name) == 'shallow' .or. trim(name) == 'mid') then
         thresh = 148.01_kind_phys
         names  = 1.0_kind_phys
      endif

      scalef = 86400.0_kind_phys

      !$acc kernels
      !$acc loop private(qmemf,qmem,qmem2,icheck)
      do i = its, itf

         if (ktop(i) <= 2) cycle

         icheck = 0
         qmemf  = 1.0_kind_phys

         !$acc loop reduction(min:qmemf)
         do k = kts, ktop(i)

            qmem = outt(i,k) * scalef

            if (qmem > thresh) then
               qmem2 = thresh / qmem
               qmemf = min(qmemf, qmem2)
               icheck = 1
            endif

            if (qmem < -0.5_kind_phys * thresh * names) then
               qmem2 = -0.5_kind_phys * names * thresh / qmem
               qmemf = min(qmemf, qmem2)
               icheck = 2
            endif

         enddo

         do k = kts, ktop(i)
            outq(i,k)  = outq(i,k)  * qmemf
            outt(i,k)  = outt(i,k)  * qmemf
            outu(i,k)  = outu(i,k)  * qmemf
            outv(i,k)  = outv(i,k)  * qmemf
            outqc(i,k) = outqc(i,k) * qmemf
         enddo

         pret(i) = pret(i) * qmemf

      enddo
      !$acc end kernels

      !-------------------------------------------------------------------
      ! 2. Prevent negative water vapor by scaling the same tendencies.
      !-------------------------------------------------------------------

      thresh = 1.0e-32_kind_phys

      !$acc kernels
      !$acc loop private(qmemf,qmem,qmem1,qmem2,qtest,icheck)
      do i = its, itf

         if (ktop(i) <= 2) cycle

         qmemf = 1.0_kind_phys

         !$acc loop reduction(min:qmemf)
         do k = kts, ktop(i)

            qmem = outq(i,k)

            if (abs(qmem) > 0.0_kind_phys .and. q(i,k) > 1.0e-6_kind_phys) then

               qtest = q(i,k) + outq(i,k) * dt

               if (qtest < thresh) then
                  qmem1 = abs(outq(i,k))
                  qmem2 = abs((thresh - q(i,k)) / dt)

                  if (qmem1 > 0.0_kind_phys) then
                     qmemf = min(qmemf, qmem2 / qmem1)
                     qmemf = max(0.0_kind_phys, qmemf)
                  endif
               endif

            endif

         enddo

         do k = kts, ktop(i)
            outq(i,k)  = outq(i,k)  * qmemf
            outt(i,k)  = outt(i,k)  * qmemf
            outu(i,k)  = outu(i,k)  * qmemf
            outv(i,k)  = outv(i,k)  * qmemf
            outqc(i,k) = outqc(i,k) * qmemf
         enddo

         pret(i) = pret(i) * qmemf

      enddo
      !$acc end kernels

   end subroutine neg_check

 end module modNegCheck
 
