!>\file cu_c3_driver_ccpp.F90
!! CCPP-facing C3 driver that preserves the existing CCPP interface
!! while replacing the legacy deep/shallow entry points with a direct call to
!! the shared CUP_C3 core used by MPAS.
!!
!! This is first attempt does the following
!!   * keep the current CCPP argument list
!!   * keep the current CCPP setup / packing logic
!!   * insert a direct CUP_C3 call where the old shallow/deep path began
!!   * keep the existing CCPP postprocessing arrays / outputs as much as possible
!!


module cu_c3_driver_ccpp

   use mod_cu_kinds, only: kind_phys
   use modNegCheck, only: neg_check
   use progsigma, only: progsigma_calc
   use module_cu_c3, only: maxiens, deep, shal, mid, nmp,                    &
        icumulus_gf, cumulus_type, closure_choice, cum_entr_rate,            &
        cum_use_excess, cum_use_smooth_tend, cum_min_cloud_depth,            &
        cum_hei_down_land, cum_hei_down_ocean,                               &
        cum_hei_updf_land, cum_hei_updf_ocean,                               &
        cum_max_edt_land, cum_max_edt_ocean,                                 &
        cum_ave_layer, cum_t_star, cum_fr_min_entr,                          &
        c0_deep, c0_shal, c0_mid, autoconv, downdraft,                       &
        use_momentum_transp, use_tracer_transp, convection_tracer,           &
        use_sub3d, use_pass_cloudvol, liq_ice_number_conc,                   &
        Hcts, CUP_C3, calc_lcl, set_Tq_pertub, FractLiqF, fct1d3,            &
        initModConvParGF, modConvParGF_initialized

   implicit none

   private
   public :: cu_c3_driver_ccpp_run, progsigma_calc

contains

!> This is the Grell-Freitas convection scheme driver module.
!! \section arg_table_cu_c3_driver_ccpp_run Argument Table
!! \htmlinclude cu_c3_driver_ccpp_run.html
!!
!>\section gen_c3_driver Grell-Freitas Cumulus Scheme Driver General Algorithm
      subroutine cu_c3_driver_ccpp_run(mype,ntracer,garea,im,km,dt,flag_init,flag_restart, &
           do_ca,progsigma,cactiv,cactiv_m,g,cp,fv,r_d,xlv,r_v,tsfc,xlon,xlat, &
           tke_pbl,ten_t_pbl,ten_q_pbl,forcet,                              &
           forceqv_spechum,phil,prslk,delp,raincv,tmf,qmicro,sigmain,       &
           betascu,betamcu,betadcu,qv_spechum,t,cld1d,us,vs,t2di,omega,     &
           qv2di_spechum,p2di,psuri,                                        &
           hbot,htop,kcnv,xland,hfx2,qfx2,aod_gf,cliw,clcw,ca_deep,rainevap,&
           pbl,ud_mf,dd_mf,dt_mf,cnvw_moist,cnvc,imfshalcnv,                &
           flag_for_scnv_generic_tend,flag_for_dcnv_generic_tend,           &
           dtend,dtidx,ntqv,ntiw,ntcw,index_of_temperature,index_of_x_wind, &
           index_of_y_wind,index_of_process_scnv,index_of_process_dcnv,     &
           fhour,fh_dfi_radar,ix_dfi_radar,num_dfi_radar,cap_suppress,      &
           dfi_radar_max_intervals,ldiag3d,qci_conv,do_cap_suppress,        &
           sigmaout,maxupmf,maxMF,do_mynnedmf,ichoice_in,ichoicem_in,       &
           ichoice_s_in,ten_t,ten_u,ten_v,ten_q,dcliw,dclcw,errmsg,errflg)

      implicit none
      integer, parameter :: ensdim=16
      integer            :: imid_gf=1
      integer, parameter :: ideep=1
      integer            :: ichoice=0
      integer            :: ichoicem=13
      integer            :: ichoice_s=3
      integer            :: init_status
      
      logical, intent(in) :: do_cap_suppress
      real(kind=kind_phys), parameter :: aodc0=0.14_kind_phys
      real(kind=kind_phys), parameter :: aodreturn=30._kind_phys
      integer, parameter :: dicycle=0
      integer, parameter :: dicycle_m=0
      integer :: ishallow_g3

      integer :: its,ite,jts,jte,kts,kte
      integer, intent(in) :: im,km,ntracer,mype
      integer, intent(in) :: ichoice_in,ichoicem_in,ichoice_s_in
      logical, intent(in) :: flag_init, flag_restart, do_mynnedmf
      logical, intent(in) :: flag_for_scnv_generic_tend,flag_for_dcnv_generic_tend, do_ca
      real(kind=kind_phys), intent(in) :: g,cp,fv,r_d,xlv,r_v,betascu,betamcu,betadcu
      logical, intent(in) :: ldiag3d
      logical, intent(in) :: progsigma
      real(kind=kind_phys), intent(inout), optional :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:), index_of_x_wind, index_of_y_wind,     &
           index_of_temperature, index_of_process_scnv, index_of_process_dcnv,  &
           ntqv, ntcw, ntiw
      
      real(kind=kind_phys), dimension(:,:), intent(in), optional :: forcet, forceqv_spechum
      real(kind=kind_phys), dimension(:,:), intent(in) :: ten_t_pbl,tke_pbl
      real(kind=kind_phys), dimension(:,:), intent(in) :: ten_q_pbl
      real(kind=kind_phys), dimension(:,:), intent(in) :: omega, phil, delp, prslk
      real(kind=kind_phys), dimension(:,:), intent(in), optional :: sigmain, qmicro
      real(kind=kind_phys), dimension(:,:), intent(in) :: t, us, vs
      real(kind=kind_phys), dimension(:,:), intent(inout), optional :: qci_conv
      real(kind=kind_phys), dimension(:,:), intent(out) :: cnvw_moist, cnvc
      real(kind=kind_phys), dimension(:,:), intent(out), optional :: sigmaout
      real(kind=kind_phys), dimension(:,:), intent(in) :: cliw, clcw
      real(kind=kind_phys), dimension(:,:,:), intent(in) :: tmf

      real(kind=kind_phys), allocatable :: clcw_save(:,:), cliw_save(:,:)

      integer, intent(in) :: dfi_radar_max_intervals
      real(kind=kind_phys), intent(in) :: fhour, fh_dfi_radar(:)
      integer, intent(in) :: num_dfi_radar, ix_dfi_radar(:)
      real(kind=kind_phys), intent(in), optional :: cap_suppress(:,:)

      integer, dimension(:), intent(out) :: hbot, htop, kcnv
      integer, dimension(:), intent(in) :: xland
      real(kind=kind_phys), dimension(:), intent(in) :: pbl
      real(kind=kind_phys), dimension(:), intent(in), optional :: maxMF
      integer, dimension(im) :: tropics

      real(kind=kind_phys), dimension(:), intent(in) :: hfx2, qfx2, psuri, tsfc, xlat, xlon
      real(kind=kind_phys), dimension(:), intent(in), optional :: ca_deep
      real(kind=kind_phys), dimension(:,:), intent(out), optional :: ud_mf
      real(kind=kind_phys), dimension(:,:), intent(out) :: dd_mf, dt_mf
      real(kind=kind_phys), dimension(:), intent(out) :: raincv, cld1d, rainevap
      real(kind=kind_phys), dimension(:), intent(out), optional :: maxupmf
      real(kind=kind_phys), dimension(:,:), intent(in) :: t2di, p2di

      real(kind=kind_phys), dimension(:,:), intent(in) :: qv2di_spechum, qv_spechum
      real(kind=kind_phys), dimension(:), intent(inout), optional :: aod_gf
      real(kind=kind_phys), dimension(im,km) :: qv2di, qv, forceqv, cnvw

      real(kind=kind_phys), dimension(:), intent(in) :: garea
      real(kind=kind_phys), intent(in) :: dt

      integer, intent(in) :: imfshalcnv
      integer, dimension(:), intent(inout), optional :: cactiv, cactiv_m
      real(kind_phys), dimension(:,:), intent(out) :: ten_t, ten_u, ten_v, dcliw, dclcw
      real(kind_phys), dimension(:,:,:), intent(out) :: ten_q
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! --- Original CCPP local arrays retained ---
      integer, dimension(im) :: k22_shallow,kbcon_shallow,ktop_shallow
      real(kind=kind_phys), dimension(im) :: rand_mom,rand_vmas
      real(kind=kind_phys), dimension(im,4) :: rand_clos
      real(kind=kind_phys), dimension(im,km,11) :: gdc,gdc2
      real(kind=kind_phys), dimension(im) :: ht, ccn_gf, ccn_m, dx, frhm, frhd
      real(kind=kind_phys), dimension(im,km) :: outt,outq,outqc,phh,subm,cupclw,cupclws
      real(kind=kind_phys), dimension(im,km) :: dhdt,zu,zus,zd,phf,zum,zdm,outum,outvm
      real(kind=kind_phys), dimension(im,km) :: outts,outqs,outqcs,outu,outv,outus,outvs
      real(kind=kind_phys), dimension(im,km) :: outtm,outqm,outqcm,submm,cupclwm
      real(kind=kind_phys), dimension(im,km) :: cnvwt,cnvwts,cnvwtm
      real(kind=kind_phys), dimension(im,km) :: hco,hcdo,zdo,zdd,hcom,hcdom,zdom,tmfq
      real(kind=kind_phys), dimension(km) :: zh
      real(kind=kind_phys), dimension(im) :: tau_ecmwf,edt,edtm,edtd,ter11,aa0,xlandi
      real(kind=kind_phys), dimension(im) :: pret,prets,pretm,hexec
      real(kind=kind_phys), dimension(im,10) :: forcing,forcing2
      integer, dimension(im) :: kbcon,ktop,ierr,ierrs,ierrm,kpbli
      integer, dimension(im) :: k22s,kbcons,ktops,k22,jmin,jminm
      integer, dimension(im) :: kbconm,ktopm,k22m
      real(kind=kind_phys), dimension(im,km) :: rho_dryar
      integer, parameter :: ipn = 0
      real(kind=kind_phys), dimension(im,km) :: zo,t2d,q2d,po,p2d,rhoi,clw_ten,new_qv_spechum,new_cliw,new_clcw
      real(kind=kind_phys), dimension(im,km) :: tn,qo,tshall,qshall,dz8w,omeg
      real(kind=kind_phys), dimension(im) :: z1,psur,cuten,cutens,cutenm
      real(kind=kind_phys), dimension(im) :: umean,vmean,pmean,mc_thresh
      real(kind=kind_phys), dimension(im) :: xmbs,xmbs2,xmb,xmbm,xmb_dumm,mconv
      integer :: i,j,k,icldck,ipr,jpr,jpr_deep,ipr_deep,uidx,vidx,tidx,qidx
      integer :: itf,jtf,ktf,iss,jss,nbegin,nend,cliw_idx,clcw_idx
      integer :: high_resolution
      real(kind=kind_phys) :: clwtot,clwtot1,excess,tcrit,tscl_kf,dp,dq,sub_spread,subcenter
      real(kind=kind_phys) :: dsubclw,dsubclws,dsubclwm,dtime_max,ztm,ztq,hfm,qfm,rkbcon,rktop
      real(kind=kind_phys), dimension(km) :: massflx,trcflx_in1,clw_in1,po_cup
      real(kind=kind_phys), dimension(im) :: flux_tun,tun_rad_mid,tun_rad_shall,tun_rad_deep
      character(len=50) :: ierrc(im),ierrcm(im),ierrcs(im)
      real(kind=kind_phys), dimension(im) :: hfx,qfx
      real(kind=kind_phys) :: tem,tem1,tf,tcr,tcrf,psum,arg_deep,arg_mid,arg_shal
      real(kind=kind_phys) :: cliw_shal,clcw_shal,tem_shal, cliw_both, weight_sum
      real(kind=kind_phys) :: cliw_deep,clcw_deep,tem_deep, clcw_both
      real(kind=kind_phys) :: t_tend, qv_tend, u_tend, v_tend
      real(kind=kind_phys) :: gdc_cloud
      integer :: cliw_deep_idx, clcw_deep_idx, cliw_shal_idx, clcw_shal_idx
      real(kind=kind_phys) :: cap_suppress_j(im)
      integer :: itime, do_cap_suppress_here
      logical :: exit_func
      real(kind=kind_phys) :: ccnclean

      ! --- New shared-C3 arrays patterned after the extracted MPAS driver ---
      integer, parameter :: mtp = 1
      integer, parameter :: ON = 1, OFF = 0
      integer :: plume, ii_plume
      integer :: use_excess, use_smooth_tend
      real(kind=kind_phys) :: pten, pqen, paph, zrho, total_dp
      real(kind=kind_phys) :: pahfs, pqhfl, zkhvfl, pgeoh
      real(kind=kind_phys) :: hei_down_land, hei_down_ocean, hei_updf_land, hei_updf_ocean
      real(kind=kind_phys) :: max_edt_land, max_edt_ocean, fr_min_entr, T_star, ave_layer
      real(kind=kind_phys) :: c0, min_cloud_depth
      real(kind=kind_phys) :: FSCAV(mtp)
      real(kind=kind_phys), dimension(im) :: zws, ccn, ztexec, zqexec, cum_ztexec, cum_zqexec
      real(kind=kind_phys), dimension(im) :: zlcl_sfc, plcl_sfc, tlcl_sfc, tke_mean, wlpool_c3
      real(kind=kind_phys), dimension(im) :: col_sat, stochastic_sig, tsur
      integer, dimension(im) :: max_ktop
      real(kind=kind_phys), dimension(im) :: h_sfc_flux, le_sfc_flux
      real(kind=kind_phys), dimension(im) :: xlons, xlats
      real(kind=kind_phys), dimension(im) :: lightn_dens

      !-----------------------------------------------------------------------
      ! MPAS-style thermodynamic states reconstructed in CCPP driver
      ! (i,k) ordering consistent with CCPP physics
      !-----------------------------------------------------------------------
      real(kind=kind_phys), dimension(im,km) :: temp_old, temp_new
      real(kind=kind_phys), dimension(im,km) :: temp_new_BL, temp_new_ADV
      
      real(kind=kind_phys), dimension(im,km) :: qv_old, qv_new, ten_q_pbl_mr
      real(kind=kind_phys), dimension(im,km) :: qv_new_BL, qv_new_ADV
      
      ! CUP_C3 follows the MPAS wrapper orientation: (vertical, horizontal).
      ! Keep these separate from the legacy CCPP arrays, which are mostly (im,km).
      real(kind=kind_phys), dimension(km,im) :: c3_rho, c3_t, c3_q, c3_tn, c3_qo
      real(kind=kind_phys), dimension(km,im) :: c3_tn_bl, c3_qo_bl, c3_tn_adv, c3_qo_adv
      real(kind=kind_phys), dimension(km,im) :: c3_zo, c3_po, c3_piexner, c3_us, c3_vs
      real(kind=kind_phys), dimension(km,im) :: c3_dm2d, c3_dhdt, c3_buoyx, c3_cnvcf, c3_turb_len_scale
      real(kind=kind_phys), dimension(km,im,1) :: c3_omeg
      real(kind=kind_phys), dimension(nmp,km,im) :: mpqi, mpql, mpcf
      real(kind=kind_phys), dimension(nmp,km,im,maxiens) :: outmpqi, outmpql, outmpcf
      real(kind=kind_phys), dimension(mtp,km,im) :: se_chem
      real(kind=kind_phys), dimension(mtp,km,im,maxiens) :: out_chem
      real(kind=kind_phys), dimension(km,im,maxiens) :: outt_c3, outq_c3, outqc_c3, outu_c3, outv_c3
      real(kind=kind_phys), dimension(km,im,maxiens) :: subten_q_c3, subten_t_c3, subten_u_c3, subten_v_c3
      real(kind=kind_phys), dimension(km,im,maxiens) :: outnliq_c3, outnice_c3, outbuoy_c3
      integer, dimension(im,maxiens) :: ierr4d, jmin4d, klcl4d, k224d, kbcon4d, ktop4d, kstabi4d, kstabm4d
      real(kind=kind_phys), dimension(im,maxiens) :: cprr4d, xmb4d, edt4d, pwav4d, sigma4d
      real(kind=kind_phys), dimension(km,im,maxiens) :: pcup5d, up_massentr5d, up_massdetr5d
      real(kind=kind_phys), dimension(km,im,maxiens) :: dd_massentr5d, dd_massdetr5d, zup5d, zdn5d
      real(kind=kind_phys), dimension(km,im,maxiens) :: prup5d, prdn5d, clwup5d, tup5d, conv_cld_fr5d

      integer :: kstop

      parameter (tf=258.16_kind_phys, tcr=273.16_kind_phys, tcrf=1.0_kind_phys/(tcr-tf))
            
      errmsg = ''
      errflg = 0

      init_status = initModConvParGF()
      
      ten_t = 0.0_kind_phys
      ten_u = 0.0_kind_phys
      ten_v = 0.0_kind_phys
      ten_q = 0.0_kind_phys
      dcliw = 0.0_kind_phys
      dclcw = 0.0_kind_phys
      new_clcw = clcw
      new_cliw = cliw

      ichoice   = ichoice_in
      ichoicem  = ichoicem_in
      ichoice_s = ichoice_s_in

      ! --- retain the existing CCPP setup logic ---
      ! --- this section is intentionally close to the current cu_c3_driver_ccpp.F90 ---

     if(do_cap_suppress) then
!$acc serial
       do itime=1,num_dfi_radar
         if(ix_dfi_radar(itime)<1) cycle
         if(fhour<fh_dfi_radar(itime)) cycle
         if(fhour>=fh_dfi_radar(itime+1)) cycle
         exit
       enddo
!$acc end serial
     endif
     if(do_cap_suppress .and. itime<=num_dfi_radar) then
        do_cap_suppress_here = 1
!$acc kernels
        cap_suppress_j(:) = cap_suppress(:,itime)
!$acc end kernels
     else
        do_cap_suppress_here = 0
!$acc kernels
        cap_suppress_j(:) = 0
!$acc end kernels
     endif

     if(ldiag3d) then
       if(flag_for_dcnv_generic_tend) then
         cliw_deep_idx=0
         clcw_deep_idx=0
       else
         cliw_deep_idx=dtidx(100+ntiw,index_of_process_dcnv)
         clcw_deep_idx=dtidx(100+ntcw,index_of_process_dcnv)
       endif
       if(flag_for_scnv_generic_tend) then
         cliw_shal_idx=0
         clcw_shal_idx=0
       else
         cliw_shal_idx=dtidx(100+ntiw,index_of_process_scnv)
         clcw_shal_idx=dtidx(100+ntcw,index_of_process_scnv)
       endif
       if(cliw_deep_idx>=1 .or. clcw_deep_idx>=1 .or. &
            cliw_shal_idx>=1 .or.  clcw_shal_idx>=1) then
         allocate(clcw_save(im,km), cliw_save(im,km))
!$acc enter data create(clcw_save,cliw_save)
!$acc kernels
         clcw_save(:,:)=clcw(:,:)
         cliw_save(:,:)=cliw(:,:)
!$acc end kernels
       endif
     endif

! TODO these should be coming in from outside
!
!    cactiv(:)      = 0
     rand_mom(:)    = 0._kind_phys
     rand_vmas(:)   = 0._kind_phys
     rand_clos(:,:) = 0._kind_phys
     lightn_dens(:) = 0._kind_phys
     kcnv(:) = 0
     
!$acc end kernels
!
     its=1
     ite=im
     itf=ite
     jts=1
     jte=1
     jtf=jte
     kts=1
     kte=km
     ktf=kte-1
!$acc kernels
! 
     tropics(:)=0
!
!> - Set tuning constants for radiation coupling
!
     tun_rad_shall(:)=.012
     tun_rad_mid(:)=.15 !.02
     tun_rad_deep(:)=.3 !.065
     edt(:)=0._kind_phys
     edtm(:)=0._kind_phys
     edtd(:)=0._kind_phys
     zdd(:,:)=0._kind_phys
     flux_tun(:)=5.
! dx for scale awareness
!$acc end kernels

     if (imfshalcnv == 5) then
      ishallow_g3 = 1
     else
      ishallow_g3 = 0
     end if
     high_resolution=0
     subcenter=0._kind_phys

!$acc kernels
     ud_mf(:,:) =0._kind_phys
     dd_mf(:,:) =0._kind_phys
     dt_mf(:,:) =0._kind_phys
     tau_ecmwf(:)=0._kind_phys
!$acc end kernels

!$acc kernels
     ht(:)=phil(:,1)/g
!$acc loop private(zh)
     do i=its,ite
      cld1d(i)=0._kind_phys
      zo(i,:)=phil(i,:)/g
      dz8w(i,1)=zo(i,2)-zo(i,1)
      zh(1)=0._kind_phys
      kpbli(i)=2
      do k=kts+1,ktf
       dz8w(i,k)=zo(i,k+1)-zo(i,k)
      enddo
!$acc loop seq
      do k=kts+1,ktf
       zh(k)=zh(k-1)+dz8w(i,k-1)
       if(zh(k).gt.pbl(i))then
        kpbli(i)=max(2,k)
        exit
       endif
      enddo
     enddo
!$acc end kernels

!$acc kernels
     do i= its,itf
      forcing(i,:)=0._kind_phys
      forcing2(i,:)=0._kind_phys
      ccn_gf(i) = 0._kind_phys
      ccn_m(i) = 0._kind_phys

      ! set aod and ccn
      if (flag_init .and. .not.flag_restart) then
        aod_gf(i)=aodc0
      else
        if((cactiv(i).eq.0) .and. (cactiv_m(i).eq.0))then
          if(aodc0>aod_gf(i)) aod_gf(i)=aod_gf(i)+((aodc0-aod_gf(i))*(dt/(aodreturn*60)))
          if(aod_gf(i)>aodc0) aod_gf(i)=aodc0
        endif
      endif

      ccn_gf(i)=max(5., (aod_gf(i)/0.0027)**(1/0.640))
      ccn_m(i)=ccn_gf(i)

      ccnclean=max(5., (aodc0/0.0027)**(1/0.640))

      hbot(i)  =kte
      htop(i)  =kts
      raincv(i)=0._kind_phys
      xlandi(i)=real(xland(i))
!     if(abs(xlandi(i)-1.).le.1.e-3) tun_rad_shall(i)=.15
!     if(abs(xlandi(i)-1.).le.1.e-3) flux_tun(i)=1.5
     enddo
     do i= its,itf
      mconv(i)=0._kind_phys
     enddo
     do k=kts,kte
      do i= its,itf
       omeg(i,k)=0._kind_phys
       zu(i,k)=0._kind_phys
       zum(i,k)=0._kind_phys
       zus(i,k)=0._kind_phys
       zd(i,k)=0._kind_phys
       zdm(i,k)=0._kind_phys
      enddo
     enddo

     psur(:)=0.01*psuri(:)
     do i=its,itf
      ter11(i)=max(0.,ht(i))
     enddo
     do k=kts,kte
      do i=its,ite
       cnvw(i,k)=0._kind_phys
       cnvc(i,k)=0._kind_phys
       gdc(i,k,1)=0._kind_phys
       gdc(i,k,2)=0._kind_phys
       gdc(i,k,3)=0._kind_phys
       gdc(i,k,4)=0._kind_phys
       gdc(i,k,7)=0._kind_phys
       gdc(i,k,8)=0._kind_phys
       gdc(i,k,9)=0._kind_phys
       gdc(i,k,10)=0._kind_phys
       gdc2(i,k,1)=0._kind_phys
      enddo
     enddo

      do k=kts,kte
         do i=its,ite
            tmfq(i,k)=tmf(i,k,1)
         enddo
      enddo

     ierr(:)=-99
     ierrm(:)=-99
     ierrs(:)=-99
     cuten(:)=0._kind_phys
     cutenm(:)=0._kind_phys
     cutens(:)=0._kind_phys
!$acc end kernels
     ierrc(:)=" "
!$acc kernels
     

     kbcon(:)=0
     kbcons(:)=0
     kbconm(:)=0

     ktop(:)=0
     ktops(:)=0
     ktopm(:)=0

     xmb(:)=0._kind_phys
     xmb_dumm(:)=0._kind_phys
     xmbm(:)=0._kind_phys
     xmbs(:)=0._kind_phys
     xmbs2(:)=0._kind_phys

     k22s(:)=0
     k22m(:)=0
     k22(:)=0

     jmin(:)=0
     jminm(:)=0

     pret(:)=0._kind_phys
     prets(:)=0._kind_phys
     pretm(:)=0._kind_phys

     umean(:)=0._kind_phys
     vmean(:)=0._kind_phys
     pmean(:)=0._kind_phys

     cupclw(:,:)=0._kind_phys
     cupclwm(:,:)=0._kind_phys
     cupclws(:,:)=0._kind_phys

     cnvwt(:,:)=0._kind_phys
     cnvwts(:,:)=0._kind_phys
     cnvwtm(:,:)=0._kind_phys

     hco(:,:)=0._kind_phys
     hcom(:,:)=0._kind_phys
     hcdo(:,:)=0._kind_phys
     hcdom(:,:)=0._kind_phys

     outt(:,:)=0._kind_phys
     outts(:,:)=0._kind_phys
     outtm(:,:)=0._kind_phys

     outu(:,:)=0._kind_phys
     outus(:,:)=0._kind_phys
     outum(:,:)=0._kind_phys

     outv(:,:)=0._kind_phys
     outvs(:,:)=0._kind_phys
     outvm(:,:)=0._kind_phys

     outq(:,:)=0._kind_phys
     outqs(:,:)=0._kind_phys
     outqm(:,:)=0._kind_phys

     outqc(:,:)=0._kind_phys
     outqcs(:,:)=0._kind_phys
     outqcm(:,:)=0._kind_phys

     subm(:,:)=0._kind_phys
     dhdt(:,:)=0._kind_phys

     frhm(:)=0._kind_phys
     frhd(:)=0._kind_phys

     temp_old(:,:)     = 0._kind_phys
     temp_new(:,:)     = 0._kind_phys
     temp_new_BL(:,:)  = 0._kind_phys
     temp_new_ADV(:,:) = 0._kind_phys
     
     qv_old(:,:)       = 0._kind_phys
     qv_new(:,:)       = 0._kind_phys
     qv_new_BL(:,:)    = 0._kind_phys
     qv_new_ADV(:,:)   = 0._kind_phys

     !-----------------------------------------------------------------------
     ! LB: Convert FV3/CCPP specific humidity to dry-air mixing ratio.
     !
     ! CCPP/FV3 inputs:
     !   qv_spechum        = current specific humidity
     !   qv2di_spechum     = specific humidity after dynamics, before physics
     !   forceqv_spechum   = dynamics tendency of specific humidity
     !   ten_q_pbl        = PBL tendency of specific humidity
     !
     ! C3 expects water vapor as dry-air mixing ratio.
     ! For q = r/(1+r), r = q/(1-q).
     ! For tendencies, dr/dt = dq/dt / (1-q)^2.
     !-----------------------------------------------------------------------

     do i = its, itf
        do k = kts, ktf
           qv2di(i,k) = qv2di_spechum(i,k) / &
                (1.0_kind_phys - qv2di_spechum(i,k))
           
           qv(i,k) = qv_spechum(i,k) / &
                (1.0_kind_phys - qv_spechum(i,k))
           
           forceqv(i,k) = forceqv_spechum(i,k) / &
                (1.0_kind_phys - qv2di_spechum(i,k))**2
           
           ten_q_pbl_mr(i,k) = ten_q_pbl(i,k) / &
                (1.0_kind_phys - qv2di_spechum(i,k))**2
        end do
     end do
     
     !-----------------------------------------------------------------------
     ! LB: Reconstruct MPAS-style thermodynamic states for shared C3.
     !
     ! GFS/CCPP calling sequence:
     !   dynamics -> radiation -> surface -> PBL -> gravity wave drag -> convection
     !
     ! CCPP inputs:
     !   t2di          = state after dynamics, before physics
     !   t             = current state entering convection
     !   forcet        = dynamics temperature tendency
     !   forceqv       = dynamics water-vapor tendency, dry mixing-ratio basis
     !   ten_t_pbl    = PBL temperature tendency
     !   ten_q_pbl    = PBL water-vapor tendency
     !
     ! MPAS-style equivalents:
     !   temp_old      = previous timestep / before dynamics
     !   temp_new      = state including all tendencies up to convection
     !   temp_new_ADV  = dynamics-only state
     !   temp_new_BL   = PBL-only state reconstructed from temp_old
     !-----------------------------------------------------------------------
     
     do i = its, itf
        do k = kts, ktf
           
           ! Previous timestep / before dynamics
           temp_old(i,k) = t2di(i,k) - forcet(i,k) * dt
           qv_old(i,k)   = max(1.e-16_kind_phys,qv2di(i,k) - forceqv(i,k) * dt)
           
           ! Current state entering convection
           temp_new(i,k) = t(i,k)
           qv_new(i,k)   = max(1.e-16_kind_phys,qv(i,k))
           
           ! Dynamics-only state
           temp_new_ADV(i,k) = t2di(i,k)
           qv_new_ADV(i,k)   = max(1.e-16_kind_phys,qv2di(i,k))
           
           ! PBL-only state reconstructed from old state
           temp_new_BL(i,k) = temp_old(i,k) + ten_t_pbl(i,k) * dt
           qv_new_BL(i,k)   = max(1.e-16_kind_phys,qv_old(i,k)) + ten_q_pbl_mr(i,k) * dt
           
           ! Total moist-static-energy forcing 
           dhdt(i,k) = cp  * ((temp_new(i,k) - temp_old(i,k)) / dt) + &
                xlv * ((qv_new(i,k)   - qv_old(i,k))   / dt)
           
        end do
     end do
     
     do k = kts, ktf
        do i = its, itf

           ! Pressure in hPa
           p2d(i,k) = 0.01 * p2di(i,k)
           po(i,k)  = p2d(i,k)
           
           rhoi(i,k) = 100.*p2d(i,k) / &
                (r_d*(temp_new_ADV(i,k)*(1.+fv*qv_new_ADV(i,k))))
                      
           ! shallow/PBL state
           tshall(i,k) = temp_new_BL(i,k)
           qshall(i,k) = max(1.e-16,qv_new_BL(i,k))

        end do
     end do

     !Turbulent kinetic energy
     do i = its, itf
        tke_mean(i) = 0._kind_phys
        total_dp = 0._kind_phys
        do k = kts+1, kpbli(i)-1
           dp = -0.5*(po(i,k+1)-po(i,k-1))
           tke_mean(i) = tke_mean(i) + (tke_pbl(i,k)/2.0)*dp
           !TODO:TKE division by two above because how tke_pbl is saved in the MYNN scheme.
           !Need to make this less scheme dependent in the future
           total_dp = total_dp + dp
        end do
        tke_mean(i) = max(1.e-5, tke_mean(i)/(total_dp+1.e-6))
     end do
     
      
!$acc end kernels

!$acc kernels
     do i=its,itf
      do k=kts,kpbli(i)
         tshall(i,k)=t(i,k)
         qshall(i,k)=max(1.e-16,qv(i,k))
      enddo
     enddo
!
! converting hfx2 and qfx2 to w/m2
!    hfx=cp*rho*hfx2
!    qfx=xlv*qfx2
     do i=its,itf
      hfx(i)=hfx2(i)*cp*rhoi(i,1)
      qfx(i)=qfx2(i)*xlv*rhoi(i,1)
      dx(i) = sqrt(garea(i))
      mc_thresh(i)=3.25/dx(i)
     enddo    

     do i=its,itf
      do k=kts,kpbli(i)
       !tn(i,k)=t(i,k)
       qo(i,k)=max(1.e-16,qv(i,k))
      enddo
     enddo
     nbegin=0
     nend=0

!$acc loop collapse(2) independent private(dp)
     do k=  kts+1,ktf-1
      do i = its,itf
       if((p2d(i,1)-p2d(i,k)).gt.150.and.p2d(i,k).gt.300)then
         dp=-.5*(p2d(i,k+1)-p2d(i,k-1))
!$acc atomic
         umean(i)=umean(i)+us(i,k)*dp
!$acc atomic
         vmean(i)=vmean(i)+vs(i,k)*dp
!$acc atomic
         pmean(i)=pmean(i)+dp
       endif
      enddo
     enddo
     do i = its,itf
       psum=0._kind_phys
       do k=kts,ktf-3
        if (clcw(i,k) .gt. -999.0 .and. clcw(i,k+1) .gt. -999.0 )then
           dp=(p2d(i,k)-p2d(i,k+1))
           psum=psum+dp
           clwtot = cliw(i,k) + clcw(i,k)
           if(clwtot.lt.1.e-32)clwtot=0.
           forcing(i,7)=forcing(i,7)+clwtot*dp
        endif
       enddo
       if(psum.gt.0)forcing(i,7)=forcing(i,7)/psum
       forcing2(i,7)=forcing(i,7)
     enddo
     do k=kts,ktf-1
      do i = its,itf
        omeg(i,k)= omega(i,k)
      enddo
   enddo
   
     !-- Initialize output:
     ierr4d(:,:)       = 0
     jmin4d(:,:)       = kts
     klcl4d(:,:)       = kts
     k224d(:,:)        = kts
     kbcon4d(:,:)      = kts
     ktop4d(:,:)       = kts
     kstabi4d(:,:)     = kts
     kstabm4d(:,:)     = kts
     
     cprr4d(:,:)       = 0.0_kind_phys
     xmb4d(:,:)        = 0.0_kind_phys
     edt4d(:,:)        = 0.0_kind_phys
     pwav4d(:,:)       = 0.0_kind_phys
     sigma4d(:,:)      = 0.0_kind_phys
     
     pcup5d(:,:,:)        = 0.0_kind_phys
     up_massentr5d(:,:,:) = 0.0_kind_phys
     up_massdetr5d(:,:,:) = 0.0_kind_phys
     dd_massentr5d(:,:,:) = 0.0_kind_phys
     dd_massdetr5d(:,:,:) = 0.0_kind_phys
     zup5d(:,:,:)         = 0.0_kind_phys
     zdn5d(:,:,:)         = 0.0_kind_phys

     !Initialize input fields:
     FSCAV(:)          = 0.0_kind_phys     ! okay if tracer/scavenge inactive
     stochastic_sig(:) = 1.0_kind_phys     ! okay neutral/default
     col_sat(:)        = 0.0_kind_phys     ! CHECK: may need real column saturation/flag
     se_chem(:,:,:)    = 0.0_kind_phys     ! okay if chemistry inactive
     
     do i = its, itf
         
        ! --- surface pressure (Pa -> hPa) ---
        psur(i) = 0.01_kind_phys * psuri(i)
        
        ! --- terrain height ---
        ter11(i) = max(0.0_kind_phys, ht(i))
        
        ! --- land/ocean mask (IMPORTANT: MPAS-style mapping) ---
        ! MPAS: 0 = land, 1 = ocean
        if (xland(i) < 1.5_kind_phys) then
           xlandi(i) = 0.0_kind_phys   ! land
        else
           xlandi(i) = 1.0_kind_phys   ! ocean
        end if
        
        ! --- lat/lon (used if popilating sounding, should pass in) ---
        xlons(i) = xlon(i)
        xlats(i) = xlat(i)
        
        ! --- surface temperature ---
        tsur(i) = tsfc(i)
        
     end do
     
     ! Pack CCPP column arrays into the shared-C3 working arrays.
     ! CCPP/FV3 arrays are (im,km); CUP_C3 expects the MPAS wrapper orientation (km,im).
     do i=its,itf

        wlpool_c3(i) = 0.0_kind_phys 
        ! first level = half layer above ground
        c3_zo(kts,i) = ter11(i) + 0.5_kind_phys * dz8w(i,kts)
        
        do k=kts,kte
           c3_t(k,i)      = temp_old(i,k)
           c3_q(k,i)      = qv_old(i,k)  
           c3_tn(k,i)     = temp_new(i,k)
           c3_qo(k,i)     = qv_new(i,k) 
           c3_tn_bl(k,i)  = temp_new_BL(i,k)
           c3_qo_bl(k,i)  = qv_new_BL(i,k)
           c3_tn_adv(k,i) = temp_new_ADV(i,k)
           c3_qo_adv(k,i) = qv_new_ADV(i,k)
           c3_dhdt(k,i)   = dhdt(i,k)
           c3_rho(k,i)            = rhoi(i,k)
           c3_po(k,i)             = po(i,k)
           c3_piexner(k,i)        = prslk(i,k)
           c3_us(k,i)             = us(i,k)
           c3_vs(k,i)             = vs(i,k)
           c3_dm2d(k,i)           = rhoi(i,k)
           c3_buoyx(k,i)          = 0.0_kind_phys !Not used in current config
           c3_cnvcf(k,i)          = 0.0_kind_phys !Not used in current config
           c3_turb_len_scale(k,i) = 0.0_kind_phys !Not used in current config
           c3_omeg(k,i,1)         = omeg(i,k)
           mpqi(:,k,i)            = 0.0_kind_phys
           mpql(:,k,i)            = 0.0_kind_phys
           mpcf(:,k,i)            = 0.0_kind_phys
        end do   
     end do
     do i = its, itf
        ! Special solution for height above ground
        ! first level = half layer above ground
        c3_zo(kts,i) = ter11(i) + 0.5_kind_phys * dz8w(i,kts)
        ! build upward
        do k = kts+1, ktf
           c3_zo(k,i) = c3_zo(k-1,i) + &
                0.5_kind_phys * (dz8w(i,k-1) + dz8w(i,k))
        end do
     end do
     
     do i = its, itf
        ! --- max cloud top (~50 hPa level) ---
        max_ktop(i) = kts - 1 + minloc(abs(c3_po(kts:ktf,i) - 50.0_kind_phys), 1)
     enddo
     
     !Use MPAS driver logic to compute temperature humidity
     !excess of the convective source parcel (relative to the environment)

     do i = its, itf
        
         pten = temp_old(i,1)
         pqen = qv_old(i,1)
         paph = 100.0_kind_phys * psur(i)
         
         zrho = paph / (r_d * pten * &
              (1.0_kind_phys + fv*pqen))
         
         h_sfc_flux(i)  = hfx(i)
         le_sfc_flux(i) = qfx(i)
         
         ! CCPP hfx/qfx are positive upward in W m-2.
         ! MPAS uses pahfs = -sflux_t and pqhfl = -sflux_r in the W* calculation.
         ! pahfs is Wm-2 and pqhfl is kg m2 s-1
         pahfs = -hfx(i)
         pqhfl = -qfx2(i) * rhoi(i,1)   ! kg m-2 s-1 since qfx2 is kinematic kg/kg m/s
         
         zkhvfl = (pahfs / cp + 0.608_kind_phys * pten * pqhfl) / zrho
         
         ! First-layer depth estimate
         pgeoh = 2.0_kind_phys * (zo(i,1) - ter11(i)) * g
         
         zws(i) = max(0.0_kind_phys, &
              0.001_kind_phys - 1.5_kind_phys*0.41_kind_phys*zkhvfl*pgeoh/pten)
         
         call calc_lcl(pten, paph, pqen, tlcl_sfc(i), plcl_sfc(i), zlcl_sfc(i))
         zlcl_sfc(i) = max(zlcl_sfc(i), 0.0_kind_phys)
         
         ztexec(i) = 0.0_kind_phys
         zqexec(i) = 0.0_kind_phys
         
         if (zws(i) > tiny(pgeoh)) then
            zws(i) = 1.2_kind_phys * zws(i)**0.3333_kind_phys
            
            ztexec(i) = max(0.0_kind_phys, &
                 -1.5_kind_phys*pahfs/(zrho*zws(i)*cp))
            
            zqexec(i) = max(0.0_kind_phys, &
                 -1.5_kind_phys*pqhfl/(zrho*zws(i)))
         end if

         ! Shallow closure W* using PBL depth
         pgeoh = pbl(i) * g
         zws(i) = max(0.0_kind_phys, &
              0.001_kind_phys - 1.5_kind_phys*0.41_kind_phys*zkhvfl*pgeoh/pten)
         zws(i) = 1.2_kind_phys * zws(i)**0.3333_kind_phys
         
      end do

      ! Main plume loop starts here:
      
      do ii_plume = 1, maxiens
         
         if(ii_plume == 1) then
            plume = deep
            c0 = c0_deep
         endif
         if(ii_plume == 2) then
            plume = shal
            c0 = c0_shal
            if(icumulus_gf(shal) == 2 ) where(ierr4d(:,deep) == 0) ierr4d(:,shal) = -99
            !if(icumulus_gf(shal) == 2 ) c0 = 0.0
         endif
         if(ii_plume == 3) then
            plume = mid
            c0 = c0_mid
            if(icumulus_gf(mid) == 2 ) where(ierr4d(:,deep) == 0) ierr4d(:,mid) = -99
         endif
         
         if(icumulus_gf(plume) == OFF) cycle

         use_smooth_tend = cum_use_smooth_tend(plume)
         min_cloud_depth = cum_min_cloud_depth(plume)
         hei_down_land   = cum_hei_down_land(plume)
         hei_down_ocean  = cum_hei_down_ocean(plume)
         hei_updf_land   = cum_hei_updf_land(plume)
         hei_updf_ocean  = cum_hei_updf_ocean(plume)
         max_edt_land    = cum_max_edt_land(plume)
         max_edt_ocean   = cum_max_edt_ocean(plume)
         use_excess      = cum_use_excess(plume)
         ave_layer       = cum_ave_layer(plume)
         T_star          = cum_t_star(plume)
         fr_min_entr     = cum_fr_min_entr(plume)

         call set_Tq_pertub(use_excess,its,ite,itf,xlandi,ztexec,zqexec,cum_ztexec,cum_zqexec)
                           
         call CUP_C3(its,ite,kts,kte,itf,ktf,mtp,nmp,FSCAV,                          &
              cumulus_type(plume), closure_choice(plume),                            &
              cum_entr_rate(plume), cum_use_excess(plume),                           &
              ! --- new plume-dependent inputs ---                                   &
              use_smooth_tend, min_cloud_depth,                                      &
              hei_down_land, hei_down_ocean,                                         &
              hei_updf_land, hei_updf_ocean,                                         &
              max_edt_land, max_edt_ocean,c0,                                        &
              use_excess, ave_layer, T_star, fr_min_entr,                            &
              ! input data                                                           &
              dx, stochastic_sig, col_sat, tke_mean,                                 &
              wlpool_c3, max_ktop, dt, kpbli,                                        &
              cum_ztexec, cum_zqexec, ccn_gf,                                        &
              c3_rho, c3_omeg, c3_t, c3_q,                                           &
              ter11, h_sfc_flux, le_sfc_flux, zlcl_sfc, xlons, xlats,                &
              xlandi, c3_tn, c3_qo, c3_tn_bl, c3_qo_bl,                              &
              c3_tn_adv, c3_qo_adv, c3_zo, c3_po, c3_piexner,                        &
              tsur, psur, c3_us, c3_vs, c3_dm2d,                                     &
              se_chem, zws, c3_dhdt, c3_buoyx, c3_cnvcf,                             &
              c3_turb_len_scale, mpqi, mpql, mpcf,                                   &
              ! output data                                                          &
              outt_c3(:,:,plume), outq_c3(:,:,plume), outqc_c3(:,:,plume),           &
              outu_c3(:,:,plume), outv_c3(:,:,plume),                                &
              subten_q_c3(:,:,plume), subten_t_c3(:,:,plume),                        &
              subten_u_c3(:,:,plume), subten_v_c3(:,:,plume),                        &
              outnliq_c3(:,:,plume), outnice_c3(:,:,plume),                          &
              outbuoy_c3(:,:,plume),                                                &
              outmpqi(:,:,:,plume), outmpql(:,:,:,plume),                            &
              outmpcf(:,:,:,plume), out_chem(:,:,:,plume),                           &
              ! convective transport / plume diagnostics                            &
              ierr4d(:,plume), jmin4d(:,plume), klcl4d(:,plume),                     &
              k224d(:,plume), kbcon4d(:,plume), ktop4d(:,plume),                     &
              kstabi4d(:,plume), kstabm4d(:,plume),                                  &
              cprr4d(:,plume), xmb4d(:,plume), edt4d(:,plume),                       &
              pwav4d(:,plume), sigma4d(:,plume),                                     &
              pcup5d(:,:,plume), up_massentr5d(:,:,plume),                           &
              up_massdetr5d(:,:,plume), dd_massentr5d(:,:,plume),                    &
              dd_massdetr5d(:,:,plume), zup5d(:,:,plume), zdn5d(:,:,plume),          &
              prup5d(:,:,plume), prdn5d(:,:,plume), clwup5d(:,:,plume),              &
              tup5d(:,:,plume), conv_cld_fr5d(:,:,plume),                            &
              lightn_dens, cld1d, rainevap)
         
      
      end do !plume
      ! Map shared-C3 outputs, which are (k,i,plume), back to the legacy CCPP
      ! local arrays, which are (i,k).
      do i = its, itf
         do k = kts, kte
            
            ! deep
            outt(i,k)  = outt_c3(k,i,deep)
            outq(i,k)  = outq_c3(k,i,deep)
            outqc(i,k) = outqc_c3(k,i,deep)
            outu(i,k)  = outu_c3(k,i,deep)
            outv(i,k)  = outv_c3(k,i,deep)
            
            ! shallow
            outts(i,k)  = outt_c3(k,i,shal)
            outqs(i,k)  = outq_c3(k,i,shal)
            outqcs(i,k) = outqc_c3(k,i,shal)
            outus(i,k)  = outu_c3(k,i,shal)
            outvs(i,k)  = outv_c3(k,i,shal)
            
            ! mid / congestus
            outtm(i,k)  = outt_c3(k,i,mid)
            outqm(i,k)  = outq_c3(k,i,mid)
            outqcm(i,k) = outqc_c3(k,i,mid)
            outum(i,k)  = outu_c3(k,i,mid)
            outvm(i,k)  = outv_c3(k,i,mid)
            
         end do
      end do
      
      !-----------------------------------------------------------------------
      ! Map C3 plume-indexed diagnostics (k,i,plume) back to legacy CCPP arrays (i,k)
      !
      ! NOTE:
      ! - clwup5d, zup5d, zdn5d are returned per plume (deep/shallow/mid)
      ! - Legacy CCPP driver expects separate arrays for each plume type
      ! - We transpose (k,i) -> (i,k) here
      !
      ! ASSUMPTIONS (TO VERIFY):
      ! - clwup5d represents updraft condensate profile per plume
      ! - zup5d / zdn5d represent updraft/downdraft mass flux or velocity profiles
      ! - cnvwt and cupclw:
      !     * cnvwt  -> used to compute convective condensate tendencies (scaled by xmb*dt)
      !     * cupclw -> used for diagnostics/radiation coupling (e.g., gdc)
      !   For now, both are mapped to clwup5d pending further validation
      !-----------------------------------------------------------------------

      do i = its, itf
         do k = kts, ktf
            
            ! Convective condensate profile (per plume)
            ! Used later to compute cnvw = cnvwt * xmb * dt
            cnvwt(i,k)  = clwup5d(k,i,deep)
            cnvwts(i,k) = clwup5d(k,i,shal)
            cnvwtm(i,k) = clwup5d(k,i,mid)
            
            ! Convective cloud water for diagnostics/radiation (gdc)
            ! Currently assumed same as condensate profile
            cupclw(i,k)  = clwup5d(k,i,deep)
            cupclws(i,k) = clwup5d(k,i,shal)
            cupclwm(i,k) = clwup5d(k,i,mid)
            
            ! Updraft profiles (mass flux or velocity depending on C3 definition)
            zu(i,k)  = zup5d(k,i,deep)
            zus(i,k) = zup5d(k,i,shal)
            zum(i,k) = zup5d(k,i,mid)
            
            ! Downdraft profiles
            ! Note: no shallow downdraft in legacy driver
            zd(i,k)  = zdn5d(k,i,deep)
            zdm(i,k) = zdn5d(k,i,mid)
            
         end do
      end do
      
      !-----------------------------------------------------------------------
      ! Map plume-column diagnostics from C3 back to legacy CCPP arrays.
      !
      ! C3 stores these as plume-indexed arrays:
      !   *_4d(i,plume)
      !
      ! Legacy CCPP driver expects separate arrays for:
      !   deep    : kbcon,  ktop,  ierr,  xmb,  edt,  pret
      !   shallow : kbcons, ktops, ierrs, xmbs,       prets
      !   mid     : kbconm, ktopm, ierrm, xmbm, edtm, pretm
      !
      ! TODO (Lisa): confirm whether final precipitation should use
      ! individual plume contributions or summed cprr4d over all plumes.
      !-----------------------------------------------------------------------
      
      do i = its, itf

          ! Deep default
         kbcon(i) = 0
         ktop(i)  = 0
         ierr(i)  = ierr4d(i,deep)
         xmb(i)   = 0.0_kind_phys
         edt(i)   = 0.0_kind_phys
         pret(i)  = 0.0_kind_phys
         ! Deep plume
         if(ierr4d (i,deep )==0 ) then
            kbcon(i) = kbcon4d(i,deep)
            ktop(i)  = ktop4d(i,deep)
            ierr(i)  = ierr4d(i,deep)
            xmb(i)   = xmb4d(i,deep)
            edt(i)   = edt4d(i,deep)
            pret(i)  = cprr4d(i,deep)
         endif

         !Shallow default
         kbcons(i) = 0
         ktops(i)  = 0
         ierrs(i)  = ierr4d(i,shal)
         xmbs(i)   = 0.0_kind_phys
         edtm(i)   = 0.0_kind_phys
         prets(i)  = 0.0_kind_phys
         ! Shallow plume
         if(ierr4d (i,shal )==0 ) then
            kbcons(i) = kbcon4d(i,shal)
            ktops(i)  = ktop4d(i,shal)
            ierrs(i)  = ierr4d(i,shal)
            xmbs(i)   = xmb4d(i,shal)
            edtm(i)   = edt4d(i,shal)
            prets(i)  = cprr4d(i,shal)
         endif

         !Mid default
          kbconm(i) = 0
          ktopm(i)  = 0
          ierrm(i)  = ierr4d(i,mid)
          xmbm(i)   = 0.0_kind_phys
          edtm(i)   = 0.0_kind_phys
          pretm(i)  = 0.0_kind_phys
         ! Mid/congestus plume
         if(ierr4d (i,mid )==0 ) then
            kbconm(i) = kbcon4d(i,mid)
            ktopm(i)  = ktop4d(i,mid)
            ierrm(i)  = ierr4d(i,mid)
            xmbm(i)   = xmb4d(i,mid)
            edtm(i)   = edt4d(i,mid)
            pretm(i)  = cprr4d(i,mid)
         endif
      end do

      !LB: This is a bit of a hard-coded bandaid and could be done better:
      !Need to understand why a MYNN MF results in so much deep CU.
      do i = its,itf
         if(mconv(i).lt.0.)mconv(i)=0.
         if(do_mynnedmf) then
            if((dx(i)<6500.).and.(maxMF(i).gt.0.))ierr(i)=555
         endif
      enddo
      
      if (dx(its)<6500.) then
         imid_gf=0
      endif
      
      ! LB TODO:
      !  * add chemistry and microphysics coupling 



      !This post processing code is from the original
      !CCPP cu_c3_driver_ccpp.F90 slightly cleaned up here:

      do i = its, itf
         
         kcnv(i) = 0
         
         ! Mid/congestus
         if (pretm(i) > 0.0_kind_phys .and. ierrm(i) == 0) then
            cutenm(i) = 1.0_kind_phys
            kcnv(i)   = 1
         else
            kbconm(i) = 0
            ktopm(i)  = 0
            xmbm(i)   = 0.0_kind_phys
            pretm(i)  = 0.0_kind_phys
            cutenm(i) = 0.0_kind_phys
         endif
         
         ! Deep takes precedence over mid, as in the legacy driver
         if (pret(i) > 0.0_kind_phys .and. ierr(i) == 0) then
            cuten(i)  = 1.0_kind_phys
            cutenm(i) = 0.0_kind_phys
            pretm(i)  = 0.0_kind_phys
            xmbm(i)   = 0.0_kind_phys
            kcnv(i)   = 1
            ktopm(i)  = 0
            kbconm(i) = 0
         else
            kbcon(i) = 0
            ktop(i)  = 0
            xmb(i)   = 0.0_kind_phys
            pret(i)  = 0.0_kind_phys
            cuten(i) = 0.0_kind_phys
         endif
         
         ! Shallow
         if (prets(i) > 0.0_kind_phys .and. ierrs(i) == 0) then
            cutens(i) = 1.0_kind_phys
            kcnv(i)   = 1
         else
            kbcons(i) = 0
            ktops(i)  = 0
            xmbs(i)   = 0.0_kind_phys
            prets(i)  = 0.0_kind_phys
            cutens(i) = 0.0_kind_phys
         endif
         
      end do

      ! Limit excessive tendencies and prevent negative qv before applying updates.
      if (icumulus_gf(deep) /= OFF) then
         call neg_check('deep', 1, dt, qv, outq, outt, outu, outv, outqc, pret, &
              its, ite, kts, kte, itf, ktf, ktop)
      endif
      
      if (icumulus_gf(shal) /= OFF) then
         call neg_check('shallow', 1, dt, qv, outqs, outts, outus, outvs, outqcs, prets, &
              its, ite, kts, kte, itf, ktf, ktops)
      endif
      
      if (icumulus_gf(mid) /= OFF) then
         call neg_check('mid', 1, dt, qv, outqm, outtm, outum, outvm, outqcm, pretm, &
              its, ite, kts, kte, itf, ktf, ktopm)
      endif
      
      do i=its,itf
         massflx(:)=0.
         trcflx_in1(:)=0.
         clw_in1(:)=0.
         do k=kts,ktf
            clw_ten(i, k)=0.
         enddo
         po_cup(:)=0.
         kstop=kts
         if(ktopm(i).gt.kts .or. ktop(i).gt.kts)kstop=max(ktopm(i),ktop(i))
         if(ktops(i).gt.kts)kstop=max(kstop,ktops(i))
         kstop = min(kstop, ktf)
         if(kstop.gt.2)then
            htop(i)=kstop
            if(kbcon(i).gt.2 .or. kbconm(i).gt.2)then
               hbot(i)=max(kbconm(i),kbcon(i)) !jmin(i)
            endif

            dtime_max=dt
            forcing2(i,3)=0.
            do k=kts,kstop
               arg_deep = 1.0_kind_phys
               arg_mid  = 1.0_kind_phys
               arg_shal = 1.0_kind_phys
               t_tend  = 0.0_kind_phys
               qv_tend = 0.0_kind_phys
               u_tend  = 0.0_kind_phys
               v_tend  = 0.0_kind_phys
               gdc_cloud = 0.0_kind_phys
               gdc(i,k,1)=0._kind_phys
               gdc(i,k,2)=0._kind_phys
               gdc(i,k,3)=0._kind_phys
               gdc(i,k,4)=0._kind_phys
               gdc(i,k,7)=0._kind_phys
               gdc(i,k,8)=0._kind_phys
               gdc(i,k,9)=0._kind_phys
               gdc(i,k,10)=0._kind_phys
               gdc2(i,k,1)=0._kind_phys
               
               if (xmb(i)  > 0.0_kind_phys) arg_deep = max(1.0_kind_phys, 1.0_kind_phys + 675.0_kind_phys*zu(i,k)*xmb(i))
               if (xmbm(i) > 0.0_kind_phys) arg_mid  = max(1.0_kind_phys, 1.0_kind_phys + 675.0_kind_phys*zum(i,k)*xmbm(i))
               if (xmbs(i) > 0.0_kind_phys) arg_shal = max(1.0_kind_phys, 1.0_kind_phys + 675.0_kind_phys*zus(i,k)*xmbs(i))
               
               cnvc(i,k) = 0.04_kind_phys * (log(arg_deep) + log(arg_mid) + log(arg_shal))
               cnvc(i,k) = min(cnvc(i,k), 0.6_kind_phys)
               cnvc(i,k) = max(cnvc(i,k), 0.0_kind_phys)

               if (xmb(i) > 0.0_kind_phys) then
                  cnvw(i,k) = cnvw(i,k) + cnvwt(i,k) * xmb(i) * dt
                  ud_mf(i,k)=cuten(i)*zu(i,k)*xmb(i)*dt
                  dd_mf(i,k)=cuten(i)*zd(i,k)*edt(i)*xmb(i)*dt
               endif
               
               if (xmbs(i) > 0.0_kind_phys) then
                  cnvw(i,k) = cnvw(i,k) + cnvwts(i,k) * xmbs(i) * dt
               endif
               
               if (xmbm(i) > 0.0_kind_phys) then
                  cnvw(i,k) = cnvw(i,k) + cnvwtm(i,k) * xmbm(i) * dt
               endif

               if (cutens(i) > 0.0_kind_phys) then
                  t_tend  = t_tend  + outts(i,k)
                  qv_tend = qv_tend + outqs(i,k)
                  u_tend  = u_tend  + outus(i,k)
                  v_tend  = v_tend  + outvs(i,k)
                  gdc_cloud = gdc_cloud + tun_rad_shall(i) * cupclws(i,k)
                  gdc(i,k,1) = max(0.0_kind_phys, tun_rad_shall(i) * cupclws(i,k))
                  gdc(i,k,4) = outts(i,k) * 86400.0_kind_phys
                  gdc(i,k,8) = gdc(i,k,8) + outqs(i,k)
               endif
               
               if (cutenm(i) > 0.0_kind_phys) then
                  t_tend  = t_tend  + outtm(i,k)
                  qv_tend = qv_tend + outqm(i,k)
                  u_tend  = u_tend  + outum(i,k)
                  v_tend  = v_tend  + outvm(i,k)
                  gdc_cloud = gdc_cloud + tun_rad_mid(i) * cupclwm(i,k)
                  gdc(i,k,3) = outtm(i,k) * 86400.0_kind_phys
                  gdc(i,k,8) = gdc(i,k,8) + outqm(i,k)
               endif

               if (cuten(i) > 0.0_kind_phys) then
                  t_tend  = t_tend  + outt(i,k)
                  qv_tend = qv_tend + outq(i,k)
                  u_tend  = u_tend  + outu(i,k)
                  v_tend  = v_tend  + outv(i,k)
                  gdc_cloud = gdc_cloud + tun_rad_deep(i) * cupclw(i,k)
                  gdc(i,k,2) = outt(i,k) * 86400.0_kind_phys
                  gdc(i,k,8) = gdc(i,k,8) + outq(i,k)
               endif

               ten_t(i,k) = t_tend
               ten_u(i,k) = u_tend
               ten_v(i,k) = v_tend
               qv(i,k) = qv(i,k) + dt * qv_tend

               gdc(i,k,8) = gdc(i,k,8) * 86400.0_kind_phys * xlv / cp
               gdc(i,k,9) = gdc(i,k,2) + gdc(i,k,3) + gdc(i,k,4)
               gdc(i,k,7) = -(gdc(i,k,7) - sqrt(us(i,k)**2 + vs(i,k)**2)) / dt

               gdc2(i,k,1) = max(0.0_kind_phys, gdc_cloud)
               qci_conv(i,k) = gdc2(i,k,1)
               
               !> - FCT treats subsidence effect to cloud ice/water (begin)
               dp = 100.0_kind_phys * (p2d(i,k) - p2d(i,k+1))
               dtime_max = min(dtime_max, 0.5_kind_phys * dp)
               po_cup(k) = 0.5_kind_phys * (p2d(i,k) + p2d(i,k+1))
               
               if (clcw(i,k) > -999.0_kind_phys .and. clcw(i,k+1) > -999.0_kind_phys) then
                  
                  clwtot  = cliw(i,k)   + clcw(i,k)
                  clwtot1 = cliw(i,k+1) + clcw(i,k+1)
                  
                  if (clwtot  < 1.0e-32_kind_phys) clwtot  = 0.0_kind_phys
                  if (clwtot1 < 1.0e-32_kind_phys) clwtot1 = 0.0_kind_phys
                  
                  clw_in1(k) = clwtot
                  
                  massflx(k) = 0.0_kind_phys
                  
                  if (xmb(i) > 0.0_kind_phys) then
                     massflx(k) = massflx(k) - xmb(i) * (zu(i,k) - edt(i)*zd(i,k))
                  endif
                  
                  if (xmbm(i) > 0.0_kind_phys) then
                     massflx(k) = massflx(k) - xmbm(i) * (zum(i,k) - edtm(i)*zdm(i,k))
                  endif
                  
                  if (xmbs(i) > 0.0_kind_phys) then
                     massflx(k) = massflx(k) - xmbs(i) * zus(i,k)
                  endif
                  
                  trcflx_in1(k) = massflx(k) * 0.5_kind_phys * (clwtot + clwtot1)
                  forcing2(i,3) = forcing2(i,3) + clwtot
                  
               endif
               
            enddo
            
            massflx(1)    = 0.0_kind_phys
            trcflx_in1(1) = 0.0_kind_phys
            
            call fct1d3(kstop, kte, dtime_max, po_cup, &
                 clw_in1, massflx, trcflx_in1, clw_ten(i,:))
            
            do k = kts, kstop
               
               tem = clw_ten(i,k)
               
               if (cutens(i) > 0.0_kind_phys) tem = tem + outqcs(i,k)
               if (cuten(i)  > 0.0_kind_phys) tem = tem + outqc(i,k)
               if (cutenm(i) > 0.0_kind_phys) tem = tem + outqcm(i,k)
               
               tem = dt * tem
               
               if (tem /= tem) tem = 0.0_kind_phys
               
               tem1 = max(0.0_kind_phys, min(1.0_kind_phys, (tcr - t(i,k))*tcrf))
               
               new_cliw(i,k) = 0.0_kind_phys
               new_clcw(i,k) = 0.0_kind_phys
               dcliw(i,k)    = 0.0_kind_phys
               dclcw(i,k)    = 0.0_kind_phys
               
               if (clcw(i,k) > -999.0_kind_phys) then
                  new_cliw(i,k) = max(0.0_kind_phys, cliw(i,k) + tem * tem1)
                  new_clcw(i,k) = max(0.0_kind_phys, clcw(i,k) + tem * (1.0_kind_phys - tem1))
                  dcliw(i,k) = (new_cliw(i,k) - cliw(i,k)) / dt
                  dclcw(i,k) = (new_clcw(i,k) - clcw(i,k)) / dt
               else
                  new_cliw(i,k) = max(0.0_kind_phys, cliw(i,k) + tem)
                  dcliw(i,k) = (new_cliw(i,k) - cliw(i,k)) / dt
               endif
               
            enddo

            !LB: legacy code:
            !gdc(i,1,10)=forcing(i,1)
            !gdc(i,2,10)=forcing(i,2)
            !gdc(i,3,10)=forcing(i,3)
            !gdc(i,4,10)=forcing(i,4)
            !gdc(i,5,10)=forcing(i,5)
            !gdc(i,6,10)=forcing(i,6)
            !gdc(i,7,10)=forcing(i,7)
            !gdc(i,8,10)=forcing(i,8)
            !gdc(i,10,10)=xmb(i)
            !gdc(i,11,10)=xmbm(i)
            !gdc(i,12,10)=xmbs(i)
            !gdc(i,13,10)=hfx(i)
            !gdc(i,15,10)=qfx(i)
            !gdc(i,16,10)=pret(i)*3600.

            !LB: forcing arrays not populated
            maxupmf(i)=0.
            if(forcing2(i,6).gt.0.)then
               maxupmf(i)=maxval(xmb(i)*zu(i,kts:ktf)/forcing2(i,6))
            endif
            
            if(ktop(i).gt.2 .and.pret(i).gt.0.)dt_mf(i,ktop(i)-1)=ud_mf(i,ktop(i))
         endif
      enddo

!$acc kernels
      do i = its, itf

         raincv(i)   = 0.0_kind_phys
         cactiv(i)   = 0
         cactiv_m(i) = 0
         
         ! Old logic: deep activates total convective precip.
         if (cuten(i) > 0.0_kind_phys .and. pret(i) > 0.0_kind_phys) then
            
            raincv(i) = raincv(i) + 0.001_kind_phys * pret(i) * dt
            cactiv(i) = 1
            
            ! In old code, shallow contribution was included only in the deep branch.
            if (cutens(i) > 0.0_kind_phys .and. prets(i) > 0.0_kind_phys) then
               raincv(i) = raincv(i) + 0.001_kind_phys * prets(i) * dt
            endif
            
            if (cutenm(i) > 0.0_kind_phys .and. pretm(i) > 0.0_kind_phys) then
               raincv(i) = raincv(i) + 0.001_kind_phys * pretm(i) * dt
            endif
            
         else
            
            ! Old fallback: mid only if deep is absent.
            if (cutenm(i) > 0.0_kind_phys .and. pretm(i) > 0.0_kind_phys) then
               raincv(i) = raincv(i) + 0.001_kind_phys * pretm(i) * dt
            endif
            
         endif

         if (cutenm(i) > 0.0_kind_phys .and. pretm(i) > 0.0_kind_phys) then
            cactiv_m(i) = 1
         endif
         
         ! Unify CCN
         if (ccn_m(i) == ccn_m(i) .and. ccn_gf(i) == ccn_gf(i)) then
            if (ccn_m(i) < ccn_gf(i)) ccn_gf(i) = ccn_m(i)
         endif
         
         if (ccn_gf(i) /= ccn_gf(i) .or. ccn_gf(i) < 0.0_kind_phys) then
            ccn_gf(i) = 0.0_kind_phys
         endif
         
         ! Convert CCN back to AOD
         aod_gf(i) = 0.0027_kind_phys * (ccn_gf(i)**0.64_kind_phys)
         
         if (aod_gf(i) < 0.007_kind_phys) then
            aod_gf(i) = 0.007_kind_phys
            ccn_gf(i) = (aod_gf(i)/0.0027_kind_phys)**(1.0_kind_phys/0.640_kind_phys)
         elseif (aod_gf(i) > aodc0) then
            aod_gf(i) = aodc0
            ccn_gf(i) = (aod_gf(i)/0.0027_kind_phys)**(1.0_kind_phys/0.640_kind_phys)
         endif
         
      enddo
!$acc end kernels
      !
            
! Scale dry mixing ratios for water wapor and cloud water to specific humidy / moist mixing ratios
        do i=its,itf
           do k=kts,ktf
            new_qv_spechum(i,k) = qv(i,k)/(1.0_kind_phys+qv(i,k))
            cnvw_moist(i,k) = cnvw(i,k)/(1.0_kind_phys+qv(i,k))
            ten_q(i,k,ntqv) = (new_qv_spechum(i,k) - qv_spechum(i,k))/dt
          end do
        end do
!
! Diagnostic tendency updates
!
        if(ldiag3d) then
          if(ishallow_g3.eq.1 .and. .not.flag_for_scnv_generic_tend) then
            uidx=dtidx(index_of_x_wind,index_of_process_scnv)
            vidx=dtidx(index_of_y_wind,index_of_process_scnv)
            tidx=dtidx(index_of_temperature,index_of_process_scnv)
            qidx=dtidx(100+ntqv,index_of_process_scnv)
            if(uidx>=1) then
!$acc kernels
              do k=kts,ktf
                dtend(:,k,uidx) = dtend(:,k,uidx) + cutens(:)*outus(:,k) * dt
              enddo
!$acc end kernels
            endif
            if(vidx>=1) then
!$acc kernels
              do k=kts,ktf
                dtend(:,k,vidx) = dtend(:,k,vidx) + cutens(:)*outvs(:,k) * dt
              enddo
!$acc end kernels
            endif
            if(tidx>=1) then
!$acc kernels
              do k=kts,ktf
                dtend(:,k,tidx) = dtend(:,k,tidx) + cutens(:)*outts(:,k) * dt
              enddo
!$acc end kernels
            endif
            if(qidx>=1) then
!$acc kernels
              do k=kts,ktf
                do i=its,itf
                  tem = cutens(i)*outqs(i,k)* dt
                  tem = tem/(1.0_kind_phys+tem)
                  dtend(i,k,qidx) = dtend(i,k,qidx) + tem
                enddo
              enddo
!$acc end kernels
            endif
          endif
          if((ideep.eq.1. .or. imid_gf.eq.1) .and. .not.flag_for_dcnv_generic_tend) then
            uidx=dtidx(index_of_x_wind,index_of_process_dcnv)
            vidx=dtidx(index_of_y_wind,index_of_process_dcnv)
            tidx=dtidx(index_of_temperature,index_of_process_dcnv)
            if(uidx>=1) then
!$acc kernels
              do k=kts,ktf
                dtend(:,k,uidx) = dtend(:,k,uidx) + (cuten*outu(:,k)+cutenm*outum(:,k)) * dt
              enddo
!$acc end kernels
            endif
            if(vidx>=1) then
!$acc kernels
              do k=kts,ktf
                dtend(:,k,vidx) = dtend(:,k,vidx) + (cuten*outv(:,k)+cutenm*outvm(:,k)) * dt
              enddo
!$acc end kernels
            endif
            if(tidx>=1) then
!$acc kernels
              do k=kts,ktf
                dtend(:,k,tidx) = dtend(:,k,tidx) + (cuten*outt(:,k)+cutenm*outtm(:,k)) * dt
              enddo
!$acc end kernels
            endif

            qidx=dtidx(100+ntqv,index_of_process_dcnv)
            if(qidx>=1) then
!$acc kernels
              do k=kts,ktf
                do i=its,itf
                  tem = (cuten(i)*outq(i,k) + cutenm(i)*outqm(i,k))* dt
                  tem = tem/(1.0_kind_phys+tem)
                  dtend(i,k,qidx) = dtend(i,k,qidx) + tem
                enddo
              enddo
!$acc end kernels
            endif
          endif
          if(allocated(clcw_save)) then
!$acc parallel loop collapse(2) private(tem_shal,tem_deep,tem,tem1,weight_sum,cliw_both,clcw_both)
            do k=kts,ktf
              do i=its,itf
                tem_shal = dt*(outqcs(i,k)*cutens(i)+outqcm(i,k)*cutenm(i))
                tem_deep = dt*(outqc(i,k)*cuten(i)+clw_ten(i,k))
                tem  = tem_shal+tem_deep
                tem1 = max(0.0, min(1.0, (tcr-t(i,k))*tcrf))
                weight_sum = abs(tem_shal)+abs(tem_deep)
                if(weight_sum<1e-12) then
                  cycle
                endif

                if (clcw_save(i,k) .gt. -999.0) then
                  cliw_both = max(0.,cliw_save(i,k) + tem * tem1) - cliw_save(i,k)
                  clcw_both = max(0.,clcw_save(i,k) + tem) - clcw_save(i,k)
                else if(cliw_idx>=1) then
                  cliw_both = max(0.,cliw_save(i,k) + tem) - cliw_save(i,k)
                  clcw_both = 0
                endif
                if(cliw_deep_idx>=1) then
                  dtend(i,k,cliw_deep_idx) = dtend(i,k,cliw_deep_idx) + abs(tem_deep)/weight_sum*cliw_both
                endif
                if(clcw_deep_idx>=1) then
                  dtend(i,k,clcw_deep_idx) = dtend(i,k,clcw_deep_idx) + abs(tem_deep)/weight_sum*clcw_both
                endif
                if(cliw_shal_idx>=1) then
                  dtend(i,k,cliw_shal_idx) = dtend(i,k,cliw_shal_idx) + abs(tem_shal)/weight_sum*cliw_both
                endif
                if(clcw_shal_idx>=1) then
                  dtend(i,k,clcw_shal_idx) = dtend(i,k,clcw_shal_idx) + abs(tem_shal)/weight_sum*clcw_both
                endif
              enddo
            enddo
!$acc end parallel
          endif
        endif
   end subroutine cu_c3_driver_ccpp_run
!>@}
end module cu_c3_driver_ccpp
