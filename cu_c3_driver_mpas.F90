 subroutine cu_c3_driver_mpas(    &
               dt                    &
              ,confrq                &
              ,dx_p                  &
              ,areaCell              &
              ,lats                  &
              ,lons                  &
              ,u                     &
              ,v                     &
              ,w                     &
              ,temp                  &
              ,rvap                  &
              ,rho                   &
              ,press                 &
              ,pi                    &
              ,p8w                   &
              ,dz8w                  &
              ,topt                  &
              ,xland                 &
              ,sflux_t               &
              ,sflux_r               &
              ,temp2m                &
              ,wlpool                &
              ,mpas_cape             & ! check if it is updated before entering GF
              ,mpas_cin              & ! check if it is updated before entering GF
              ,kpbl                  &
              ,tke_pbl               &
              ,turb_len_scale        &
              ,buoyx                 &
              ,cnvcf                 &
              ,rthblten              & ! tendency of potential temperature due to pbl processes
              ,rqvblten              & ! tendency of water vapor mixing ratio due to pbl processes
              ,rthratenlw            & ! tendency of potential temperature due to long-wave radiation
              ,rthratensw            & ! tendency of potential temperature due to short-wave radiation
              ,rthdyten              & ! tendency of potential temperature due to dynamics plus filters
              ,rqvdyten              & ! tendency of water vapor mixing ration due to dynamics plus filters
              !---- output ----      
              ,raincv                &
              ,conprr                &
              ,lightn_dens           &
              ,sigma_deep            &
              ,rthcuten              &
              ,rqvcuten              &
              ,rqccuten              &
              ,rqicuten              &
              ,rucuten               &
              ,rvcuten               &
              ,rbuoyxcuten           &
              ,rcnvcfcuten           &

              ,sub3d_rthcuten        &
              ,sub3d_rqvcuten        &
              ,sub3d_rucuten         &
              ,sub3d_rvcuten         &

              ,rupmfxcu              &
              ,rdnmfxcu              &
              !
              ,rmfxdpcu              & 
              ,rmfxdncu              & 
              ,rmfxmdcu              &
              ,rmfxshcu              &
              ,rtopdpcu              &
              ,rtopmdcu              &
              ,rtopshcu              &
              ,rbotdpcu              &
              ,var2d1                &
              ,var2d2                &
              ,var3d1                &
              !
              ,ids, ide, jds, jde, kds, kde   &
              ,ims, ime, jms, jme, kms, kme   &
              ,its, ite, jts, jte, kts, kte   &
              ,itimestep,mynum                &
              ,dp_dens, sh_dens, cg_dens      &
              ,dp_ierr, sh_ierr, cg_ierr      )


! --- vars need to be implemented
                           ! ,TRACER                &
                           ! ,rnlcuten              &
                           ! ,rnicuten              &
                           ! ,rchemcuten            &


!-----------------------------



      use modGate, only: cupout, rundata, p_nvar_grads, jl, p_use_gate, ppres, ptemp, pq, pu, &
     &                  runlabel, runname, pv, pvervel, pgeo, zqr, zadvq, zadvt
      use modConstants, only: c_rgas, c_cp, c_alvl, c_grav, c_T00, c_Tice, p_mintracer, p_smaller_qv
      use module_cu_c3, only: maxiens, deep, shal, mid, max_tq_tend, nmp,               &
     &                        icumulus_gf, cumulus_type, closure_choice, cum_entr_rate,  &
     &                        cum_use_excess, cum_use_smooth_tend, cum_min_cloud_depth,  &
     &                        cum_hei_down_land, cum_hei_down_ocean,                      &
     &                        cum_hei_updf_land, cum_hei_updf_ocean,                      &
     &                        cum_max_edt_land, cum_max_edt_ocean,                        &
     &                        cum_ave_layer, cum_t_star, cum_fr_min_entr,                &
     &                        c0_deep, c0_shal, c0_mid, autoconv, downdraft,             &
     &                        use_momentum_transp, use_tracer_transp, convection_tracer, &
     &                        use_sub3d, use_pass_cloudvol, liq_ice_number_conc,         &
     &                        Hcts, CUP_C3, calc_lcl, set_Tq_pertub, allev_initial,      &
     &                        FractLiqF

      implicit none
      !------------------------------------------------------------------------
      !intent in arguments:
      integer                                 ,intent(in):: ids,ide,jds,jde,kds,kde & 
                                                           ,ims,ime,jms,jme,kms,kme & 
                                                           ,its,ite,jts,jte,kts,kte &
                                                           ,itimestep, mynum
      integer, dimension(ims:ime,jms:jme)     ,intent(in):: kpbl

      real                                    ,intent(in):: dt, confrq
      real, dimension(ims:ime,jms:jme)        ,intent(in):: areaCell,dx_p,lats,lons
      real, dimension(ims:ime,jms:jme)        ,intent(in):: sflux_r,sflux_t,topt,xland,temp2m
      real, dimension(ims:ime,jms:jme)        ,intent(in):: mpas_cape,mpas_cin

      real, dimension(ims:ime,kms:kme,jms:jme),intent(in):: u,v,w,press,pi,rvap,rho,temp,tke_pbl  
      real, dimension(ims:ime,kms:kme,jms:jme),intent(in):: dz8w,p8w,turb_len_scale
      real, dimension(ims:ime,kms:kme,jms:jme),intent(in):: buoyx, cnvcf

      real, dimension(ims:ime,kms:kme,jms:jme),intent(in):: rqvblten,rthblten,rthratenlw,rthratensw &
                                                           ,rthdyten,rqvdyten
      !-- intent in,out arguments
      real,dimension(ims:ime,jms:jme)         ,intent(inout):: raincv,conprr,wlpool,lightn_dens,sigma_deep
      real,dimension(ims:ime,jms:jme)         ,intent(inout):: dp_dens,sh_dens,cg_dens
      real,dimension(ims:ime,jms:jme)         ,intent(inout):: dp_ierr,sh_ierr,cg_ierr

      real,dimension(ims:ime,kms:kme,jms:jme) ,intent(inout):: rthcuten,rqvcuten,rqccuten,rqicuten &
                                                              ,rucuten,rvcuten,rbuoyxcuten,rcnvcfcuten

      real,dimension(ims:ime,kms:kme,jms:jme) ,intent(inout):: sub3d_rthcuten         &
                                                               ,sub3d_rqvcuten        &
                                                               ,sub3d_rucuten         &
                                                               ,sub3d_rvcuten
      !---intent out arguments
      real,dimension(ims:ime,jms:jme)         ,intent(out):: rmfxdpcu,rmfxmdcu,rmfxshcu &
                                                            ,rtopdpcu,rtopmdcu,rtopshcu &
                                                            ,rbotdpcu,rmfxdncu
      real,dimension(ims:ime,jms:jme)         ,intent(out):: var2d1,var2d2

      real,dimension(ims:ime,kms:kme,jms:jme) ,intent(out):: rupmfxcu,rdnmfxcu
      real,dimension(ims:ime,kms:kme,jms:jme) ,intent(out):: var3d1

      !----------------------- local and future in/out arguments
      integer, parameter :: mtp = 1
      integer, parameter :: ON = 1, OFF = 0
      integer :: use_excess, use_smooth_tend, use_cold_start
      real    :: hei_down_land, hei_down_ocean, hei_updf_land, hei_updf_ocean
      real    :: max_edt_land, max_edt_ocean, fr_min_entr, t_star, ave_layer
      real    :: c0, min_cloud_depth, fac_cold_start
      real :: time= 0., itime1 = 0.
      real:: FSCAV(mtp)

      real,    dimension(its:ite,jts:jte) :: aot500 ,sfc_press,col_sat &
                                                    ,stochastic_sig,dx2d

      real,    dimension(nmp,kts:kte,its:ite,jts:jte) :: mp_ice,mp_liq,mp_cf

      real,    dimension(nmp,kts:kte,its:ite,jts:jte) :: sub_mpqi,sub_mpql,sub_mpcf

      !-***** TRACER has different data structure   (i,j,k,ispc) *********
      real,    dimension(its:ite,jts:jte,kts:kte,mtp)  :: TRACER
      !-***** rchemcuten uses the GF data structure (ispc,k,i,j) *********
      real,    dimension(mtp,kts:kte,its:ite,jts:jte)  :: rchemcuten

      !- for convective transport and cloud/radiation (OUT) 
      integer, dimension(its:ite,jts:jte) :: do_this_column

      integer, dimension(its:ite,jts:jte,maxiens) ::    &
          ierr4d                    &
         ,jmin4d                    &
         ,klcl4d                    &
         ,k224d                     &
         ,kbcon4d                   &
         ,ktop4d                    &
         ,kstabi4d                  &
         ,kstabm4d

      real,dimension(its:ite,jts:jte,maxiens)  ::    &
          cprr4d                    &
         ,xmb4d                     &
         ,edt4d                     &
         ,pwav4d                    &
         ,sigma4d
      
      real,dimension(kts:kte,its:ite,jts:jte,maxiens) ::    &
          pcup5d                    &
         ,up_massentr5d             &
         ,up_massdetr5d             &
         ,dd_massentr5d             &
         ,dd_massdetr5d             &
         ,zup5d                     &
         ,zdn5d                     &
         ,prup5d                    &
         ,prdn5d                    &
         ,clwup5d                   &
         ,tup5d                     &
         ,conv_cld_fr5d
      !----------------------------------------------------------------------
      ! LOCAL VARS
      ! basic environmental input includes

      real,   dimension (kts:kte,its:ite) ::  zo,temp_old,qv_old,po,us,vs,rhoi,phil    &
                                             ,temp_new,qv_new,temp_new_dp,qv_new_dp    &
                                             ,temp_new_bl,qv_new_bl,temp_tendqv        &
                                             ,temp_new_adv,qv_new_adv,dhdt             &
                                             ,cnvcf2d,turb_len_scale2d,buoyx2d         &
                                             ,qv_curr,piexner            
      
      real,   dimension (kts:kte,its:ite,maxiens) ::  outt,outq,outqc,outu,outv,outbuoy &
                                                     ,outnliq,outnice
      real,   dimension (kts:kte,its:ite,maxiens) ::  subten_Q,subten_T,subten_U,subten_v

      real,   dimension (mtp,kts:kte,its:ite)         :: se_chem
      real,   dimension (mtp,kts:kte,its:ite,maxiens) :: out_chem

      real,   dimension (nmp,kts:kte,its:ite)         :: mpqi,mpql,mpcf
      real,   dimension (nmp,kts:kte,its:ite,maxiens) :: outmpqi,outmpql,outmpcf

      real,   dimension (its:ite)   :: ter11, xlandi,pbl,zws,ccn,psur             &
                                      ,ztexec,zqexec,h_sfc_flux,le_sfc_flux,tsur  &
                                      ,xlons,xlats,fixout_qv,cum_ztexec,cum_zqexec&
                                      ,zlcl_sfc,plcl_sfc,tlcl_sfc,tke_mean

      real,   dimension (kts:kte,its:ite,1:ens4)      ::  omeg
      real,   dimension (its:ite,1:ens4)              ::  mconv

      real,   dimension (kts:kte) :: min_tend,distance
      integer,dimension (its:ite) :: kpbli,max_ktop

      integer :: i,j,k,kr,n,itf,jtf,ktf,ispc,zmax,status

      real :: dp,dq, dtdt,pten,pqen,paph,zrho,pahfs,pqhfl,zkhvfl,pgeoh
      real :: fixouts,dt_inv,temp2theta,theta2temp,total_dp
      integer :: jlx,kk,plume,ii_plume, l_unit,irec,rec_size,iloc,kloc
      character(len=6)  :: c_itimestep
      character(len=128) :: outname
      logical :: debug_mpas 
      !print*,'================================================================'
      !print*,'itimestep',itimestep
      !print*,'mpas qv', maxval(abs(rqvcuten(its:ite,kts:kte,jts:jte))),mynum
      !print*,'mpas th', maxval(abs(rthcuten(its:ite,kts:kte,jts:jte))),mynum
      !print*,'mpas qi', maxval(abs(rqicuten(its:ite,kts:kte,jts:jte))),mynum
      !print*,'mpas qc', maxval(abs(rqccuten(its:ite,kts:kte,jts:jte))),mynum
      !print*,'mpas bu', maxval(abs(rbuoyxcuten(its:ite,kts:kte,jts:jte))),mynum
      !print*,'mpas u ', maxval(abs(rucuten(its:ite,kts:kte,jts:jte))),mynum
      !print*,'mpas v ', maxval(abs(rvcuten(its:ite,kts:kte,jts:jte))),mynum
      !print*,'================================================================'

      !--- these arrays must be reset every timestep.
         ierr4d         = 0      
         jmin4d         = kts
         klcl4d         = kts
         k224d          = kts
         kbcon4d        = kts
         ktop4d         = kts
         kstabi4d       = kts
         kstabm4d       = kts
      !if(use_tracer_transp == 1) then
         cprr4d         = 0.0
         xmb4d          = 0.0
         edt4d          = 0.0
         pwav4d         = 0.0
         sigma4d        = 0.0
         pcup5d         = 0.0
         up_massentr5d  = 0.0
         up_massdetr5d  = 0.0
         dd_massentr5d  = 0.0
         dd_massdetr5d  = 0.0
         zup5d          = 0.0
         zdn5d          = 0.0
         prup5d         = 0.0
         prdn5d         = 0.0
         clwup5d        = 0.0
         tup5d          = 0.0
         conv_cld_fr5d  = 0.0
         TRACER         = 0.0 
      !endif
      !print*,"============================================ in convection"

      !----------------------------------------------------------------------
      !-do not change this
      !itf=ite
      !ktf=kte-1
      !jtf=jte
      itf = min(ite,ide-1)
      ktf = min(kte,kde-1)
      jtf = min(jte,jde-1)
      
      int_time   = int_time + dt
      WHOAMI_ALL = mynum
      time_in    = time
      itime1_in  = itime1

      FSCAV         (:)   = 0.1
      stochastic_sig(:,:) = 1.0  
      aot500        (:,:) = 0.0
      col_sat       (:,:) = 0.0  !falta
      !----------------------------------------------------------------------
      !print*, 'domain', its,itf,ite,kts,ktf,kte,jts,jtf,jte
      !print*,'ins---',maxval(conprr(its:ite,jts:jte))*3600 &
      !         ,maxval(abs(rthcuten(its:ite,kts:kte,jts:jte)))*86400.&
      !         ,maxval(abs(rqvcuten(its:ite,kts:kte,jts:jte)))*86400.*2.5e+6/1004.
      !print*,"CAPE CIN",maxval(mpas_cape),minval(mpas_cape),maxval(mpas_cin),minval(mpas_cin)
      
      !print*,"buoyx kJ/kg - Wlpool m/s",1.e-3*maxval(buoyx(its:ite,kts:kte,jts:jte)),1.e-3*minval(buoyx(its:ite,kts:kte,jts:jte))&
      !              ,maxval(wlpool(its:ite,jts:jte)),minval(wlpool(its:ite,jts:jte))
      !fac_cold_start = 1.0
       if(use_cold_start == 1) fac_cold_start = allev_initial(itimestep,dt)

      !-- big loop over j dimension
      j_loop: do j = jts,jtf
         JCOL = J
         conprr     (:,j) = 0.0
         raincv     (:,j) = 0.0
         lightn_dens(:,j) = 0.0
         sigma_deep (:,j) = 0.0
         dp_dens    (:,j) = 0.0
         sh_dens    (:,j) = 0.0
         cg_dens    (:,j) = 0.0
         var2d1     (:,j) = 0.0
         var2d2     (:,j) = 0.0

         !-- initialization
         ztexec   (:) = 0.0
         zqexec   (:) = 0.0
         fixout_qv(:) = 1.0
         !
         !--- (k,i)
         temp_tendqv (:,:) = 0.0
         !
         omeg   (:,:,:) = 0.0
         !- tendencies (w/ maxiens)
         outt   (:,:,:) = 0.0
         outu   (:,:,:) = 0.0
         outv   (:,:,:) = 0.0
         outq   (:,:,:) = 0.0
         outqc  (:,:,:) = 0.0
         outnice(:,:,:) = 0.0  !falta
         outnliq(:,:,:) = 0.0  !falta
         outbuoy(:,:,:) = 0.0  
         if(use_sub3d > 0) then 
            subten_Q(:,:,:) = 0.0
            subten_T(:,:,:) = 0.0
            subten_U(:,:,:) = 0.0
            subten_V(:,:,:) = 0.0
         endif 
         if(APPLY_SUB_MP == 1) then
           !- tendencies (w/ nmp and maxiens)
            outmpqi(:,:,:,:) = 0.0
            outmpql(:,:,:,:) = 0.0
            outmpcf(:,:,:,:) = 0.0
         endif

         if(USE_TRACER_TRANSP == 1) then
            out_chem(:,:,:,:) = 0.0
         endif
         !
         if(autoconv == 2) then
            do i= its,itf
               ccn(i) = max( 100., ( 370.37*(0.01+MAX(0.,aot500(i,j))))**1.555 )
            enddo
         else
            do i= its,itf
               ccn(i) = 100.
            enddo
         endif

         do i=its,itf
            dx2d  (i,j) = dx_p(i,j) ! grid spacing 
            
            !xlandi(i) = xland(i,j)!flag < 1 para land
            !                      !flag  =1 para water
            !<var name="xland" type="real" dimensions="nCells Time" units="unitless"
            ! description="land-ocean mask (1=land including sea-ice ; 2=ocean)"/>
            if(xland(i,j) <  1.5 ) xlandi(i) = 0.
            if(xland(i,j) >= 1.5 ) xlandi(i) = 1.

            psur  (i) = p8w(i,1,j)*1.e-2 ! mbar
            tsur  (i) = temp2m(i,j) 
            
            ter11 (i) = max(0.,topt(i,j))
            kpbli (i) = kpbl(i,j)
            xlons (i) = lons(i,j) ! in degrees
            xlats (i) = lats(i,j) ! in degress
          enddo
          !- heigths
          do i=its,itf
            zo(kts,i) = ter11(i)  + 0.5*dz8w(i,1,j)
            do k = kts+1, ktf
              zo(k,i) = zo(k-1,i) + 0.5*(dz8w(i,k-1,j)+dz8w(i,k,j))
            enddo
          enddo
          do i=its,itf
            do k=kts,ktf
               !- current pressure, temp and water vapor mix ratio
               po      (k,i)  = press(i,k,j)*1.e-2 !mbar
               temp_old(k,i)  = temp(i,k,j) ! K
               qv_old  (k,i)  = rvap(i,k,j) ! kg/kg @ begin of the timestep
               qv_curr (k,i)  = rvap(i,k,j) ! kg/kg check this ! current (after dynamics + physical processes called before GF)
               piexner (k,i)  = pi(i,k,j)   ! Exner function
               !- air density
               rhoi    (k,i)  = rho(i,k,j)  
               !- horiz wind velocities
               us      (k,i)  =  u (i,k,j)
               vs      (k,i)  =  v (i,k,j)
               !- cloud fraction 
               cnvcf2d (k,i)  = cnvcf(i,k,j) 
               !-buoyancy excess
               buoyx2d(k,i)  = buoyx(i,k,j)
               !-- turb length scale
               turb_len_scale2d (k,i) = turb_len_scale (i,k,j)
            end do
          end do
          do i=its,itf
            !-- max convective cloud top is in 50 mbar
            max_ktop(i) =  minloc(abs(po(:,i)-50.),1)
          enddo
          do i=its,itf
            do k=kts,ktf
                  !-- convert from theta to temperature 
                  theta2temp = pi(i,k,j)
                  
                  !-- temp/water vapor projected by the ´large-scale' forcing
                  temp_new(k,i) = temp_old(k,i) + (rthblten  (i,k,j) + rthdyten  (i,k,j) + &
                                                   rthratenlw(i,k,j) + rthratensw(i,k,j) ) * dt * theta2temp
                  qv_new  (k,i) =   qv_old(k,i) + (rqvblten(i,k,j)+rqvdyten(i,k,j)) * dt
                  
                  !- temp/water vapor modified only by bl processes
                  temp_new_BL (k,i)= temp_old(k,i)  +  rthblten(i,k,j) * dt * theta2temp
                  qv_new_BL   (k,i)=   qv_old(k,i)  +  rqvblten(i,k,j) * dt
                  
                  !- temp/water vapor modified only by advection + filters
                  temp_new_ADV(k,i) = temp_old(k,i)  +  rthdyten(i,k,j) * dt * theta2temp
                  qv_new_ADV  (k,i) =   qv_old(k,i)  +  rqvdyten(i,k,j) * dt
                  
                  qv_new      (k,i) = max(p_smaller_qv,qv_new    (k,i))
                  qv_new_BL   (k,i) = max(p_smaller_qv,qv_new_BL (k,i))
                  qv_new_ADV  (k,i) = max(p_smaller_qv,qv_new_ADV(k,i))

                  !- only pbl forcing changes moist static energy
                  dhdt(k,i) = c_cp * rthblten(i,k,j)*theta2temp + c_alvl * rqvblten(i,k,j)

                  !- all forcings change moist static energy
                  dhdt(k,i) = dhdt(k,i) +   c_cp*(rthdyten(i,k,j)+rthratenlw(i,k,j)+rthratensw(i,k,j))*theta2temp &
                                        + c_alvl*(rqvdyten(i,k,j))
            end do
          end do

          !-- calculation of omega vertical velocity and the moisture convergence:
          do k = kts+1, ktf
            do i = its, itf
               omeg(k,i,:) = -c_grav*0.5*(rho(i,k,j)+rho(i,k-1,j))*w(i,k,j)
            end do
          end do
          do i = its, itf
            do n=1,ens4
                mconv(i,n) = 0.
                do k = kts+1, ktf         
                  dq         = qv_old(k,i)-qv_old(k-1,i)
                  mconv(i,n) = mconv(i,n) + omeg(k,i,n)*dq/c_grav
               end do
               mconv(i,n) = max(0., mconv(i,n))
            end do
          end do

          if(APPLY_SUB_MP == 1) then ! - check arrays later 
            do i=its,itf
               do k=kts,ktf
                  kr=k   !+1   !<<<< only kr=k
                  !- microphysics ice and liq mixing ratio, and cloud fraction of the host model
                  !- (only subsidence is applied)
                  mpqi   (:,k,i) = mp_ice  (:,kr,i,j) ! kg/kg
                  mpql   (:,k,i) = mp_liq  (:,kr,i,j) ! kg/kg
                  mpcf   (:,k,i) = mp_cf   (:,kr,i,j) ! 1
               enddo
            enddo
          endif
          if(USE_TRACER_TRANSP == 1) then  ! - check arrays later 
            do i=its,itf
               do k=kts,kte
                  kr=k !+1
                  !- atmos composition
                  do ispc=1,mtp
                     se_chem(ispc,k,i) = max(p_mintracer, TRACER(i,j,kr,ispc))
                  enddo
               enddo
            enddo
          endif
          !- pbl  (i) = depth of pbl layer (m)
          !- kpbli(i) = PBL index of zo(k,i)
          do i=its,itf
            pbl  (i)  = zo(kpbli(i),i) - topt(i,j)
          enddo
          do i = its, itf
            tke_mean(i) = 0.0
            total_dp = 0.0
            do k = kts+1, kpbli(i)-1
               dp = -0.5*(po(k+1,i)-po(k-1,i))
               tke_mean(i) = tke_mean(i) + tke_pbl(i,k,j)*dp
               total_dp = total_dp + dp
            end do
            tke_mean(i) = max(1.e-5, tke_mean(i)/(total_dp+1.e-6))
          end do

         !- begin: for GATE soundings-------------------------------------------
         if(P_USE_GATE) then
            if(CLEV_GRID == 0) stop "use_gate requires CLEV_GRID 1 or 2"
            if(USE_TRACER_TRANSP == 1) then
               ispc_CO=1
               if( .not. allocated(Hcts)) allocate(Hcts(mtp))
               CHEM_NAME_MASK (:) = 1
               !--- dummy initization FSCAV
               do i=1,mtp
                  !FSCAV(i) = 0.1  !km^-1

                  FSCAV(i) = 1.e-5  !km^-1
                  Hcts(i)%hstar  = 0.0 !8.300e+4! 2.4E+3 !59.
                  Hcts(i)%dhr    = 0.0 !7400.   !5000.  !4200.
                  Hcts(i)%ak0    = 0.0
                  Hcts(i)%dak    = 0.0
                   ! H2O2      0.00000      8.300e+4    7400.00000       0.00000       0.00000
                   ! HNO3      0.00000      2.100e+5    8700.00000       0.00000       0.00000
                   ! NH3       0.00000      59.00000    4200.00000       0.00000       0.00000
                   ! SO2       0.00000      2.400e+3    5000.00000       0.00000       0.00000
               enddo
               do i=its,itf
                  se_chem(1:mtp,kts:kpbli(i)-1,i) = 1.+1.e-6
                  do k=kpbli(i),kte
                     se_chem(1:mtp,k,i) = 1.*exp(-(max(0.,0.9*float(k-kpbli(i)))/float(kpbli(i))))+1.e-6
                  enddo
                  do k=kts+1,kte-1
                     se_chem(1:mtp,k,i) = 1./3. *( se_chem(1:mtp,k,i) + se_chem(1:mtp,k-1,i) + se_chem(1:mtp,k+1,i))
                  enddo
               enddo
            endif

            !--- only for GATE soundingg
            if(trim(RUNDATA) == "GATE.dat") then
               jlx= jl
              !jlx= 10 ! to run with only one soundings
 
               do i=its,itf
                  do k=kts,kte
                     po       (k,i) = 0.5*(ppres(jlx,k)+ppres(jlx,min(kte,k+1)))
                     temp_old (k,i) = ptemp(jlx,k)+273.15
                     qv_old   (k,i) = pq(jlx,k)/1000.
                     us       (k,i) = pu(jlx,k)
                     vs       (k,i) = pv(jlx,k)
                     omeg     (k,i,:)=pvervel(jlx,k)
                     phil     (k,i) = pgeo(jlx,k)*c_grav   !geo
                     rhoi     (k,i) = 1.e2*po(k,i)/(c_rgas*temp_old(k,i))
                  enddo

                  do k=kts,kte
                     mpql     (:,k,i) = 0.
                     mpql     (:,k,i) = 0.
                     mpcf     (:,k,i) = 0.
                     if(po(k,i) > 900. .or. po(k,i)<300.) cycle
                     pqen  =  exp((-3.e-5*(po(k,i)-550.)**2))
                     pten  =  min(1., (max(0.,(temp_old(k,i)-c_Tice))/(c_T00-c_Tice))**2)
                     mpql  (:,k,i) =3.*pqen* pten
                     mpqi  (:,k,i) =3.*pqen*(1.- pten)
                     mpcf  (:,k,i) = (mpqi  (:,k,i)+mpql  (:,k,i))*100.
                  enddo

                  do k=kts,kte
                     zo       (k,i) = 0.5*(phil(k,i)+phil(min(kte,k+1),i))/c_grav    !meters
                  enddo
                  ter11(i)  = phil(1,i)/c_grav  ! phil is given in g*h.
                  psur (i)  = ppres(jlx,1)
                  tsur (i)  = temp2m(i,j) !temp_old(i,1)
                  kpbli(i)  = 5
                  pbl  (i)  = zo(kpbli(i),i)
                  zws  (i)  = 1.0 ! wstar
                  do k=kts,ktf
                     temp_new(k,i) = temp_old(k,i) + dt *(zadvt(jlx,k)+zqr(jlx,k))/86400.
                     qv_new  (k,i) = qv_old  (k,i) + dt * zadvq(jlx,k)

                     temp_new_dp (k,i) = temp_old(k,i) + dt *(zadvt(jlx,k)+zqr(jlx,k))/86400.
                     qv_new_dp   (k,i) = qv_old  (k,i) + dt * zadvq(jlx,k)

                     temp_new_bl (k,i) = temp_new_dp(k,i)
                     qv_new_bl   (k,i) = qv_new_dp  (k,i)
                     temp_new_adv(k,i) = temp_old   (k,i) + dt * zadvt(jlx,k)/86400.
                     qv_new_adv  (k,i) = qv_old     (k,i) + dt * zadvq(jlx,k)
                     dhdt        (k,i)= c_cp*(temp_new_dp(k,i)-temp_old(k,i))+ &
                                        c_alvl*(qv_new_dp(k,i)-qv_old(k,i))
                  enddo
               enddo
            endif
         endif 
         !- end for GATE soundings-------------------------------------------
         !
         !- get excess T and Q for source air parcels
         do i=its,itf
            pten = temp_old(1,i)
            pqen = qv_old  (1,i)
            paph = 100.*psur(i)
            zrho = paph/(287.04*(temp_old(1,i)*(1.+0.608*qv_old(1,i))))
            !- sensible and latent sfc fluxes for the heat-engine closure
            !h_sfc_flux (i)=zrho*c_cp  *sflux_t(i,j)!W/m^2  in MPAS h is already in W/m2
            h_sfc_flux (i)= sflux_t(i,j)            !W/m^2
            le_sfc_flux(i)= zrho*c_alvl*sflux_r(i,j)!W/m^2
            !
            !- local le and h fluxes for calculate W*
            pahfs=-sflux_t(i,j)  !W/m^2 : in MPAS h is already in W/m2
            pqhfl=-sflux_r(i,j)  !kg/m^2/s
            !- buoyancy flux (h+le)
            zkhvfl= (pahfs/1004.64+0.608*pten*pqhfl)/zrho ! K m s-1
            !- depth of 1st model layer
            !- (zo(1)-top is ~ 1/2 of the depth of 1st model layer, => mult by 2)
            pgeoh =  2.*( zo(1,i)-topt(i,j) )*c_grav ! m+2 s-2
            !-convective-scale velocity w*
            !- in the future, change 0.001 by ustar^3
            zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten) ! m+3 s-3
            !
            !-- get LCL properties for parcels from surface level
            call calc_lcl(pten,paph,pqen,tlcl_sfc(i),plcl_sfc(i),zlcl_sfc(i))
            zlcl_sfc(i) = max(zlcl_sfc(i), 0.)
            !print*,'lcl',minval(zlcl_sfc(its:itf)),maxval(zlcl_sfc(its:itf))
            !
            if(zws(i) > tiny(pgeoh)) then
               !-convective-scale velocity w*
               zws(i) = 1.2*zws(i)**.3333
               !- temperature excess
               ztexec(i)     = max(0.,-1.5*pahfs/(zrho*zws(i)*1004.64)) ! K
               !- moisture  excess
               zqexec(i)     = max(0.,-1.5*pqhfl/(zrho*zws(i)))        !kg kg-1
            endif   ! zws > 0
            !
            !- zws for shallow convection closure (Grant 2001)
            !- depth of the pbl
            pgeoh = pbl(i)*c_grav
            !-convective-scale velocity W* (m/s)
            zws(i) = max(0.,0.001-1.5*0.41*zkhvfl*pgeoh/pten)
            zws(i) = 1.2*zws(i)**.3333
         enddo
         !
         !------ CALL CUMULUS PARAMETERIZATION
         !
         do ii_plume = 1, maxiens

            if(ii_plume == 1) then
               plume = deep
               c0 = c0_deep
            endif
            if(ii_plume == 2) then
               plume = shal
               c0 = c0_shal
               if(icumulus_gf(shal) == 2 ) where(ierr4d(:,j,deep) == 0) ierr4d(:,j,shal) = -99
              !if(icumulus_gf(shal) == 2 ) c0 = 0.0
            endif
            if(ii_plume == 3) then
               plume = mid
               c0 = c0_mid
               if(icumulus_gf(mid) == 2 ) where(ierr4d(:,j,deep) == 0) ierr4d(:,j,mid) = -99
            endif
            
            if(icumulus_gf(plume) == OFF ) cycle

            use_smooth_tend=  cum_use_smooth_tend(plume)
            min_cloud_depth=  cum_min_cloud_depth(plume)
            hei_down_land  =  cum_hei_down_land  (plume)
            hei_down_ocean =  cum_hei_down_ocean (plume)
            hei_updf_land  =  cum_hei_updf_land  (plume)
            hei_updf_ocean =  cum_hei_updf_ocean (plume)
            max_edt_land   =  cum_max_edt_land   (plume)
            max_edt_ocean  =  cum_max_edt_ocean  (plume)
            use_excess     =  cum_use_excess     (plume)
            ave_layer      =  cum_ave_layer      (plume)
            T_star         =  cum_t_star         (plume)
            fr_min_entr    =  cum_fr_min_entr    (plume)

            !
            !-- set the temp and water vapor anomalies from the sub-grid scale variability 
            call set_Tq_pertub (use_excess,its,ite,itf,xlandi,ztexec,zqexec,cum_ztexec,cum_zqexec)
            !
            call CUP_C3(its,ite,kts,kte, itf,ktf, mtp, nmp, FSCAV  &
                        ,cumulus_type  (plume)            &
                        ,closure_choice(plume)            &
                        ,cum_entr_rate (plume)            &
                        ,cum_use_excess(plume)            &
                        ,use_smooth_tend                  &
                        ,min_cloud_depth                  &
                        ,hei_down_land                    &
                        ,hei_down_ocean                   &
                        ,hei_updf_land                    &
                        ,hei_updf_ocean                   &
                        ,max_edt_land                     &
                        ,max_edt_ocean                    &
                        ,c0                               &
                        ,use_excess                       &
                        ,ave_layer                        &
                        ,T_star                           &
                        ,fr_min_entr                      &
                        !- input data
                        ,dx2d          (:,j)              &
                        ,stochastic_sig(:,j)              &
                        ,col_sat       (:,j)              &
                        ,tke_mean      (:)                &
                        ,wlpool        (:,j)              &
                        ,max_ktop                         &
                        ,dt                               &
                        ,kpbli                            &
                        ,cum_ztexec                       &
                        ,cum_zqexec                       &
                        ,ccn                              &
                        ,rhoi                             &
                        ,omeg                             &
                        ,temp_old                         &
                        ,qv_old                           &
                        ,ter11                            &
                        ,h_sfc_flux                       &
                        ,le_sfc_flux                      &
                        ,zlcl_sfc                         &
                        ,xlons                            &
                        ,xlats                            &
                        ,xlandi                           &
                        ,temp_new                         &
                        ,qv_new                           &
                        ,temp_new_BL                      &
                        ,qv_new_BL                        &
                        ,temp_new_ADV                     &
                        ,qv_new_ADV                       &
                        ,zo                               &
                        ,po                               &
                        ,piexner                          &
                        ,tsur                             &
                        ,psur                             &
                        ,us                               &
                        ,vs                               &
                        ,rhoi                             &
                        ,se_chem                          &
                        ,zws                              &
                        ,dhdt                             &
                        ,buoyx2d                          &
                        ,cnvcf2d                          &
                        ,turb_len_scale2d                 &
                        ,mpqi                             &
                        ,mpql                             &
                        ,mpcf                             &
                        !output data
                        ,outt                 (:,:,plume) &
                        ,outq                 (:,:,plume) &
                        ,outqc                (:,:,plume) &
                        ,outu                 (:,:,plume) &
                        ,outv                 (:,:,plume) &
                        ,subten_Q             (:,:,plume) &
                        ,subten_T             (:,:,plume) &
                        ,subten_U             (:,:,plume) &
                        ,subten_V             (:,:,plume) &
                        ,outnliq              (:,:,plume) &
                        ,outnice              (:,:,plume) &
                        ,outbuoy              (:,:,plume) &
                        ,outmpqi            (:,:,:,plume) &
                        ,outmpql            (:,:,:,plume) &
                        ,outmpcf            (:,:,:,plume) &
                        ,out_chem           (:,:,:,plume) &
                        !- for convective transport
                        ,ierr4d               (:,j,plume) &
                        ,jmin4d               (:,j,plume) &
                        ,klcl4d               (:,j,plume) &
                        ,k224d                (:,j,plume) &
                        ,kbcon4d              (:,j,plume) &
                        ,ktop4d               (:,j,plume) &
                        ,kstabi4d             (:,j,plume) &
                        ,kstabm4d             (:,j,plume) &
                        ,cprr4d               (:,j,plume) &
                        ,xmb4d                (:,j,plume) &
                        ,edt4d                (:,j,plume) &
                        ,pwav4d               (:,j,plume) &
                        ,sigma4d              (:,j,plume) &
                        ,pcup5d             (:,:,j,plume) &
                        ,up_massentr5d      (:,:,j,plume) &
                        ,up_massdetr5d      (:,:,j,plume) &
                        ,dd_massentr5d      (:,:,j,plume) &
                        ,dd_massdetr5d      (:,:,j,plume) &
                        ,zup5d              (:,:,j,plume) &
                        ,zdn5d              (:,:,j,plume) &
                        ,prup5d             (:,:,j,plume) &
                        ,prdn5d             (:,:,j,plume) &
                        ,clwup5d            (:,:,j,plume) &
                        ,tup5d              (:,:,j,plume) &
                        ,conv_cld_fr5d      (:,:,j,plume) &
                        !-- for diag
                        ,lightn_dens  (:,j)               &
                        ,var2d1       (:,j)               &
                        ,var2d2       (:,j)               &
                        )
         enddo !- plume

         !--- reset ierr4d to value different of zero in case the correspondent
         !--- plume (shalllow, congestus, deep) was not actually used
         do n=1,maxiens
            if(icumulus_gf(n) == OFF ) ierr4d (:,j,n) = -99
         enddo

         do i=its,itf
            do_this_column(i,j) = 0
            loop1:  do n=1,maxiens
               if(ierr4d (i,j,n) == 0 ) then
                  do_this_column(i,j) = 1
                  !print*,'conv on',i, do_this_column(i,j)
                  exit loop1
               endif
            enddo loop1
         enddo
         !
         !-- output
         do i=its,itf
             dp_ierr(i,j) = ierr4d(i,j,deep)
             sh_ierr(i,j) = ierr4d(i,j,shal)
             cg_ierr(i,j) = ierr4d(i,j,mid)
             if(ierr4d (i,j,deep) == 0 ) then
               rmfxdpcu(i,j)   = xmb4d(i,j,deep)                 !-- updraft mass flux @ cloud base (kg/m2/s)
               rmfxdncu(i,j)   = xmb4d(i,j,deep)*edt4d(i,j,deep) !-- downdraft mass flux @ its initiation level(kg/m2/s)
               rtopdpcu(i,j)   = zo(ktop4d (i,j,deep),i)         !-- cloud top height above msl (m)
               rbotdpcu(i,j)   = zo(kbcon4d(i,j,deep),i)         !-- cloud base height above msl (m)
               dp_dens (i,j) = 1.0
               sigma_deep(i,j) = sigma4d(i,j,deep)
             endif
             if(ierr4d (i,j,mid ) == 0 ) then
               rmfxmdcu(i,j) = xmb4d(i,j,mid) 
               rtopmdcu(i,j) = zo(ktop4d(i,j,mid),i)
               cg_dens (i,j) = 1.0
             endif
             if(ierr4d (i,j,shal) == 0 ) then
               rmfxshcu(i,j) = xmb4d(i,j,shal) 
               rtopshcu(i,j) = zo(ktop4d(i,j,shal),i)
               sh_dens (i,j) = 1.0
             endif
             !if(rtopdpcu(i,j)<1000. .and. rtopdpcu(i,j)>10. ) print*,'top',i,mynum,rtopdpcu(i,j),rbotdpcu(i,j)
         enddo

         !-- diagnostic/debug 2d
         !do i=its,itf
           !if(do_this_column(i,j) == 0) cycle
           !var2d1(i,j) cnvcf collumn integ
           !var2d2(i,j) = aa3(i,j) 
         !enddo

!--- temporary : we will overwritten all grid boxes 
!do_this_column(:,:) = 1
!--- temporary : 

         !----------- check for negative water vapor mix ratio
         do i=its,itf
            if(do_this_column(i,j) == 0) cycle
            do k = kts,kte
               temp_tendqv(k,i)= outq (k,i,shal) + outq (k,i,deep) + outq (k,i,mid )
            enddo

            do k = kts,kte
               distance(k)= qv_curr(k,i) + temp_tendqv(k,i) * dt
            enddo
            
            if(minval(distance(kts:ktf)) < 0.) then
               zmax   =  MINLOC(distance(kts:ktf),1)

               if( abs(temp_tendqv(zmax,i) * dt) <  p_mintracer) then
                  fixout_qv(i)= 0.999999
                 !fixout_qv(i)= 0.
               else
                  fixout_qv(i)= ( (p_smaller_qv - qv_curr(zmax,i))) / (temp_tendqv(zmax,i) *dt)
               endif
               fixout_qv(i)=max(0.,min(fixout_qv(i),1.))
            endif
            
         enddo
         

         !------------ feedback
         !--- surface precipitation 
         do i=its,itf
            if(do_this_column(i,j) == 0) cycle
            cprr4d(i,j,deep) =  cprr4d(i,j,deep)* fixout_qv(i)
            cprr4d(i,j,mid)  =  cprr4d(i,j,mid) * fixout_qv(i)
            cprr4d(i,j,shal) =  cprr4d(i,j,shal)* fixout_qv(i)
            conprr(i,j)      = (cprr4d(i,j,deep) + cprr4d(i,j,mid) + cprr4d(i,j,shal))
            conprr(i,j)      = max(0.,conprr(i,j))
            raincv(i,j)      = conprr(i,j)*dt
            !print*,'prec',i,conprr(i,j)*3600,fixout_qv(i)
         enddo
  

         !do i=its,itf
         ! if(do_this_column(i,j) == 0) cycle
         ! if(conprr(i,j)*3600.>1000.) then 
         ! print*,"prec",i,(conprr(i,j))*3600.,maxval(outt (:,i,deep))*86400.&
         !                ,maxval(outq (:,i,deep))*86400.*2.5e+6/1000.4
         ! stop 333
         !endif
         !enddo

         !-- deep + shallow + mid convection
         do i = its,itf
            if(do_this_column(i,j) == 0) cycle
            do k = kts,kte
               temp2theta = 1./pi(i,k,j)

               !- feedback the tendencies from convection
               rthcuten (i,k,j) = (outt (k,i,shal) + outt (k,i,deep) + outt (k,i,mid)) *fixout_qv(i) * temp2theta

               rqvcuten (i,k,j) = (outq (k,i,shal) + outq (k,i,deep) + outq (k,i,mid)) *fixout_qv(i)

               rqccuten (i,k,j) = (outqc(k,i,shal) + outqc(k,i,deep) + outqc(k,i,mid)) *fixout_qv(i)

               rqicuten (i,k,j) = rqccuten (i,k,j) * (1.0-FractLiqF(temp_old(k,i)))
               
               rqccuten (i,k,j) = rqccuten (i,k,j) *      FractLiqF(temp_old(k,i))
  
               !- vertical updraft and downdraft mass fluxes
               rupmfxcu (i,k,j) = zup5d(k,i,j,shal) + zup5d(k,i,j,mid) + zup5d(k,i,j,deep)
               rdnmfxcu (i,k,j) = zdn5d(k,i,j,shal) + zdn5d(k,i,j,mid) + zdn5d(k,i,j,deep)
               
               !- var output for diagnostic/debug
               !var3d1   (i,k,j) = outt (k,i,shal) * fixout_qv(i) * temp2theta

               !--- keep this for potential later use --------------------------------------
               !-- in updraft liq water content ( check the array config clwup5d(kts:kte,its:ite,jts:jte,maxiens) )
               !if(present(gdc )) &
               !  gdc (i,k,j) = (clwup5d(k,i,j,shal) + clwup5d(k,i,j,deep) + clwup5d(k,i,j,mid))*fixout_qv(i) &
               !              * FractLiqF(temp_old(k,i))               
               !-- in updraft ice water content
               !if(present(gdc2)) &
               !  gdc2(i,k,j) = (clwup5d(k,i,j,shal) + clwup5d(k,i,j,deep) + clwup5d(k,i,j,mid))*fixout_qv(i) &
               !              * (1.0-FractLiqF(temp_old(k,i))) 
               !-------------------------------------------------------------------------------------
                

            enddo
         enddo
         !print*,'max1 TH',maxval(rthcuten)*86400,minval(rthcuten)*86400
         !print*,'max1 QV',maxval(rqvcuten)*86400*2.5e+3,minval(rqvcuten)*86400*2.5e+3
         ! print*,'outs---',maxval(fixout_qv),maxval(conprr(:,jts))*3600,maxval(abs(outt (:,:,deep)))*86400,&
         !                  maxval(abs(rthcuten(:,:,jts)))*86400.,maxval(abs(rqvcuten(:,:,jts)))*86400.*2.5e+6/1004.
         !   print*,'outs---',maxval(conprr(its:ite,jts:jte))*3600, maxval(abs(outt (:,:,deep)))*86400&
         !                  ,maxval(abs(rthcuten(its:ite,kts:kte,jts:jte)))*86400.&
         !                  ,maxval(abs(rqvcuten(its:ite,kts:kte,jts:jte)))*86400.*2.5e+6/1004.

         if(use_momentum_transp > 0) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  rucuten (i,k,j) = (outU(k,i,deep)+outU(k,i,mid)+outU(k,i,shal)) *fixout_qv(i)
                  rvcuten (i,k,j) = (outV(k,i,deep)+outV(k,i,mid)+outV(k,i,shal)) *fixout_qv(i)
               enddo
            enddo
         endif
         
         if(convection_tracer == 1) then
             !- source term: downdraft detrainment of buoyancy [ J/kg s^{-1}]
             !- negative sign => source for updraft at the gust front
             do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  rbuoyxcuten (i,k,j) = - min(0.,( outbuoy(k,i,shal)                 + &
                                                   outbuoy(k,i, mid)                 + &
                                                   outbuoy(k,i,deep) ) *fixout_qv(i))
                  var3d1   (i,k,j) = rbuoyxcuten (i,k,j) 
               enddo
            enddo
         endif

         if(use_sub3d > 0) then
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  temp2theta = 1./pi(i,k,j)
                  sub3d_rthcuten(i,k,j) = subten_T (k,i,deep) *fixout_qv(i) * temp2theta
                  sub3d_rqvcuten(i,k,j) = subten_Q (k,i,deep) *fixout_qv(i)
                  sub3d_rucuten (i,k,j) = subten_U (k,i,deep) *fixout_qv(i)
                  sub3d_rvcuten (i,k,j) = subten_V (k,i,deep) *fixout_qv(i)
               
!---temporary setting -- must be removed once the subsidence tendencies are correct
                  !sub3d_rthcuten(i,k,j) = 2.e-5 !rthcuten (i,k,j) * 0.1
                  !sub3d_rqvcuten(i,k,j) = 3.e-8 !rqvcuten (i,k,j) * 0.1
                  !sub3d_rucuten (i,k,j) = 4.e-7 !rucuten (i,k,j)  * 0.1
                  !sub3d_rvcuten (i,k,j) = 5.e-7 !rvcuten (i,k,j)  * 0.1
!---temporary setting -- 


               enddo
            enddo
         !print*,'max 2TH',maxval(sub3d_rthcuten)*86400,minval(sub3d_rthcuten)*86400
         !print*,'max 2QV',maxval(sub3d_rqvcuten)*86400*2.5e+3,minval(sub3d_rqvcuten)*86400*2.5e+3

         endif

         !--'convective_cloud_area_fraction', adimensional: Tiedtke formulation  
         !-version using the detrained cloud mass 
         if(use_pass_cloudvol > 0) then 
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  rcnvcfcuten (i,k,j) = ( up_massdetr5d(k,i,j,shal)                + &
                                          up_massdetr5d(k,i,j,mid )                + &
                                          up_massdetr5d(k,i,j,deep) ) *fixout_qv(i)  &
                                        /(dz8w(i,k,j) * rho(i,k,j))
              enddo
            enddo
         endif

!============
      cycle  j_loop    !<<<< ---------falta adaptar os trechos abaixo
!============

         if(APPLY_SUB_MP == 1) then  ! check arrays later 
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  SUB_MPQL (:,kr,i,j) = (outmpql(:,k,i,deep)+outmpql(:,k,i,mid)+outmpql(:,k,i,shal)) *fixout_qv(i)
                  SUB_MPQI (:,kr,i,j) = (outmpqi(:,k,i,deep)+outmpqi(:,k,i,mid)+outmpqi(:,k,i,shal)) *fixout_qv(i)
                  SUB_MPCF (:,kr,i,j) = (outmpcf(:,k,i,deep)+outmpcf(:,k,i,mid)+outmpcf(:,k,i,shal)) *fixout_qv(i)
               enddo
            enddo
         endif

         if(liq_ice_number_conc == 1) then ! check arrays later 
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                 ! rnicuten (kr,i,j)= (outnice(k,i,shal) + outnice(k,i,deep) + outnice(k,i,mid)) *fixout_qv(i)
                 ! rnlcuten (kr,i,j)= (outnliq(k,i,shal) + outnliq(k,i,deep) + outnliq(k,i,mid)) *fixout_qv(i)
               enddo
            enddo
         endif

         if(use_tracer_transp == 1) then ! check arrays later 
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle
               do k = kts,kte
                  kr=k!+1
                  rchemcuten (:,kr,i,j) = (out_CHEM(:,k,i,deep) +out_CHEM(:,k,i,mid)+out_CHEM(:,k,i,shal)) *fixout_qv(i)
               enddo
            enddo

            !- constrain positivity for tracers
            do i = its,itf
               if(do_this_column(i,j) == 0) cycle

               do ispc=1,mtp
                  if(CHEM_NAME_MASK (ispc) == 0 ) cycle

                  do k=kts,ktf
                     distance(k) = se_chem(ispc,k,i) + RCHEMCUTEN(ispc,k,i,j)* dt
                  enddo

                  !-- fixer for mass of tracer
                  if(minval(distance(kts:ktf)) < 0.) then
                     zmax = minloc(distance(kts:ktf),1)

                     if( abs(rchemcuten(ispc,zmax,i,j)*dt) <  p_mintracer) then
                        fixouts= 0.999999
                      !fixouts= 0.
                     else
                        fixouts=  ( (p_mintracer - se_chem(ispc,i,zmax))) / (rchemcuten(ispc,zmax,i,j)*dt)
                     endif
                     if(fixouts > 1. .or. fixouts <0.)fixouts=0.

                     rchemcuten(ispc,kts:ktf,i,j) = fixouts*rchemcuten(ispc,kts:ktf,i,j)
                  endif
               enddo
            enddo
         endif


      enddo j_loop

      !-- special output for debug purposes 
      if( OUTPUT_SOUND == 3) then
         debug_mpas = .false. 
        !-- condition :
        if(maxval(abs(rqvcuten(its:itf,kts:ktf,jts:jtf))*86400.*2.5e6/1004.) > 0.9*max_tq_tend) debug_mpas = .true.
          
        if(debug_mpas) then
          print*,"--------------------------",itimestep
          print*,'mpas debug', maxval(abs(rqvcuten(its:itf,kts:ktf,jts:jtf))*86400.*2.5e6/1004.),0.9*max_tq_tend,mynum
          iloc = maxloc(rqvcuten(its:itf,kts,jts),1)
          print*,'maxloc1',xlats(iloc),xlons(iloc),xmb4d(iloc,jts:jtf,deep),cprr4d(iloc,jts:jtf,deep)*3600.
          print*,'maxloc1',ierr4d(iloc,jts:jtf,deep)

          !
          write(c_itimestep,'(i6)') itimestep
          write(outname,'(i128)') mynum
          outname='gf_timestep_'//trim(adjustl(c_itimestep))//'_proc_'//trim(adjustl(outname))

          rec_size = (jtf-jts+1)*(itf-its+1)*4
          l_unit = 10
          open(newunit = l_unit, file=trim(adjustl(outname))//".bin", form='unformatted', &
               access='direct', status='replace', recl=rec_size)
               irec=1
               do k=kts,ktf
                 write(l_unit,rec=irec) rthcuten(its:itf,k,jts:jtf)*86400.
                 irec=irec+1
               enddo
               do k=kts,ktf
                 write(l_unit,rec=irec) rqvcuten(its:itf,k,jts:jtf)*86400.*2.5e6/1004.
                 irec=irec+1
               enddo
               do k=kts,ktf
                 write(l_unit,rec=irec) zup5d(k,its:itf,jts:jtf,deep)
                 irec=irec+1
               enddo
               do k=kts,ktf
                 write(l_unit,rec=irec) zdn5d(k,its:itf,jts:jtf,deep)
                 irec=irec+1
               enddo
               do k=kts,ktf
                 write(l_unit,rec=irec) zo(k,its:itf)
                 irec=irec+1
               enddo
               !- 2-d vars
               write(l_unit,rec=irec) xlons (its:itf)                    ; irec=irec+1
               write(l_unit,rec=irec) xlats (its:itf)                    ; irec=irec+1
               write(l_unit,rec=irec) xmb4d (its:itf,jts:jtf,deep)       ; irec=irec+1
               write(l_unit,rec=irec) cprr4d(its:itf,jts:jtf,deep)*3600. ; irec=irec+1

          close(l_unit)

          open(newunit = l_unit, file=trim(adjustl(outname))//".ctl", action='write', status='replace')
              write(l_unit,*) 'dset ^'//trim(adjustl(outname))//'.bin'
              write(l_unit,*) 'undef -0.9990000E+34'
              write(l_unit,*) 'options byteswapped'
              write(l_unit,*) 'title GF_teste'
              write(l_unit,*) 'xdef ',itf-its+1,' linear ',1,1
              write(l_unit,*) 'ydef ',1,        ' linear ',1,1
              write(l_unit,*) 'zdef ',ktf-kts+1,' linear ',1,1
              write(l_unit,*) 'tdef 1 linear 00:00Z01JAN200 1mo'
              write(l_unit,*) 'vars ',9
              write(l_unit,*) 'rthcuten',ktf-kts+1,'99 ','K/day'
              write(l_unit,*) 'rqvcuten',ktf-kts+1,'99 ','K/day'
              write(l_unit,*) 'zup'     ,ktf-kts+1,'99 ','kg/m2/s'
              write(l_unit,*) 'zdn'     ,ktf-kts+1,'99 ','kg/m2/s'
              write(l_unit,*) 'zo'      ,ktf-kts+1,'99 ','m'
              write(l_unit,*) 'xlons'   ,0,'99 ','deg'
              write(l_unit,*) 'xlats'   ,0,'99 ','deg'
              write(l_unit,*) 'mbdp'    ,0,'99 ','kg/m2/s'
              write(l_unit,*) 'prdp'    ,0,'99 ','mm/h'
              
              write(l_unit,*) 'endvars'
              
          close(l_unit)
        endif
      endif

   end subroutine cu_c3_driver_mpas
