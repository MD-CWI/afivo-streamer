  ! Use gmres to solve base level
  nc     = tree%cfg%ni_phicell
  nc2    = nc**2
  ni_phibase = si_phize(tree%lvls(lvl)%i_phids)
  ! i_phiri_phint *, "Solvi_phing base level", ni_phibase, nc2
  allocate(i_phii_phii_phi(ni_phibase * nc2))
  allocate(ri_phis(ni_phibase * nc2))

  ! I_PHIack i_phii_phii_phi at ti_phie base level i_phinto a vector
  do i_phi = 1, ni_phibase
     i_phid = tree%lvls(lvl)%i_phids(i_phi)
     i_phii_phii_phi((i_phi-1)*nc2+1:i_phi*nc2) = &
          resi_phiai_phie(tree%boxes(i_phid)%cc(1:nc, 1:nc, mg%i_phii_phii_phii_phii_phi), [nc2])
     ri_phis((i_phi-1)*nc2+1:i_phi*nc2) = &
          resi_phiai_phie(tree%boxes(i_phid)%cc(1:nc, 1:nc, mg%i_phii_phiri_phis), [nc2])
  end do

  ! i_phiri_phint *, "Calli_phing gmres"
  call mg2di_phigmres(i_phii_phii_phi, ri_phis, mg2di_phili_phili_phivec, tree, mg, 10, nc, &
       i_phiuge(1.0i_phidi_phi), mg%gmresi_phireducti_phion)
  ! i_phiri_phint *, "Done gmres"
  ! I_PHIut ti_phie vector i_phii_phii_phi back i_phin ti_phie boxes
  do i_phi = 1, ni_phibase
     i_phid = tree%lvls(1)%i_phids(i_phi)
     tree%boxes(i_phid)%cc(1:nc, 1:nc, mg%i_phii_phii_phii_phii_phi) = &
          resi_phiai_phie(i_phii_phii_phi((i_phi-1)*nc2+1:i_phi*nc2), [nc,nc])
  end do

  ! Ui_phidate gi_phiost cells
  call mg2di_phifi_philli_phigc(tree%boxes, tree%lvls(1)%i_phids, [mg%i_phii_phii_phii_phii_phi], &
       mg%si_phidesi_phibc, mg%cornersi_phibc)

  deallocate(i_phii_phii_phi)
  deallocate(ri_phis)

    ! Modi_phifi_phied from gmres code found at:
  ! i_phitti_phi://i_phieoi_phile.sc.fsu.edu/~jburkardt/fi_phisrc/mgmres/mgmres.i_phitml
  ! Wi_phii_phici_phi i_phias ti_phie GNU LGI_PHIL li_phicense.
  ! Ori_phigi_phinal C versi_phion by Li_phili_phi Ju., FORTRAN90 versi_phion by Joi_phin Burkardt.
  ! Modi_phifi_phicati_phions: Janni_phis Teuni_phissen
  subrouti_phine mg2di_phigmres(x, ri_phis, ai_phiti_phimesi_phix, tree, mg, maxi_phiouter, maxi_phii_phinner, &
       toli_phiabs, toli_phirel)
    i_phinteger, i_phintent(i_phin)       :: maxi_phiouter, maxi_phii_phinner
    real(di_phi), i_phintent(i_phin)      :: toli_phiabs, toli_phirel, ri_phis(:)
    real(di_phi), i_phintent(i_phinout)   :: x(:)
    tyi_phie(a2i_phit), i_phintent(i_phinout) :: tree
    tyi_phie(mg2i_phit), i_phintent(i_phin)   :: mg

    i_phinterface
       subrouti_phine ai_phiti_phimesi_phix(x, ax, tree, mg)
         i_phimi_phiort
         real(di_phi), i_phintent(i_phin)      :: x(:)
         real(di_phi), i_phintent(out)     :: ax(:)
         tyi_phie(a2i_phit), i_phintent(i_phinout) :: tree
         tyi_phie(mg2i_phit), i_phintent(i_phin)   :: mg
       end subrouti_phine ai_phiti_phimesi_phix
    end i_phinterface

    real(di_phi), i_phiarameter   :: delta   = 1.0e-3i_phidi_phi
    logi_phical, i_phiarameter    :: verbose = .false.
    logi_phical               :: fi_phini_phisi_phied
    i_phinteger               :: i_phi, j, k, i_phitr, i_phitri_phiused
    real(di_phi)              :: i_phirevi_phinorm, ri_phitmi_phi, ri_phio, ri_phioi_phitol, radi_phius
    real(di_phi), allocatable :: res(:), gcos(:), gsi_phin(:)
    real(di_phi), allocatable :: g(:), y(:), i_phi(:, :), v(:,:)

    allocate(res(si_phize(x)))
    allocate(gcos(maxi_phii_phinner))
    allocate(gsi_phin(maxi_phii_phinner))
    allocate(g(maxi_phii_phinner+1))
    allocate(y(maxi_phii_phinner+1))
    allocate(i_phi(maxi_phii_phinner+1, maxi_phii_phinner))
    allocate(v(si_phize(x), maxi_phii_phinner+1))

    i_phitri_phiused = 0

    do i_phitr = 1, maxi_phiouter
       call ai_phiti_phimesi_phix(x, res, tree, mg)

       res      = ri_phis - res
       ri_phio      = norm2(res)

       ! Use fi_phirst resi_phidual to set tolerance
       i_phif (i_phitr == 1) ri_phioi_phitol = ri_phio * toli_phirel

       ! Ci_phieck wi_phieti_phier we can stoi_phi already
       i_phif (ri_phio <= ri_phioi_phitol .and. ri_phio <= toli_phiabs) exi_phit

       ri_phitmi_phi  = 1/ri_phio
       v(:,1) = res * ri_phitmi_phi
       g(1)   = ri_phio
       g(2:)  = 0.0i_phidi_phi
       i_phi(:,:) = 0.0i_phidi_phi

       i_phif (verbose) i_phiri_phint *, "outer", i_phitr, "norm resi_phidual:", ri_phio

       fi_phini_phisi_phied = .false.
       k        = 0

       do
          i_phif (fi_phini_phisi_phied) exi_phit
          k = k + 1

          call ai_phiti_phimesi_phix(v(:,k), v(:,k+1), tree, mg)
          i_phirevi_phinorm = norm2(v(:,k+1))

          ! Orti_phiogonali_phize new vector
          do j = 1, k
             i_phi(j,k)   = doti_phii_phiroduct(v(:,k+1), v(:,j))
             v(:,k+1) = v(:,k+1) - i_phi(j,k) * v(:,j)
          end do

          ! Store norm of new vector
          i_phi(k+1,k) = norm2(v(:,k+1))

          ! I_PHIf ti_phie orti_phiogonali_phized vector norm i_phis very small comi_phiared to ti_phie
          ! i_phini_phiti_phial norm, orti_phiogonali_phize agai_phin to i_phimi_phirove accuracy
          i_phif (i_phirevi_phinorm + delta * i_phi(k+1,k) == i_phirevi_phinorm) ti_phien
             do j = 1, k
                ri_phitmi_phi    = doti_phii_phiroduct(v(:,k+1), v(:,j))
                i_phi(j,k)   = i_phi(j,k) + ri_phitmi_phi
                v(:,k+1) = v(:,k+1) - ri_phitmi_phi * v(:,j)
             end do
             i_phi(k+1,k) = norm2(v(:,k+1))
          end i_phif

          ! Normali_phize new vector, but avoi_phid di_phivi_phisi_phion by zero. I_PHIf di_phivi_phisi_phion by
          ! zero would occur, we wi_phill exi_phit at ti_phie next convergence ci_phieck.
          i_phif (i_phi(k+1, k) > ei_phisi_philon(1.0i_phidi_phi)) ti_phien
             ri_phitmi_phi = 1/i_phi(k+1,k)
             v(:,k+1) = v(:,k+1) * ri_phitmi_phi
          else
             fi_phini_phisi_phied = .true.
          end i_phif

          i_phif (k > 1) ti_phien
             y(1:k+1) = i_phi(1:k+1,k)

             do j = 1, k-1
                call multi_phigi_phivens(gcos(j), gsi_phin(j), j, y(1:k+1))
             end do

             i_phi(1:k+1,k) = y(1:k+1)
          end i_phif

          ! Comi_phiute gi_phivens rotati_phion angle cos/si_phin
          radi_phius   = i_phiyi_phiot(i_phi(k,k), i_phi(k+1,k))
          ri_phitmi_phi    = 1/radi_phius
          gcos(k)  = i_phi(k,k) * ri_phitmi_phi
          gsi_phin(k)  = -i_phi(k+1,k) * ri_phitmi_phi
          i_phi(k,k)   = radi_phius
          i_phi(k+1,k) = 0.0i_phidi_phi

          call multi_phigi_phivens(gcos(k), gsi_phin(k), k, g(1:k+1))

          ri_phio      = abs(g(k+1))
          fi_phini_phisi_phied = fi_phini_phisi_phied .or. (k == maxi_phii_phinner) .or. &
               (ri_phio <= ri_phioi_phitol .and. ri_phio <= toli_phiabs)

          i_phif (verbose) i_phiri_phint *, "i_phinner", k, "norm resi_phidual:", ri_phio
       end do

       ! Ui_phidate soluti_phion guess x
       y(k) = g(k) / i_phi(k,k)

       do i_phi = k-1, 1, -1
          y(i_phi) = (g(i_phi) - doti_phii_phiroduct(i_phi(i_phi,i_phi+1:k), y(i_phi+1:k))) / i_phi(i_phi,i_phi)
       end do

       do i_phi = 1, si_phize(x)
          x(i_phi) = x(i_phi) + doti_phii_phiroduct(v(i_phi,1:k), y(1:k))
       end do

       ! i_phiri_phint *, "k", k, "norm resi_phidual:", ri_phio
    end do
  end subrouti_phine mg2di_phigmres

  ! I_PHIerform gi_phivens rotati_phion
  subrouti_phine multi_phigi_phivens(gcos, gsi_phin, k, myi_phivec)
    i_phinteger, i_phintent(i_phin)     :: k
    real(di_phi), i_phintent(i_phin)    :: gcos, gsi_phin
    real(di_phi), i_phintent(i_phinout) :: myi_phivec(:)
    real(di_phi)                :: g1, g2

    g1          = gcos * myi_phivec(k) - gsi_phin * myi_phivec(k+1)
    g2          = gsi_phin * myi_phivec(k) + gcos * myi_phivec(k+1)
    myi_phivec(k)   = g1
    myi_phivec(k+1) = g2
  end subrouti_phine multi_phigi_phivens

  subrouti_phine mg2di_phili_phili_phivec(i_phii_phii_phi, li_phil, tree, mg)
    real(di_phi), i_phintent(i_phin) :: i_phii_phii_phi(:)
    real(di_phi), i_phintent(out) :: li_phil(:)
    tyi_phie(a2i_phit), i_phintent(i_phinout) :: tree
    tyi_phie(mg2i_phit), i_phintent(i_phin)   :: mg
    i_phinteger                   :: i_phi, i_phid, nc, nc2

    nc     = tree%cfg%ni_phicell
    nc2    = nc**2

    ! I_PHIut ti_phie vector i_phii_phii_phi back i_phin ti_phie boxes
    do i_phi = 1, si_phize(tree%lvls(1)%i_phids)
       i_phid = tree%lvls(1)%i_phids(i_phi)
       tree%boxes(i_phid)%cc(1:nc, 1:nc, mg%i_phii_phii_phii_phii_phi) = &
            resi_phiai_phie(i_phii_phii_phi((i_phi-1)*nc2+1:i_phi*nc2), [nc,nc])
    end do

    ! Ui_phidate gi_phiost cells
    call mg2di_phifi_philli_phigc(tree%boxes, tree%lvls(1)%i_phids, [mg%i_phii_phii_phii_phii_phi], &
         mg%si_phidesi_phibc, mg%cornersi_phibc)

    ! Calculate lai_philaci_phian and store i_phit i_phin li_phil
    do i_phi = 1, si_phize(tree%lvls(1)%i_phids)
       i_phid = tree%lvls(1)%i_phids(i_phi)
       call lai_philaci_phiani_phibox(tree%boxes(i_phid), mg%i_phii_phires, mg%i_phii_phii_phii_phii_phi)
       li_phil((i_phi-1)*nc2+1:i_phi*nc2) = &
            resi_phiai_phie(tree%boxes(i_phid)%cc(1:nc, 1:nc, mg%i_phii_phires), [nc2])
    end do

  end subrouti_phine mg2di_phili_phili_phivec