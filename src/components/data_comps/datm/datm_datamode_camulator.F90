module datm_datamode_camulator_mod

  !---------------------------------------------------------------------------
  ! datm_datamode_camulator.F90
  !
  ! CAMulator DATM mode — file-based coupling between CESM2/CPL7 and the
  ! CAMulator Python AI atmosphere (Chapman et al. 2025).
  !
  ! Called once per 6-hour coupling step from datm_comp_mod.F90:
  !
  !   case('CAMULATOR')
  !      call datm_datamode_camulator_run(x2a, a2x, ...)
  !
  ! Each coupling step:
  !   1. Gather SST (So_t) and ice fraction (Sf_ifrac) from x2a to master task
  !   2. Master writes camulator_sst_in.nc   (So_t, Sf_ifrac, YMD, TOD)
  !   3. Master writes camulator_go.flag      → signals Python server
  !   4. Master polls  camulator_done.flag   ← Python has finished inference
  !   5. Master reads  camulator_cam_out.nc  (winds, T, Q, P, SW, LW, precip)
  !   6. Broadcast results to all MPI ranks
  !   7. Each rank populates its a2x partition
  !
  ! Sign / unit conventions (see README_Coupling.md):
  !   - CAMulator TAUX/TAUY are provided as 10m winds (Sa_u, Sa_v).
  !     The CPL7 bulk formula computes wind stress from these.
  !     (Phase 2: bypass bulk formula by directly setting Faxx_taux/tauy)
  !   - CAMulator LW output is net LW at surface (FLNS, positive upward).
  !     We split into downward LW:  Faxa_lwdn = FLNSD (downward component).
  !     If Python provides only net LW, estimate FLNSD = FLNS_net + sigma*Ts^4.
  !   - SW: Python server reconstructs FSDS (downwelling SW) from FSNS via
  !     FSDS = FSNS / (1-alpha_sfc), where alpha_sfc uses CICE ice fraction.
  !     We split FSDS into direct/diffuse components using fixed fractions.
  !     The coupler (seq_flux_mct.F90) then applies ocean/ice albedo to get
  !     absorbed SW — this is the correct sign/magnitude convention.
  !   - Precip: PRECT in m/s liquid-water equivalent. Split rain/snow by
  !     temperature and convert to kg m-2 s-1 (multiply by rho_w=1000).
  !
  ! Timeout: if camulator_done.flag does not appear within POLL_MAX_ITER
  !          iterations (default 3600 x 1s = 1 hour), abort with clear message.
  !
  ! NetCDF variables in camulator_sst_in.nc:
  !   sst   (ngrid)  : ocean surface temp from CPL [K]
  !   ifrac (ngrid)  : sea-ice fraction from CPL [0-1]
  !   ymd            : current date YYYYMMDD (int)
  !   tod            : current time-of-day in seconds (int)
  !
  ! NetCDF variables in camulator_cam_out.nc (written by Python):
  !   u10   (ngrid)  : 10m zonal wind [m/s]       -> Sa_u
  !   v10   (ngrid)  : 10m meridional wind [m/s]   -> Sa_v
  !   tbot  (ngrid)  : near-surface temperature [K]-> Sa_tbot, Sa_ptem
  !   qbot  (ngrid)  : specific humidity [kg/kg]   -> Sa_shum
  !   pbot  (ngrid)  : near-surface pressure [Pa]  -> Sa_pbot, Sa_pslv
  !   fsds  (ngrid)  : downwelling SW [W m-2]        -> Faxa_sw* (split; coupler applies albedo)
  !   flnsd (ngrid)  : downward LW at surface [W m-2]-> Faxa_lwdn
  !   prect (ngrid)  : total precip [m s-1 liq eq] -> Faxa_rainl / Faxa_snowl
  !
  ! Author: Claude Code  (generated 2026-02-22)
  ! See also: climate/README_Coupling.md
  !---------------------------------------------------------------------------

  use shr_kind_mod , only : IN=>SHR_KIND_IN, R8=>SHR_KIND_R8, &
                             CL=>SHR_KIND_CL, CS=>SHR_KIND_CS
  use shr_sys_mod  , only : shr_sys_abort, shr_sys_flush
  use shr_mpi_mod  , only : shr_mpi_bcast
  use shr_const_mod, only : SHR_CONST_TKFRZ, SHR_CONST_STEBOL
  use shr_file_mod , only : shr_file_getUnit, shr_file_freeUnit
  use mct_mod
  use netcdf
  use mpi

  implicit none
  private

  !--- public interface ---
  public :: datm_datamode_camulator_run

  !--- flag file names (written/read relative to CESM run directory / cwd) ---
  character(len=*), parameter :: GO_FLAG   = 'camulator_go.flag'
  character(len=*), parameter :: DONE_FLAG = 'camulator_done.flag'
  character(len=*), parameter :: SST_FILE  = 'camulator_sst_in.nc'
  character(len=*), parameter :: CAM_FILE  = 'camulator_cam_out.nc'

  !--- polling parameters ---
  integer(IN), parameter :: POLL_MAX_ITER  = 3600   ! iterations before abort
  integer(IN), parameter :: POLL_SLEEP_SEC = 1      ! seconds per poll cycle

  !--- physical constants ---
  real(R8), parameter :: tKFrz  = SHR_CONST_TKFRZ  ! 273.15 K
  real(R8), parameter :: stebol = SHR_CONST_STEBOL  ! Stefan-Boltzmann [W m-2 K-4]
  real(R8), parameter :: rho_w  = 1000.0_R8         ! liquid water density [kg m-3]

  !--- SW partitioning fractions (sum to 1.0; from CORE2 DATM convention) ---
  real(R8), parameter :: sw_frac_swvdr = 0.28_R8
  real(R8), parameter :: sw_frac_swndr = 0.31_R8
  real(R8), parameter :: sw_frac_swvdf = 0.24_R8
  real(R8), parameter :: sw_frac_swndf = 0.17_R8

  !--- safety clamps ---
  real(R8), parameter :: WIND_MAX =   80.0_R8    ! m s-1
  real(R8), parameter :: QBOT_MIN =   1.0e-9_R8  ! kg kg-1 (avoid negative q)

  !--- module-level state (set on firstcall, reused each coupling step) ---
  logical,  save :: cam_initialized = .false.
  integer(IN), save :: lsize_global = 0   ! total ATM grid points (nxg*nyg)
  integer(IN), save :: npes_save    = 0   ! number of MPI ranks (cached)

  !--- x2a field indices (set on firstcall) ---
  integer(IN), save :: ksst   = -1   ! So_t
  integer(IN), save :: kifrac = -1   ! Sf_ifrac

  !--- a2x field indices (set on firstcall) ---
  integer(IN), save :: ku, kv, ktbot, kptem, kshum, kdens
  integer(IN), save :: kpbot, kpslv, kz
  integer(IN), save :: klwdn
  integer(IN), save :: kswvdr, kswndr, kswvdf, kswndf, kswnet
  integer(IN), save :: krl, krc, ksl, ksc

  !--- global gathered/broadcast arrays (master only for input; all for output) ---
  real(R8), allocatable, save :: g_sst(:)     ! global SST [K]
  real(R8), allocatable, save :: g_ifrac(:)   ! global ice fraction [0-1]
  real(R8), allocatable, save :: g_u10(:)     ! zonal wind at bottom model level [m s-1]
  real(R8), allocatable, save :: g_v10(:)     ! meridional wind at bottom model level [m s-1]
  real(R8), allocatable, save :: g_tbot(:)    ! temperature at bottom model level [K]
  real(R8), allocatable, save :: g_zbot(:)    ! height of bottom model level midpoint [m]
  real(R8), allocatable, save :: g_qbot(:)    ! specific humidity at bottom model level [kg/kg]
  real(R8), allocatable, save :: g_pbot(:)    ! surface pressure [Pa]
  real(R8), allocatable, save :: g_fsds(:)    ! downwelling SW [W m-2] (reconstructed from FSNS in Python)
  real(R8), allocatable, save :: g_flnsd(:)   ! downward LW [W m-2]
  real(R8), allocatable, save :: g_prect(:)   ! total precip [m s-1]

  !--- per-rank gather/scatter work arrays (allocated on firstcall) ---
  integer(IN), allocatable, save :: recvcounts(:)  ! MPI_Gatherv recv counts
  integer(IN), allocatable, save :: displs(:)      ! MPI_Gatherv displacements

  save

!===============================================================================
CONTAINS
!===============================================================================

  subroutine datm_datamode_camulator_run( &
       x2a, a2x, ggrid, gsmap, &
       mpicom, my_task, master_task, logunit, &
       currentYMD, currentTOD, nxg, nyg)

    !---------------------------------------------------------------------------
    ! Arguments
    !---------------------------------------------------------------------------
    type(mct_aVect), intent(inout) :: x2a        ! coupler -> DATM fields
    type(mct_aVect), intent(inout) :: a2x        ! DATM -> coupler fields
    type(mct_gGrid), intent(in)    :: ggrid      ! ATM grid
    type(mct_gsMap), intent(in)    :: gsmap      ! ATM global seg map
    integer(IN),     intent(in)    :: mpicom     ! MPI communicator
    integer(IN),     intent(in)    :: my_task    ! this MPI rank
    integer(IN),     intent(in)    :: master_task! rank 0 (master)
    integer(IN),     intent(in)    :: logunit    ! log file unit
    integer(IN),     intent(in)    :: currentYMD ! YYYYMMDD
    integer(IN),     intent(in)    :: currentTOD ! seconds since midnight
    integer(IN),     intent(in)    :: nxg, nyg   ! global ATM grid dimensions

    !---------------------------------------------------------------------------
    ! Local variables
    !---------------------------------------------------------------------------
    integer(IN) :: n, lsize, ierr, npes
    integer(IN) :: lsize_local
    real(R8), allocatable :: local_sst(:), local_ifrac(:)
    real(R8), allocatable :: local_u10(:), local_v10(:), local_tbot(:), local_zbot(:)
    real(R8), allocatable :: local_qbot(:), local_pbot(:)
    real(R8), allocatable :: local_fsds(:), local_flnsd(:), local_prect(:)
    real(R8)  :: rain_kg, snow_kg, swdn, dens
    logical   :: flag_exists
    integer(IN) :: poll_iter
    character(*), parameter :: subname = '(datm_datamode_camulator_run) '
    character(*), parameter :: F00 = "('(datm_cam) ',a)"
    character(*), parameter :: F01 = "('(datm_cam) ',a,i12)"
    character(*), parameter :: F02 = "('(datm_cam) ',a,f12.4)"

    !---------------------------------------------------------------------------
    ! First-call initialisation: cache field indices and allocate global arrays
    !---------------------------------------------------------------------------
    if (.not. cam_initialized) then
       call MPI_Comm_size(mpicom, npes, ierr)
       npes_save    = npes
       lsize_global = nxg * nyg

       !--- x2a input field indices ---
       ksst   = mct_aVect_indexRA(x2a, 'So_t',    perrWith='quiet')
       kifrac = mct_aVect_indexRA(x2a, 'Sf_ifrac', perrWith='quiet')

       if (ksst < 1) then
          if (my_task == master_task) &
               write(logunit,F00) 'WARNING: So_t not in x2a — SST will default to 273.15 K'
       end if
       if (kifrac < 1) then
          if (my_task == master_task) &
               write(logunit,F00) 'WARNING: Sf_ifrac not in x2a — ice fraction will default to 0'
       end if

       !--- a2x output field indices ---
       kz    = mct_aVect_indexRA(a2x, 'Sa_z'      )
       ku    = mct_aVect_indexRA(a2x, 'Sa_u'      )
       kv    = mct_aVect_indexRA(a2x, 'Sa_v'      )
       ktbot = mct_aVect_indexRA(a2x, 'Sa_tbot'   )
       kptem = mct_aVect_indexRA(a2x, 'Sa_ptem'   )
       kshum = mct_aVect_indexRA(a2x, 'Sa_shum'   )
       kdens = mct_aVect_indexRA(a2x, 'Sa_dens'   )
       kpbot = mct_aVect_indexRA(a2x, 'Sa_pbot'   )
       kpslv = mct_aVect_indexRA(a2x, 'Sa_pslv'   )
       klwdn = mct_aVect_indexRA(a2x, 'Faxa_lwdn' )
       kswvdr= mct_aVect_indexRA(a2x, 'Faxa_swvdr')
       kswndr= mct_aVect_indexRA(a2x, 'Faxa_swndr')
       kswvdf= mct_aVect_indexRA(a2x, 'Faxa_swvdf')
       kswndf= mct_aVect_indexRA(a2x, 'Faxa_swndf')
       kswnet= mct_aVect_indexRA(a2x, 'Faxa_swnet')
       krl   = mct_aVect_indexRA(a2x, 'Faxa_rainl')
       krc   = mct_aVect_indexRA(a2x, 'Faxa_rainc')
       ksl   = mct_aVect_indexRA(a2x, 'Faxa_snowl')
       ksc   = mct_aVect_indexRA(a2x, 'Faxa_snowc')

       !--- MPI gather setup ---
       allocate(recvcounts(npes), displs(npes))

       !--- global output arrays (all ranks) ---
       allocate(g_sst  (lsize_global)); g_sst   = tKFrz
       allocate(g_ifrac(lsize_global)); g_ifrac  = 0.0_R8
       allocate(g_u10  (lsize_global)); g_u10    = 0.0_R8
       allocate(g_v10  (lsize_global)); g_v10    = 0.0_R8
       allocate(g_tbot (lsize_global)); g_tbot   = tKFrz + 15.0_R8
       allocate(g_zbot (lsize_global)); g_zbot   = 61.0_R8   ! mid-lat default [m]
       allocate(g_qbot (lsize_global)); g_qbot   = 1.0e-3_R8
       allocate(g_pbot (lsize_global)); g_pbot   = 1.01325e5_R8
       allocate(g_fsds (lsize_global)); g_fsds   = 0.0_R8
       allocate(g_flnsd(lsize_global)); g_flnsd  = 300.0_R8
       allocate(g_prect(lsize_global)); g_prect  = 0.0_R8

       cam_initialized = .true.

       if (my_task == master_task) then
          write(logunit,F00) 'CAMULATOR mode initialised'
          write(logunit,F01) '  lsize_global = ', lsize_global
          write(logunit,F01) '  npes         = ', npes
          write(logunit,F01) '  ksst (So_t)  = ', ksst
          write(logunit,F01) '  kifrac       = ', kifrac
          call shr_sys_flush(logunit)
       end if
    end if  ! firstcall

    lsize_local = mct_aVect_lsize(x2a)

    !===========================================================================
    ! 1. Gather SST and ice fraction from all ranks to master
    !===========================================================================

    allocate(local_sst  (lsize_local))
    allocate(local_ifrac(lsize_local))

    if (ksst > 0) then
       local_sst(:) = x2a%rAttr(ksst, 1:lsize_local)
    else
       local_sst(:) = tKFrz   ! freezing point fallback
    end if

    if (kifrac > 0) then
       local_ifrac(:) = x2a%rAttr(kifrac, 1:lsize_local)
    else
       local_ifrac(:) = 0.0_R8
    end if

    !--- build MPI_Gatherv displacement array ---
    call MPI_Gather(lsize_local, 1, MPI_INTEGER, &
                    recvcounts,  1, MPI_INTEGER, &
                    master_task, mpicom, ierr)

    if (my_task == master_task) then
       displs(1) = 0
       do n = 2, npes_save
          displs(n) = displs(n-1) + recvcounts(n-1)
       end do
    end if

    call MPI_Gatherv(local_sst,   lsize_local, MPI_DOUBLE_PRECISION, &
                     g_sst,  recvcounts, displs, MPI_DOUBLE_PRECISION, &
                     master_task, mpicom, ierr)

    call MPI_Gatherv(local_ifrac, lsize_local, MPI_DOUBLE_PRECISION, &
                     g_ifrac, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                     master_task, mpicom, ierr)

    deallocate(local_sst, local_ifrac)

    !===========================================================================
    ! 2-4. Master task: write sst_in.nc, signal Python, wait for result
    !===========================================================================

    if (my_task == master_task) then

       !--- delete stale done.flag if present ---
       inquire(file=DONE_FLAG, exist=flag_exists)
       if (flag_exists) call delete_file(DONE_FLAG)

       !--- write camulator_sst_in.nc ---
       call write_sst_nc(g_sst, g_ifrac, lsize_global, currentYMD, currentTOD, logunit)

       write(logunit, F00) 'CAMULATOR: wrote '//trim(SST_FILE)//', signalling Python...'
       call shr_sys_flush(logunit)

       !--- write go.flag (empty file) ---
       call write_flag(GO_FLAG)

       !--- poll for done.flag ---
       poll_iter = 0
       flag_exists = .false.
       do while (.not. flag_exists)
          call sleep(POLL_SLEEP_SEC)
          inquire(file=DONE_FLAG, exist=flag_exists)
          poll_iter = poll_iter + 1
          if (poll_iter > POLL_MAX_ITER) then
             call shr_sys_abort(trim(subname)// &
                  'CAMULATOR: timed out waiting for camulator_done.flag. '// &
                  'Is camulator_server.py running?')
          end if
          if (mod(poll_iter,60) == 0) then
             write(logunit,F01) 'CAMULATOR: still waiting for done.flag, iter=', poll_iter
             call shr_sys_flush(logunit)
          end if
       end do

       write(logunit,F01) 'CAMULATOR: received done.flag after ', poll_iter, ' iterations'
       call shr_sys_flush(logunit)

       !--- read camulator_cam_out.nc ---
       call read_cam_nc(g_u10, g_v10, g_tbot, g_zbot, g_qbot, g_pbot, &
                        g_fsds, g_flnsd, g_prect, lsize_global, logunit)

       !--- apply safety clamps ---
       where (abs(g_u10)  > WIND_MAX ) g_u10  = sign(WIND_MAX,  g_u10)
       where (abs(g_v10)  > WIND_MAX ) g_v10  = sign(WIND_MAX,  g_v10)
       where (g_qbot < QBOT_MIN      ) g_qbot  = QBOT_MIN
       where (g_flnsd < 0.0_R8       ) g_flnsd = 0.0_R8
       where (g_fsds  < 0.0_R8       ) g_fsds  = 0.0_R8
       where (g_prect < 0.0_R8       ) g_prect = 0.0_R8

       !--- log global diagnostics ---
       write(logunit,F02) 'CAMULATOR: global mean SST [K]   = ', sum(g_sst)/lsize_global
       write(logunit,F02) 'CAMULATOR: global mean |u_bot|m/s= ', &
            sum(sqrt(g_u10**2+g_v10**2))/lsize_global
       write(logunit,F02) 'CAMULATOR: global mean zbot [m]  = ', sum(g_zbot)/lsize_global
       write(logunit,F02) 'CAMULATOR: global mean FSDS W/m2 = ', sum(g_fsds)/lsize_global
       write(logunit,F02) 'CAMULATOR: global mean FLNSD W/m2= ', sum(g_flnsd)/lsize_global
       call shr_sys_flush(logunit)

       !--- remove done.flag ---
       call delete_file(DONE_FLAG)

    end if  ! master_task

    !===========================================================================
    ! 5. Broadcast global output arrays to all MPI ranks
    !===========================================================================

    call shr_mpi_bcast(g_u10,   mpicom, 'g_u10'  )
    call shr_mpi_bcast(g_v10,   mpicom, 'g_v10'  )
    call shr_mpi_bcast(g_tbot,  mpicom, 'g_tbot' )
    call shr_mpi_bcast(g_zbot,  mpicom, 'g_zbot' )
    call shr_mpi_bcast(g_qbot,  mpicom, 'g_qbot' )
    call shr_mpi_bcast(g_pbot,  mpicom, 'g_pbot' )
    call shr_mpi_bcast(g_fsds,  mpicom, 'g_fsds' )
    call shr_mpi_bcast(g_flnsd, mpicom, 'g_flnsd')
    call shr_mpi_bcast(g_prect, mpicom, 'g_prect')

    !===========================================================================
    ! 6. Scatter from global arrays to each rank's local a2x partition
    !===========================================================================

    allocate(local_u10  (lsize_local))
    allocate(local_v10  (lsize_local))
    allocate(local_tbot (lsize_local))
    allocate(local_zbot (lsize_local))
    allocate(local_qbot (lsize_local))
    allocate(local_pbot (lsize_local))
    allocate(local_fsds (lsize_local))
    allocate(local_flnsd(lsize_local))
    allocate(local_prect(lsize_local))

    call MPI_Scatterv(g_u10,  recvcounts, displs, MPI_DOUBLE_PRECISION, &
                      local_u10,   lsize_local,   MPI_DOUBLE_PRECISION, &
                      master_task, mpicom, ierr)
    call MPI_Scatterv(g_v10,  recvcounts, displs, MPI_DOUBLE_PRECISION, &
                      local_v10,   lsize_local,   MPI_DOUBLE_PRECISION, &
                      master_task, mpicom, ierr)
    call MPI_Scatterv(g_tbot, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                      local_tbot,  lsize_local,   MPI_DOUBLE_PRECISION, &
                      master_task, mpicom, ierr)
    call MPI_Scatterv(g_zbot, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                      local_zbot,  lsize_local,   MPI_DOUBLE_PRECISION, &
                      master_task, mpicom, ierr)
    call MPI_Scatterv(g_qbot, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                      local_qbot,  lsize_local,   MPI_DOUBLE_PRECISION, &
                      master_task, mpicom, ierr)
    call MPI_Scatterv(g_pbot, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                      local_pbot,  lsize_local,   MPI_DOUBLE_PRECISION, &
                      master_task, mpicom, ierr)
    call MPI_Scatterv(g_fsds, recvcounts, displs, MPI_DOUBLE_PRECISION, &
                      local_fsds,  lsize_local,   MPI_DOUBLE_PRECISION, &
                      master_task, mpicom, ierr)
    call MPI_Scatterv(g_flnsd,recvcounts, displs, MPI_DOUBLE_PRECISION, &
                      local_flnsd, lsize_local,   MPI_DOUBLE_PRECISION, &
                      master_task, mpicom, ierr)
    call MPI_Scatterv(g_prect,recvcounts, displs, MPI_DOUBLE_PRECISION, &
                      local_prect, lsize_local,   MPI_DOUBLE_PRECISION, &
                      master_task, mpicom, ierr)

    !===========================================================================
    ! 7. Populate a2x fields from scattered local arrays
    !    (CPL7 bulk formula computes wind stress from Sa_u/Sa_v + OCN SST)
    !===========================================================================

    lsize = mct_aVect_lsize(a2x)

    do n = 1, lsize

       !--- reference height [m] — dynamic, computed in camulator_server.py ---
       ! z_bot = (Rd/g) * (-ln(p_mid_frac)) * T_bot  where p_mid_frac = 0.5*(hybi[-2]+1)
       ! CAM6 L32: hybi[-2]=0.9851 → z_bot ≈ 0.2187*T_bot
       ! Range: ~50 m (polar, T=230 K) to ~67 m (tropical, T=305 K).
       ! shr_flux_atmOcn uses zbot in alz=log(zbot/zref) for Monin-Obukhov corrections.
       a2x%rAttr(kz,    n) = local_zbot(n)

       !--- winds [m s-1] ---
       a2x%rAttr(ku,    n) = local_u10(n)
       a2x%rAttr(kv,    n) = local_v10(n)

       !--- temperature [K] ---
       a2x%rAttr(ktbot, n) = local_tbot(n)
       a2x%rAttr(kptem, n) = local_tbot(n)   ! potential temp ≈ tbot at surface

       !--- specific humidity [kg kg-1] ---
       a2x%rAttr(kshum, n) = local_qbot(n)

       !--- pressure [Pa] ---
       a2x%rAttr(kpbot, n) = local_pbot(n)
       a2x%rAttr(kpslv, n) = local_pbot(n)

       !--- air density from ideal gas law: rho = p / (Rd * Tv)
       !    Tv = T * (1 + 0.608*q) (virtual temperature)
       dens = local_pbot(n) / (287.04_R8 * local_tbot(n) * (1.0_R8 + 0.608_R8*local_qbot(n)))
       a2x%rAttr(kdens, n) = dens

       !--- downward longwave [W m-2] ---
       a2x%rAttr(klwdn, n) = local_flnsd(n)

       !--- shortwave: split downwelling SW (FSDS) into four bands ---
       !    FSDS is downwelling before surface reflection (reconstructed in Python).
       !    CPL7 seq_flux_mct.F90 applies ocean/ice albedo to these fields
       !    to compute reflected SW, so we MUST pass downwelling (not net) here.
       swdn = local_fsds(n)
       a2x%rAttr(kswvdr, n) = swdn * sw_frac_swvdr
       a2x%rAttr(kswndr, n) = swdn * sw_frac_swndr
       a2x%rAttr(kswvdf, n) = swdn * sw_frac_swvdf
       a2x%rAttr(kswndf, n) = swdn * sw_frac_swndf
       a2x%rAttr(kswnet, n) = swdn

       !--- precipitation [kg m-2 s-1] from PRECT [m s-1] * rho_w
       !    rain if tbot > freezing, snow otherwise
       rain_kg = local_prect(n) * rho_w
       snow_kg = 0.0_R8
       if (local_tbot(n) < tKFrz) then
          snow_kg = rain_kg
          rain_kg = 0.0_R8
       end if
       a2x%rAttr(krl, n) = rain_kg   ! large-scale liquid
       a2x%rAttr(krc, n) = 0.0_R8   ! convective liquid (set to 0 for now)
       a2x%rAttr(ksl, n) = snow_kg   ! large-scale snow
       a2x%rAttr(ksc, n) = 0.0_R8   ! convective snow

    end do  ! lsize

    deallocate(local_u10, local_v10, local_tbot, local_zbot, local_qbot, local_pbot)
    deallocate(local_fsds, local_flnsd, local_prect)

  end subroutine datm_datamode_camulator_run


  !============================================================================
  ! Private helper: write camulator_sst_in.nc
  !============================================================================

  subroutine write_sst_nc(sst, ifrac, ngrid, ymd, tod, logunit)

    real(R8),    intent(in) :: sst(ngrid), ifrac(ngrid)
    integer(IN), intent(in) :: ngrid, ymd, tod, logunit

    integer(IN) :: ncid, dimid, vid_sst, vid_ifrac, vid_ymd, vid_tod, ierr

    character(*), parameter :: F00 = "('(datm_cam write_sst_nc) ',a)"

    ierr = nf90_create(SST_FILE, NF90_CLOBBER, ncid)
    if (ierr /= NF90_NOERR) call shr_sys_abort('CAMULATOR: cannot create '//SST_FILE)

    ierr = nf90_def_dim(ncid, 'ngrid', ngrid, dimid)
    ierr = nf90_def_var(ncid, 'sst',   NF90_DOUBLE, [dimid], vid_sst  )
    ierr = nf90_def_var(ncid, 'ifrac', NF90_DOUBLE, [dimid], vid_ifrac)
    ierr = nf90_def_var(ncid, 'ymd',   NF90_INT,    vid_ymd )
    ierr = nf90_def_var(ncid, 'tod',   NF90_INT,    vid_tod )
    ierr = nf90_put_att(ncid, vid_sst,   'units', 'K'  )
    ierr = nf90_put_att(ncid, vid_sst,   'long_name', 'ocean surface temperature from CPL')
    ierr = nf90_put_att(ncid, vid_ifrac, 'units', '1'  )
    ierr = nf90_put_att(ncid, vid_ifrac, 'long_name', 'sea ice fraction from CPL')
    ierr = nf90_put_att(ncid, vid_ymd,   'long_name', 'date YYYYMMDD')
    ierr = nf90_put_att(ncid, vid_tod,   'long_name', 'time of day seconds')
    ierr = nf90_enddef(ncid)

    ierr = nf90_put_var(ncid, vid_sst,   sst  )
    ierr = nf90_put_var(ncid, vid_ifrac, ifrac)
    ierr = nf90_put_var(ncid, vid_ymd,   ymd  )
    ierr = nf90_put_var(ncid, vid_tod,   tod  )
    ierr = nf90_close(ncid)

    if (ierr /= NF90_NOERR) &
         call shr_sys_abort('CAMULATOR: error closing '//SST_FILE)

  end subroutine write_sst_nc


  !============================================================================
  ! Private helper: read camulator_cam_out.nc
  !============================================================================

  subroutine read_cam_nc(u10, v10, tbot, zbot, qbot, pbot, fsds, flnsd, prect, &
                          ngrid, logunit)

    real(R8),    intent(out) :: u10(ngrid), v10(ngrid), tbot(ngrid), zbot(ngrid)
    real(R8),    intent(out) :: qbot(ngrid), pbot(ngrid)
    real(R8),    intent(out) :: fsds(ngrid), flnsd(ngrid), prect(ngrid)
    integer(IN), intent(in)  :: ngrid, logunit

    integer(IN) :: ncid, varid, ierr

    ierr = nf90_open(CAM_FILE, NF90_NOWRITE, ncid)
    if (ierr /= NF90_NOERR) &
         call shr_sys_abort('CAMULATOR: cannot open '//CAM_FILE// &
                            '. Is camulator_server.py producing output?')

    call nc_get(ncid, 'u10',   u10,   ngrid)
    call nc_get(ncid, 'v10',   v10,   ngrid)
    call nc_get(ncid, 'tbot',  tbot,  ngrid)
    call nc_get(ncid, 'zbot',  zbot,  ngrid)
    call nc_get(ncid, 'qbot',  qbot,  ngrid)
    call nc_get(ncid, 'pbot',  pbot,  ngrid)
    call nc_get(ncid, 'fsds',  fsds,  ngrid)
    call nc_get(ncid, 'flnsd', flnsd, ngrid)
    call nc_get(ncid, 'prect', prect, ngrid)

    ierr = nf90_close(ncid)

  end subroutine read_cam_nc


  !============================================================================
  ! Private helper: read one double variable from an open NetCDF file
  !============================================================================

  subroutine nc_get(ncid, varname, arr, n)
    integer(IN),      intent(in)  :: ncid, n
    character(len=*), intent(in)  :: varname
    real(R8),         intent(out) :: arr(n)
    integer(IN) :: varid, ierr
    ierr = nf90_inq_varid(ncid, varname, varid)
    if (ierr /= NF90_NOERR) &
         call shr_sys_abort('CAMULATOR: variable '//trim(varname)//' not in '//CAM_FILE)
    ierr = nf90_get_var(ncid, varid, arr)
    if (ierr /= NF90_NOERR) &
         call shr_sys_abort('CAMULATOR: error reading '//trim(varname)//' from '//CAM_FILE)
  end subroutine nc_get


  !============================================================================
  ! Private helper: write an empty flag file
  !============================================================================

  subroutine write_flag(filename)
    !--- Use shr_file_getUnit so we never collide with CIME's logunit (which can
    !    be unit 99 in this config).  Using a hardcoded unit=99 was accidentally
    !    disconnecting the coupler's logunit from cpl.log, routing all coupler
    !    output (tStamp, SYPD, SUCCESSFUL TERMINATION) to fort.99 instead.
    character(len=*), intent(in) :: filename
    integer(IN) :: nu
    nu = shr_file_getUnit()
    open(unit=nu, file=filename, status='replace', action='write')
    close(nu)
    call shr_file_freeUnit(nu)
  end subroutine write_flag


  !============================================================================
  ! Private helper: delete a file (Fortran 2003 open+close trick)
  !============================================================================

  subroutine delete_file(filename)
    character(len=*), intent(in) :: filename
    integer(IN) :: nu
    nu = shr_file_getUnit()
    open(unit=nu, file=filename, status='old')
    close(nu, status='delete')
    call shr_file_freeUnit(nu)
  end subroutine delete_file


end module datm_datamode_camulator_mod
