!!****m* ABINIT/m_indabs
!! NAME
!!  m_indabs
!!
!! FUNCTION
!!  Compute indirect absorption rates
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (NC)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_indabs

 use defs_basis
 use iso_c_binding
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hide_blas
 use m_copy
 use m_ifc
 use m_ebands
 use m_wfk
 use m_ddb
 use m_ddk
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_wfd
 use m_skw
 use m_krank
 use m_lgroup
 use m_ephwg
 use m_sort
 use m_hdr
 use m_sigtk
 use m_ephtk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_nctk
 use m_rf2
 use m_dtset
 use m_dtfil
 use m_clib
 use m_mkffnl

 use defs_abitypes,    only : mpi_type
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_time,           only : cwtime, cwtime_report, timab, sec2str
 use m_bz_mesh,        only : findqg0
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, c2r, get_diag, linfit, iseven, simpson_cplx, simpson, print_arr
 use m_io_tools,       only : iomode_from_fname, file_exists, is_open, open_file
 use m_special_funcs,  only : gaussian
 use m_fftcore,        only : ngfft_seq, sphereboundary, get_kg, kgindex
 use m_cgtk,           only : cgtk_rotate
 use m_cgtools,        only : cg_zdotc, cg_real_zdotc, cg_zgemm, fxphas_seq
 use m_crystal,        only : crystal_t
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt
 use m_occ,            only : occ_fd, occ_be !occ_dfde,
 use m_kg,             only : getph, mkkpg
 use m_bz_mesh,        only : isamek
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_ioarr,          only : read_rhor
 use m_pawang,         only : pawang_type, make_angular_mesh
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawrhoij,       only : pawrhoij_type
 use m_pawfgr,         only : pawfgr_type
 use m_symtk,          only : littlegroup_q
 use m_dfpt_cgwf,      only : dfpt_cgwf
 use m_dynmat,         only : pheigvec_normalize, massmult_and_breaksym, phdispl_from_eigvec
 use m_geometry,       only : phdispl_cart2red

 implicit none

public :: indabs

contains


subroutine indabs(wfk0_path,dtfil,ngfft,ngfftf,dtset,cryst,ebands,dvdb,ifc,&
                      pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path 
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dvdb_t),target,intent(inout) :: dvdb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:) !natom,nspinor*usecprj)

 ! local variables
 character(len=500) :: msg
 integer :: my_rank,nprocs,nfft,nfftf,mgfft,mgfftf,mband,mpw
 integer :: natom,natom3,nsppol,nspinor,nspden,nkpt,n1,n2,n3,n4,n5,n6
 integer :: spin,ik,iq,sband,fband,nkibz,nkbz
 integer :: usevnl,optlocal,optnl,opt_gvnlx1,usecprj,sij_opt 
 integer :: interpolated,comm_rpt,cplex,db_iqpt,timrev_q
 integer :: ikq,nkpg,nkpg1,ik_ibz,isym_k,trev_k,isym_kq,ikq_ibz,trev_kq
 integer :: istwf_kq,npw_kq,istwf_kqibz,npw_kqibz,istwf_k,npw_k,istwf_kibz,npw_kibz
 integer :: ierr,iband,imyp,idir,ipert,my_npert
 integer :: max_umklapp,i1,i2,i3,ipw,onpw,ii
 logical :: gen_eigenpb,use_ftinterp,verbose,isirr_k,isirr_kq
 integer,parameter :: tim_getgh1c=1,berryopt0=0,istw1=1,ider0=0,idir0=0
 integer,parameter :: useylmgr1=0,useylmgr0=0,ndat1=1,istwfk1=1
 type(wfd_t) :: wfd
 type(ddkop_t) :: ddkop
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 real(dp) :: cpu_all,wall_all,gflops_all
 real(dp) :: ecut 

 ! arrays
 real(dp) :: v_bks(2,3),kpt(3),kpt_ibz(3),qpt(3),kqpt(3),kptu(3)
 real(dp) :: ylmgr_dum(1,1,1)
 integer :: symq(4,2,cryst%nsym), g0_k(3),g0_kq(3),sym_k(6),sym_kq(6)
 integer :: work_ngfft(18),gmax(3)

 integer,allocatable :: nband(:,:),wfd_istwfk(:),bz2ibz(:,:),gtmp(:,:),kg_kq(:,:),kg_k(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 real(dp),allocatable :: ph1d(:,:),cg_ket(:,:),cg_bra(:,:)
 real(dp),allocatable :: wtk(:), kibz(:,:), kbz(:,:)
 real(dp),allocatable :: grad_berry(:,:),vtrial(:,:),vlocal(:,:,:,:),h1kets_kq(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:),ffnlk(:,:,:,:),kpg1_k(:,:),kpg_k(:,:),ffnl1(:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:),ph3d(:,:,:),ph3d1(:,:,:) 
 real(dp),allocatable :: kinpw1(:),dkinpw(:),cgwork(:,:)
 complex(dpc),allocatable :: vmat(:,:,:,:,:)

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor
 nspden = dtset%nspden; nkpt = ebands%nkpt

 ! FFT meshes from input file (not necessarly equal to the ones found in the external files.
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)
 call wrtout(std_out, sjoin("FFT grid:", ltoa(ngfft)))
 

 ! For now no parallelism over perturbations. Each processor handles all perturbations
 my_npert = natom3

 ! Construct object to store final results.
 ecut = dtset%ecut ! dtset%dilatmx

 mband = dtset%mband

 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask, (mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur, (mband, nkpt, nsppol))
 keep_ur = .False.

 bks_mask = .True.

 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1
 nband = mband

 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, mband, nband, nkpt, nsppol, bks_mask,&
              dtset%nspden, nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
              dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(nband)


 mpw = maxval(wfd%npwarr)
 call wrtout(std_out,sjoin("mpw:", itoa(mpw)))
 call wrtout(std_out,sjoin("Npw for all k points",ltoa(wfd%npwarr)))
 call wrtout(std_out,sjoin("istwfk for all k points",ltoa(wfd%istwfk)))
 print *, "gmet", cryst%gmet
 ABI_MALLOC(vmat, (3,mband,mband,nkpt,nsppol))
 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 ABI_CALLOC(cg_ket, (2, mpw*ebands%nspinor))
 ABI_CALLOC(cg_bra, (2, mpw*ebands%nspinor))


 do spin = 1,nsppol   
   do ik = 1,nkpt
     npw_k = wfd%npwarr(ik); istwf_k = wfd%istwfk(ik)
     ! ik loop outside so this is only called once per k point
     call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, ebands%kptns(:,ik), istwf_k, npw_k, wfd%kdata(ik)%kg_k)
     ! Precompute all needed matrix elements in IBZ
     do sband = 1,mband
       call wfd%copy_cg(sband,ik,spin,cg_ket)
       call ddkop%apply(ebands%eig(sband,ik,spin), npw_k, nspinor, cg_ket, cwaveprj0)
       do fband = sband,mband       
         ! Don't bother computing if they have the same population
         if (abs(ebands%occ(sband,ik,spin) - ebands%occ(fband,ik,spin)) .lt. 1.0d-12) cycle
         call wfd%copy_cg(fband,ik,spin,cg_bra)
         v_bks = ddkop%get_braket(ebands%eig(sband,ik,spin),istwf_k, npw_k, nspinor, cg_bra, mode="cart")
         vmat(:,sband,fband,ik,spin) = cmplx(v_bks(1,:), v_bks(2,:),kind=dp)
         vmat(:,fband,sband,ik,spin) = cmplx(v_bks(1,:), -v_bks(2,:),kind=dp)
       end do !fband
     end do !sband
   end do !ik
 end do !spin
 

 ! if PAW, one has to solve a generalized eigenproblem
 ! BE careful here because I will need sij_opt == -1
 usecprj = 0
 gen_eigenpb = psps%usepaw == 1; sij_opt = 0; if (gen_eigenpb) sij_opt = 1

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1   ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2      ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0 ! gvnlx1 is output

 ABI_MALLOC(grad_berry, (2, nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 ! 1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 ! 2) Perform the setup needed for the non-local factors:
 !
 ! Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 ! PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 ! Get one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx, natom, n1, n2, n3, ph1d, cryst%xred)

 call init_hamiltonian(gs_hamkq, psps, pawtab, nspinor, nsppol, nspden, natom,&
  dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg,&
  comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab,&
  usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, use_gpu_cuda=dtset%use_gpu_cuda)

 ! Allocate work space arrays.
 ! vtrial and vlocal are required for Sternheimer (H0)
 ! Note nvloc in vlocal (we will select one/four spin components afterwards)
 ABI_CALLOC(vtrial, (nfftf, nspden))
 ABI_CALLOC(vlocal, (n4, n5, n6, gs_hamkq%nvloc))

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)
 
 !Get mapping between kpts in ibz and bz
 call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, ebands%kptopt, ebands%nshiftk, ebands%shiftk, &
   nkibz, kibz, wtk, nkbz, kbz, bz2ibz=bz2ibz) 

 ABI_FREE(wtk)
 ABI_FREE(kibz)
 ABI_CHECK(nkibz == ebands%nkpt, "nkibz != ebands%nkpt")
 verbose = .true.
 if (nkpt > 50) verbose = .false.
 ! Initialize v1scf to store first order potential
 cplex = 2
 ABI_MALLOC(v1scf, (cplex, nfftf, nspden, dvdb%my_npert))
 call wrtout(std_out, sjoin("MSG: nkpt", itoa(nkbz)))


 ABI_MALLOC(kg_k, (3,mpw))
 ABI_MALLOC(kg_kq, (3,mpw))
 ! Spherical Harmonics for useylm == 1.
 ABI_MALLOC(ylm_k, (mpw, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylm_kq, (mpw, psps%mpsang**2 * psps%useylm))
 ! TODO: useylmgr1 has value 0 but the name implies it should be 1
 ABI_MALLOC(ylmgr_kq, (mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))

 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 max_umklapp = 10 
 
 do ik=1,nkbz
   kpt = kbz(:, ik)
   do i3=-max_umklapp,max_umklapp
     do i2=-max_umklapp,max_umklapp
       do i1=-max_umklapp,max_umklapp
         kptu = kpt + one*[i1, i2, i3]
         ! TODO: g0 umklapp here can enter into play gmax may not be large enough!
         call get_kg(kptu, istwfk1, 1.1_dp * ecut, cryst%gmet, onpw, gtmp)
         mpw = max(mpw, onpw)
         call wrtout(std_out, sjoin("npw trial",itoa(onpw)))
         do ipw=1,onpw
           do ii=1,3
             gmax(ii) = max(gmax(ii), abs(gtmp(ii, ipw)))
           end do
         end do
         ABI_FREE(gtmp)
       end do
     end do
   end do
 end do
 call wrtout(std_out,sjoin("New mpw:", itoa(mpw)))

!! Init work_ngfft
!call ngfft_seq(work_ngfft, gmax)
!if (verbose) call wrtout(std_out,sjoin("work_ngfft:", ltoa(work_ngfft)))
 
 ! Loop over q in full BZ
 do iq=1,nkbz
   qpt = kbz(:,iq)
   if (verbose) call wrtout(std_out, sjoin("MSG: Treating qpt", ktoa(qpt)))

   ! Interpolation setup taken from m_gkk
   interpolated = 0
   if (dtset%eph_use_ftinterp /= 0) then
     !MSG_WARNING(sjoin("Enforcing FT interpolation for q-point", ktoa(qpt)))
     comm_rpt = xmpi_comm_self
     call dvdb%ftinterp_setup(dtset%ddb_ngqpt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)
     call dvdb%ftinterp_qpt(qpt, nfftf, ngfftf, v1scf, dvdb%comm_rpt)
     interpolated = 1
   else
     ! Find the index of the q-point in the DVDB.
     db_iqpt = dvdb%findq(qpt)
     if (db_iqpt /= -1) then
       !if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Found: ",ktoa(qpt)," in DVDB with index ",itoa(db_iqpt)))
       ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
       ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
       call dvdb%readsym_allv1(db_iqpt, cplex, nfftf, ngfftf, v1scf, comm)
     else
       MSG_WARNING(sjoin("Cannot find q-point:", ktoa(qpt), "in DVDB file"))
     end if
   end if

   do ik=1,nkbz
     do spin = 1,nsppol
       kpt = kbz(:,ik)
       kqpt = kpt + qpt
       ! Find the index of the k+q point
       call findqg0(ikq, g0_kq, kqpt, nkbz, kbz, [1,1,1])
       if (verbose) call wrtout(std_out, sjoin("MSG: Treating transitions from kpt", ktoa(kpt), "to k+q", ktoa(kqpt), "connected by qpt", ktoa(qpt)))

       ! Get the symmetry information mapping the k point to the IBZ
       sym_k = bz2ibz(:,ik)
       ik_ibz = sym_k(1); isym_k = sym_k(2)
       trev_k = sym_k(6); g0_k = sym_k(3:5)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       kpt_ibz = ebands%kptns(:,ik_ibz)
       if (verbose) then
         call wrtout(std_out, sjoin("Symmetry information for kpt", ktoa(kpt), ":"))
         call wrtout(std_out, sjoin("Index/kpt in ibz:", itoa(ik_ibz), ktoa(ebands%kptns(:,ik_ibz))))
         call wrtout(std_out, sjoin("Connected by symmetry operation", itoa(isym_k), itoa(trev_k)))
         call wrtout(std_out, sjoin("G0", ltoa(g0_k)))
       end if
   
       ! Get the symmetry information mapping the k point to the IBZ
       sym_kq = bz2ibz(:,ikq)
       ikq_ibz = sym_kq(1); isym_kq = sym_kq(2)
       trev_kq = sym_kq(6); g0_kq = sym_kq(3:5)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       if (verbose) then
         call wrtout(std_out, sjoin("Symmetry information for kpt", ktoa(kqpt), ":"))
         call wrtout(std_out, sjoin("Index/kpt in ibz:", itoa(ikq_ibz), ktoa(ebands%kptns(:,ikq_ibz))))
         call wrtout(std_out, sjoin("Connected by symmetry operation", itoa(isym_kq), itoa(trev_kq)))
         call wrtout(std_out, sjoin("G0", ltoa(g0_kq)))
       end if

       istwf_kibz = wfd%istwfk(ik_ibz); npw_kibz = wfd%npwarr(ik_ibz)
       ! Get npw_k, kg_k for k.
       if (isirr_k) then
         ! Copy u_k(G)
         istwf_k = istwf_kibz; npw_k = npw_kibz
         ABI_CHECK(mpw >= npw_k, "mpw < npw_k")
         kg_k(:,1:npw_k) = wfd%kdata(ik_ibz)%kg_k
       else
         ! Reconstruct u_k(G) from the IBZ image.
         istwf_k = 1
         call get_kg(kpt, istwf_k, ecut, cryst%gmet, npw_k, gtmp)
         ABI_CHECK(mpw >= npw_k, "mpw < npw_k")
         kg_k(:,1:npw_k) = gtmp(:,:npw_k)
         ABI_FREE(gtmp)
       end if


       istwf_kqibz = wfd%istwfk(ikq_ibz); npw_kqibz = wfd%npwarr(ikq_ibz)
       ! Get npw_kq, kg_kq for k+q.
       if (isirr_kq) then
         ! Copy u_kq(G)
         istwf_k = istwf_kqibz; npw_kq = npw_kqibz
         ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
         kg_kq(:,1:npw_kq) = wfd%kdata(ikq_ibz)%kg_k
       else
         ! Reconstruct u_kq(G) from the IBZ image.
         istwf_kq = 1
         call get_kg(kqpt, istwf_kq, ecut, cryst%gmet, npw_kq, gtmp)
         ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
         kg_kq(:,1:npw_kq) = gtmp(:,:npw_kq)
         ABI_FREE(gtmp)
       end if

       ! Allocate cgwork to read in cg at irreducible k
       ABI_MALLOC(cgwork, (2, npw_k*nspinor))


       ABI_MALLOC_OR_DIE(vlocal1, (cplex*n4, n5, n6, gs_hamkq%nvloc, my_npert), ierr)     

       ! Compute k+G vectors
       nkpg = 3*dtset%nloalg(3)
       ABI_MALLOC(kpg_k, (npw_k, nkpg))
       if (nkpg > 0) call mkkpg(kg_k, kpg_k, kpt, nkpg, npw_k)

       !TODO: This loop recomputes ffnls many times. No effect for my pseudopotentials but I'm not sure about more sophisticated ones
       ! Compute nonlocal form factors ffnlk at (k+G)
       ABI_MALLOC(ffnlk, (npw_k, 1, psps%lmnmax, psps%ntypat))
       call mkffnl(psps%dimekb, 1, psps%ekb, ffnlk, psps%ffspl, &
         cryst%gmet, cryst%gprimd, ider0, idir0, psps%indlmn, kg_k, kpg_k, kpt, psps%lmnmax, &
         psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg, npw_k, psps%ntypat, &
         psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_k, ylmgr_dum, &
         comm=comm)
       
       ! Compute k+q+G vectors
       nkpg1 = 3*dtset%nloalg(3) ! Determine size for kpg1_k. Either 3 (compute vectors) or 0 (don't compute, don't allocate) 
       ABI_MALLOC(kpg1_k, (npw_kq, nkpg1)) !Initialize array to store k+G vectors
       if (nkpg1 > 0) call mkkpg(kg_kq, kpg1_k, kqpt, nkpg1, npw_kq)

       ! Compute nonlocal form factors ffnl1 at (k+q+G) 
       ! TODO: comm? ffnl_request?
       ABI_MALLOC(ffnl1, (npw_kq, 1, psps%lmnmax, psps%ntypat))
       call mkffnl(psps%dimekb, 1, psps%ekb, ffnl1, psps%ffspl, cryst%gmet, cryst%gprimd, ider0, idir0, &
         psps%indlmn, kg_kq, kpg1_k, kqpt, psps%lmnmax, psps%lnmax, psps%mpsang, psps%mqgrid_ff, nkpg1, &
         npw_kq, psps%ntypat, psps%pspso, psps%qgrid_ff, cryst%rmet, psps%usepaw, psps%useylm, ylm_kq, ylmgr_kq, &
         comm=comm)

       do iband=1,mband
         ! Symmetrize k wavefunctions in the BZ from IBZ (if needed).
         if (isirr_k) then
           if (verbose) call wrtout(std_out, sjoin("Copying cg for irreducible k point",itoa(ik_ibz)))
           ! Copy u_kq(G)
           !call wfd%copy_cg(iband, ik_ibz, spin, ket_k)
         else
           ! Reconstruct u_kq(G) from the IBZ image.
           ! Use cgwork as workspace array, results stored in 
           if (verbose) call wrtout(std_out, sjoin("band", itoa(iband), "ik_ikbz", itoa(ik_ibz), "spin", itoa(spin)))
           call wfd%copy_cg(iband, ik_ibz, spin, cgwork)
          !call cgtk_rotate(cryst, kpt_ibz, isym_k, trev_k, g0_k, nspinor, ndat1, &
          !                 npw_kibz, wfd%kdata(ik_ibz)%kg_k, &
          !                 npw_k, kg_k, istwf_kibz, istwf_k, cgwork, ket_k, work_ngfft, work)
         end if
         
       end do !iband

       ! Loop over all 3*natom perturbations (Each CPU prepares its own potentials)
       ! In the inner loop, I calculate H1 * psi_k, stored in h1kets_kq on the k+q sphere.
       do imyp=1,my_npert
         ! From m_gkk: assign each perturbation a direction and perturbation number
         idir = mod(imyp-1, 3) + 1
         ipert = (imyp - idir) / 3 + 1
         if (verbose) write(msg, '(a,2i4)') " Treating ipert, idir = ", ipert, idir
         call wrtout(std_out, msg, do_flush=.True.)


         ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
         call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_hamkq%nvloc, &
           pawfgr, mpi_enreg, vtrial, v1scf(:,:,:,imyp), vlocal, vlocal1(:,:,:,:,imyp))

         ! Continue to initialize the Hamiltonian (call it here to support dfpt_cgwf Sternheimer).
         call gs_hamkq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)
         call rf_hamkq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,imyp), with_nonlocal=.true.)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         ! Allocates ph3d, kinpw1, dkinpw (not ph3d1?)
         call getgh1c_setup(gs_hamkq, rf_hamkq, dtset, psps, kpt, kqpt, idir, ipert, &  ! In
           cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, &             ! In
           npw_k, npw_kq, useylmgr1, kg_k, ylm_k, kg_kq, ylm_kq, ylmgr_kq, &         ! In
           dkinpw, nkpg, nkpg1, kpg_k, kpg1_k, kinpw1, ffnlk, ffnl1, ph3d, ph3d1, &  ! Out
           reuse_kpg_k=1, reuse_kpg1_k=1, reuse_ffnlk=1, reuse_ffnl1=1)              ! Reuse some arrays

       ABI_FREE(ph3d)
       ABI_SFREE(ph3d1)
       ABI_FREE(kinpw1)
       ABI_FREE(dkinpw)
       
       end do !imyp
   
       ABI_FREE(ffnl1)
       ABI_FREE(ffnlk)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)
       ABI_FREE(vlocal1)
       ABI_FREE(cgwork)
     end do !spin
   end do !ik
 end do !iq

 ABI_FREE(kg_kq)
 ABI_FREE(kg_k)
 ABI_FREE(vmat)
 ABI_FREE(cg_ket)
 ABI_FREE(cg_bra)
 ABI_FREE(grad_berry)
 ABI_FREE(vtrial)
 ABI_FREE(vlocal)
 ABI_FREE(v1scf)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr_kq)

 end subroutine indabs

end module m_indabs
