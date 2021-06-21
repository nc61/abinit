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
 integer :: my_rank,nprocs,nfft,nfftf,mgfft,mgfftf,mband,mpw
 integer :: natom,natom3,nsppol,nspinor,nspden,nkpt,n1,n2,n3,n4,n5,n6
 integer :: isppol,ik,iq,sband,fband,npw_k,istwf_k,nkibz,nkbz
 integer :: usevnl,optlocal,optnl,opt_gvnlx1,usecprj,sij_opt 
 logical :: gen_eigenpb 
 integer,parameter :: berryopt0=0
 type(wfd_t) :: wfd
 type(ddkop_t) :: ddkop
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 real(dp) :: cpu_all,wall_all,gflops_all
 real(dp) :: ecut 

 ! arrays
 real(dp) :: v_bks(2,3), kk(3)

 integer,allocatable :: nband(:,:),wfd_istwfk(:),bz2ibz(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 real(dp),allocatable :: ph1d(:,:),cg_ket(:,:),cg_bra(:,:)
 real(dp),allocatable :: wtk(:), kibz(:,:), kbz(:,:)
 real(dp),allocatable :: grad_berry(:,:),vtrial(:,:),vlocal(:,:,:,:)
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
 ABI_MALLOC(vmat, (3,mband,mband,nkpt,nsppol))
 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 ABI_CALLOC(cg_ket, (2, mpw*ebands%nspinor))
 ABI_CALLOC(cg_bra, (2, mpw*ebands%nspinor))


 do isppol = 1,nsppol   
   do ik = 1,nkpt
     npw_k = wfd%npwarr(ik); istwf_k = wfd%istwfk(ik)
     ! ik loop outside so this is only called once per k point
     call ddkop%setup_spin_kpoint(dtset, cryst, psps, isppol, ebands%kptns(:,ik), istwf_k, npw_k, wfd%kdata(ik)%kg_k)
     ! Precompute all needed matrix elements in IBZ
     do sband = 1,mband
       call wfd%copy_cg(sband,ik,isppol,cg_ket)
       call ddkop%apply(ebands%eig(sband,ik,isppol), npw_k, nspinor, cg_ket, cwaveprj0)
       do fband = sband,mband       
         ! Don't bother computing if they have the same population
         if (abs(ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)) .lt. 1.0d-12) cycle
         call wfd%copy_cg(fband,ik,isppol,cg_bra)
         v_bks = ddkop%get_braket(ebands%eig(sband,ik,isppol),istwf_k, npw_k, nspinor, cg_bra, mode="cart")
         vmat(:,sband,fband,ik,isppol) = cmplx(v_bks(1,:), v_bks(2,:),kind=dp)
         vmat(:,fband,sband,ik,isppol) = cmplx(v_bks(1,:), -v_bks(2,:),kind=dp)
       end do !fband
     end do !sband
   end do !ik
 end do !isppol
 

 ABI_FREE(wtk)
 ABI_FREE(kibz)
 ABI_CHECK(nkibz == ebands%nkpt, "nkibz != ebands%nkpt")

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
 ! vtrial and vlocal are required for Sternheimer (H0). DFPT routines do not need it.
 ! Note nvloc in vlocal (we will select one/four spin components afterwards)
 ABI_CALLOC(vtrial, (nfftf, nspden))
 ABI_CALLOC(vlocal, (n4, n5, n6, gs_hamkq%nvloc))

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)
 
 !Get mapping between kpts in ibz and bz
 call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, ebands%kptopt, ebands%nshiftk, ebands%shiftk, &
   nkibz, kibz, wtk, nkbz, kbz, bz2ibz=bz2ibz) 

 ! Loop over q in full BZ
 do iq=1,nkbz
      
 end do !iq

 ABI_FREE(vmat)
 ABI_FREE(cg_ket)
 ABI_FREE(cg_bra)
 ABI_FREE(grad_berry)
 ABI_FREE(vtrial)
 ABI_FREE(vlocal)

 end subroutine indabs

end module m_indabs
