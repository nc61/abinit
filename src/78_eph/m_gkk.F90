!!****m* ABINIT/m_gkk
!! NAME
!!
!! FUNCTION
!!  Tools for the computation of electron-phonon coupling matrix elements (gkk)
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (GKA, MG)
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

module m_gkk

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_dtset
 use m_ifc
 use m_ebands
 use m_ddb
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_krank
 use m_pawcprj
 use m_wfk
 use m_nctk
 use m_dtfil
 use m_dvq
 use m_sort
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_abitypes,    only : MPI_type
 use m_time,           only : cwtime, sec2str
 use m_io_tools,       only : iomode_from_fname
 use m_fstrings,       only : itoa, sjoin, ktoa, ltoa, strcat
 use m_symtk,          only : littlegroup_q
 use m_fftcore,        only : get_kg
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : findqg0
 use m_cgtools,        only : dotprod_g
 use m_kg,             only : getph
 use m_numeric_tools,  only : arth 
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
 use m_eig2d,          only : gkk_t, gkk_init, gkk_ncwrite, gkk_free
 use m_wfd,            only : wfd_init, wfd_t
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_symtk,          only : matr3inv, mati3inv

 !My added libs
 use m_hdr
 use m_wfk
 use m_io_tools,      only: get_unit
 use m_fstrings,      only: endswith
 use m_optic_tools,   only: sym2cart, pmat2cart
 use m_htetra
 use m_kpts,          only: tetra_from_kptrlatt
 use m_numeric_tools, only: arth

 implicit none

 private
!!***

 public :: eph_gkk
 public :: absrate_ind2
 public :: ncwrite_v1qnu          ! Compute \delta V_{q,nu)(r) and dump results to netcdf file.
 public :: gkq_atm_to_gkq_nu

contains  !===========================================================================

subroutine absrate_ind2(wfk0_path,dtfil,ngfft,ngfftf,dtset,cryst,ebands,dvdb,ifc,&
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

!Local variables ------------------------------
!scalars
 type(dvqop_t) :: dvqop
 type(krank_t) :: krank
 character(len=500) :: errmsg
 character(len=fnlen) :: ddkfile_1, ddkfile_2, ddkfile_3, gs_wfkpath
 complex(dpc) :: self_energy
 integer :: npw_k, istwfk,istwfkq, usecprj
 integer :: branch, ib1, ib2, ik, iq, spin
 integer :: sband, fband, iband
 integer :: mpw, my_start, my_stop, num_fband=0, num_sband=0, nw
 integer :: ierr, isym, sym, krank_prev=-2
 integer :: rank(cryst%nsym), iperm(cryst%nsym)
 real(dp) :: wmin, wmax, step 
 type(htetra_t) :: htetra
 type(wfd_t) :: wfd

 !arrays
 complex(dpc) :: pmat_1(3), pmat_2(3), pmat_mag2
 complex(dpc),allocatable :: pmat(:,:,:,:,:), numerator(:), detuning_plus(:), detuning_minus(:)
 complex(dpc),allocatable :: path_summand_plus(:), path_summand_1_plus(:), path_summand_2_plus(:), path_sum_plus(:)
 complex(dpc),allocatable :: path_summand_minus(:), path_summand_1_minus(:), path_summand_2_minus(:), path_sum_minus(:)
 integer :: g0_k(3)
 integer,allocatable :: wfd_istwfk(:),nband(:,:), kg_k(:,:), bz2ibz_indexes(:), fbands(:), sbands(:), ikq(:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 real(dp) :: qpt(3),kpt(3),kqpt(3),klatt(3,3),rlatt(3,3),symrec(3,3), kpt_new(3)
 real(dp),allocatable :: energy_fs(:) 
 real(dp),allocatable :: cg_sband(:,:), cg_iband(:,:), cg_fband(:,:), gkq(:,:), wmesh(:),weights(:,:,:,:), transrate_integral(:), transrate_total(:), phonon_populations(:)
 real(dp),allocatable :: phonon_energies(:,:)

 ! Load all optical matrix elements
 ddkfile_1 = "AlAs_2o_DS4_1WF7"
 ddkfile_2 = "AlAs_2o_DS5_1WF8"
 ddkfile_3 = "AlAs_2o_DS6_1WF9"
 gs_wfkpath = wfk0_path 

 !pmat = pmat(band 1, band 2, k, direction, spin)
 call get_opt_matel(pmat, gs_wfkpath, ddkfile_1, ddkfile_2, ddkfile_3, comm)
 print *, cryst%symrel_cart(:,:,:)
 

 ik = 6
 isym = 6 
 iband = 1
 fband = 3
 sband = 4
 krank = krank_new(ebands%nkpt, ebands%kptns)
 iperm = [(isym, isym=1,cryst%nsym)]
 print *, "iperm"
 print *, iperm

 do isym=1,cryst%nsym
   kpt = ebands%kptns(:,ik)
   symrec = cryst%symrec(:,:,isym)
   kpt_new = matmul(symrec, kpt)
   rank(isym) = krank%get_rank(kpt_new)
 end do
 call sort_int(cryst%nsym, rank, iperm)
 do isym=1,cryst%nsym
   if (rank(isym) == krank_prev) cycle
   sym = iperm(isym)
   krank_prev = rank(isym)
 end do


 ! Define (hard coded) optical mesh
 nw = 100
 wmin = zero
 wmax = 0.1
 ABI_MALLOC(wmesh, (nw))
 step = (wmax - wmin)/(nw - 1)
 wmesh = arth(wmin, step, nw)
 
 
 ABI_CALLOC(fbands, (ebands%mband))
 ABI_CALLOC(sbands, (ebands%mband))

 print *, "final bands"
   do ib1 = 1,ebands%mband
     if (all(ebands%eig(ib1,:,:) > ebands%fermie)) then 
       num_fband = num_fband + 1
       fbands(num_fband) = ib1
     end if
   end do
 print *, fbands

 print *, "Starting bands"
 do ib1 = 1,ebands%mband
   if (any(fbands == ib1)) cycle
   num_sband = num_sband + 1
   sbands(num_sband) = ib1
 end do
 print *, sbands

 ABI_MALLOC(bks_mask, (dtset%mband, ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(nband, (ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(wfd_istwfk, (ebands%nkpt))

 ! Read in all wave functions (memory intensive)
 bks_mask = .True.; keep_ur = .False.
 wfd_istwfk = 1
 nband = dtset%mband

 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, dtset%mband, nband, ebands%nkpt, ebands%nsppol, bks_mask,&
               dtset%nspden, ebands%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk,&
               ebands%kptns, ngfft, dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(nband)
 ABI_FREE(wfd_istwfk)


 call find_mpw(mpw, ebands%kptns, ebands%nsppol, ebands%nkpt, cryst%gmet,dtset%ecut,comm)
 dvqop = dvqop_new(dtset, dvdb, cryst, ifc, pawtab, psps, mpi_enreg, mpw, ngfft, ngfftf)
 
 rlatt = ebands%kptrlatt
 call matr3inv(rlatt, klatt)
 ABI_MALLOC(bz2ibz_indexes, (ebands%nkpt))
 bz2ibz_indexes = [(ik, ik=1,ebands%nkpt)]
 call htetra_init(htetra, bz2ibz_indexes, cryst%gprimd, klatt, ebands%kptns, & 
&                     ebands%nkpt, ebands%kptns, ebands%nkpt, ierr, errmsg, comm, 2)

 usecprj = 0
 ABI_MALLOC(cwaveprj0, (cryst%natom, ebands%nspinor*usecprj))

 ! ********** Allocate variables used in the loop **********
 ABI_MALLOC(cg_iband, (2, mpw*ebands%nspinor))
 ABI_MALLOC(cg_sband, (2, mpw*ebands%nspinor))
 ABI_MALLOC(cg_fband, (2, mpw*ebands%nspinor))
 ABI_MALLOC(detuning_minus, (dvqop%natom3))
 ABI_MALLOC(detuning_plus, (dvqop%natom3))
 ABI_MALLOC(energy_fs, (ebands%nkpt))
 ABI_MALLOC(ikq, (ebands%nkpt))
 ABI_MALLOC(gkq, (2,dvqop%natom3))
 ABI_MALLOC(numerator, (dvqop%natom3))
 ABI_CALLOC(path_sum_minus, (dvqop%natom3))
 ABI_CALLOC(path_sum_plus, (dvqop%natom3))
 ABI_MALLOC(path_summand_minus, (dvqop%natom3))
 ABI_MALLOC(path_summand_plus, (dvqop%natom3))
 ABI_MALLOC(path_summand_1_minus, (dvqop%natom3))
 ABI_MALLOC(path_summand_1_plus, (dvqop%natom3))
 ABI_MALLOC(path_summand_2_minus, (dvqop%natom3))
 ABI_MALLOC(path_summand_2_plus, (dvqop%natom3))
 ABI_MALLOC(phonon_energies, (dvqop%natom3,2))
 ABI_MALLOC(phonon_populations, (dvqop%natom3))
 ABI_MALLOC(transrate_integral, (nw))
 ABI_MALLOC(transrate_total, (nw))
 ABI_MALLOC(weights, (nw,dvqop%natom3,2,2))
 ! ********** End allocation **********

 call xmpi_split_work(ebands%nkpt, comm, my_start, my_stop)
 transrate_total = zero
 self_energy = (0.0, 0.0)
 my_start = 1
 my_stop = 1
 do spin=1,ebands%nsppol
   do iq=my_start,my_stop
     qpt = ebands%kptns(:,iq)
     ! Skip zero frequency phonon
     if (sum(abs(qpt)) < 1.0d-12) cycle
     call dvqop%setupq(cryst,qpt,spin,pawfgr,comm)
     phonon_populations = one/(exp(dvqop%phfreq/(9.89d-4)) - 1)

     do ik=1,ebands%nkpt
       kpt = ebands%kptns(:,ik)
       kqpt = kpt + qpt
       call findqg0(ikq(ik), g0_k, kqpt, ebands%nkpt, ebands%kptns, [1,1,1])
     end do !ik

     do ik=1,ebands%nkpt

       call dvqop%setup_spin_kpoint(dtset, cryst, psps, spin, kpt, kqpt, wfd%istwfk(ik), wfd%istwfk(ikq(ik)),&
                                    wfd%npwarr(ik), wfd%npwarr(ikq(ik)), wfd%kdata(ik)%kg_k, wfd%kdata(ikq(ik))%kg_k)

       ! Loop over starting (valence) bands
       do ib1=1,num_sband
         sband = sbands(ib1)
         ! Get sband wave function
         ABI_CHECK(mpw >= ebands%npwarr(ik), "mpw < npw_k")
         call wfd%copy_cg(ib1, ik, spin, cg_sband)

         ! Loop over final (conduction) bands
         do ib2=1,num_fband
           fband = fbands(ib2)
           if (fband > 8) cycle
           energy_fs = ebands%eig(fband,ikq,spin) - ebands%eig(sband,ik,spin)

           ! Most important line: ignore if not possible to make transition
           if (wmax + maxval(dvqop%phfreq) .lt. minval(energy_fs)) cycle

           ! Compute weights for phonon absorption (+) and emission (-)
           do branch=1,dvqop%natom3
             call htetra%get_onewk_wvals(ik, 0, nw, wmesh, one, ebands%nkpt, energy_fs - dvqop%phfreq(branch), weights(:,branch,2,1))
             call htetra%get_onewk_wvals(ik, 0, nw, wmesh, one, ebands%nkpt, energy_fs + dvqop%phfreq(branch), weights(:,branch,1,1))
             weights = weights/ebands%nkpt
           end do

           ! Get fband wave function
           ABI_CHECK(mpw >= ebands%npwarr(ikq(ik)), "mpw < npw_k")
           call wfd%copy_cg(fband, ikq(ik), spin, cg_fband)

           ! Initialize path sums to zero
           path_sum_plus = zero
           path_sum_minus = zero

           ! Loop over intermediate bands (can be at k or k+q depending on order)
           do iband=1,ebands%mband - 2
             call wfd%copy_cg(sband, ik, spin, cg_iband)

             !Phonon first
             gkq = dvqop%get_gkq(ebands%eig(iband,ikq(ik),spin), wfd%istwfk(ik), wfd%npwarr(ik), wfd%istwfk(ikq(ik)),&
                 wfd%npwarr(ikq(ik)), ebands%nspinor, cg_iband, cg_sband, cwaveprj0)
             numerator = pmat(fband,iband,ikq(ik),1,spin)*CMPLX(gkq(1,:), gkq(2,:), kind=dpc)

             detuning_plus = ebands%eig(iband,ikq(ik),spin) - ebands%eig(sband,ik,spin) + dvqop%phfreq + self_energy 
             detuning_minus = ebands%eig(iband,ikq(ik),spin) - ebands%eig(sband,ik,spin) - dvqop%phfreq + self_energy 

             path_summand_1_plus = numerator/detuning_plus
             path_summand_1_minus = numerator/detuning_minus

             !Photon first
             gkq = dvqop%get_gkq(ebands%eig(fband,ikq(ik),spin), wfd%istwfk(ik), wfd%npwarr(ik), wfd%istwfk(ikq(ik)),&
                 wfd%npwarr(ikq(ik)), ebands%nspinor, cg_fband, cg_iband, cwaveprj0)
             numerator = CMPLX(gkq(1,:), gkq(2,:), kind=dpc)*pmat(iband,sband,ikq(ik),1,spin)

             detuning_plus = ebands%eig(iband,ik,spin) - ebands%eig(fband,ikq(ik),spin) + dvqop%phfreq + self_energy 
             detuning_minus = ebands%eig(iband,ik,spin) - ebands%eig(fband,ikq(ik),spin) - dvqop%phfreq + self_energy 

             path_summand_2_plus = numerator/detuning_plus
             path_summand_2_minus = numerator/detuning_minus

             ! Accumulate plus and minus paths separately for each pert order
             path_summand_plus = path_summand_1_plus + path_summand_2_plus
             path_summand_minus = path_summand_1_plus + path_summand_2_plus

             ! Add total (both pert orders) contr. involving iband
             path_sum_plus = path_sum_plus + path_summand_plus
             path_sum_minus = path_sum_minus + path_summand_minus

           end do !iband

           ! Compute the integral with appropriate weights and phonon occs
           ! TODO: Fix electron occupations to replace factor of 2
           transrate_integral =( &
                                matmul(abs(path_sum_plus)**2,transpose(weights(:,:,1,1))) &
&                             + matmul(abs(path_sum_minus)**2*phonon_populations, transpose(weights(:,:,2,1))) &
&                              )*(2)

           if (any(transrate_integral .lt. 0)) then
           print *, "electron pops"
           print *, ebands%occ(ik, sband, spin) - ebands%occ(ikq(ik), fband, spin)

           print *, "phonon pops"
           print *, phonon_populations
           end if

           transrate_total = transrate_total + transrate_integral
          end do !fband
        end do !sband
     end do !spin
   end do !ik
 end do !iq
 call xmpi_sum(transrate_total, comm, ierr)
 if (ierr /= 0) MSG_ERROR("Error in xmpi sum for transition rate")
 print *, "total trans"
 print *, transrate_total
 
 call htetra%free
 ! ********** Free loop vars ********
 ABI_FREE(cg_iband)
 ABI_FREE(cg_fband)
 ABI_FREE(cg_sband)
 ABI_FREE(detuning_minus)
 ABI_FREE(detuning_plus)
 ABI_FREE(energy_fs)
 ABI_FREE(gkq)
 ABI_FREE(ikq)
 ABI_FREE(numerator)
 ABI_FREE(path_sum_minus)
 ABI_FREE(path_sum_plus)
 ABI_FREE(path_summand_minus)
 ABI_FREE(path_summand_plus)
 ABI_FREE(path_summand_1_minus)
 ABI_FREE(path_summand_1_plus)
 ABI_FREE(path_summand_2_minus)
 ABI_FREE(path_summand_2_plus)
 ABI_FREE(transrate_integral)
 ABI_FREE(transrate_total)
 ABI_FREE(phonon_energies)
 ABI_FREE(phonon_populations)
 ABI_FREE(weights)
 ! ********** end loop vars **********

 call pawcprj_free(cwaveprj0)

 ABI_FREE(bz2ibz_indexes)
 ABI_FREE(sbands)
 ABI_FREE(fbands)
 ABI_FREE(wmesh)

end subroutine absrate_ind2

!!***

!!****f* m_gkk/eph_gkk
!! NAME
!!  eph_gkk
!!
!! FUNCTION
!!  Compute electron-phonon coupling matrix elements.
!!
!! INPUTS
!! wk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! PARENTS
!!      m_eph_driver
!!
!! NOTES
!!
!! CHILDREN
!!      get_kg
!!
!! SOURCE

subroutine eph_gkk(wfk0_path,wfq_path,dtfil,ngfft,ngfftf,dtset,cryst,ebands_k,ebands_kq,dvdb,ifc,&
                       pawfgr,pawang,pawrad,pawtab,psps,mpi_enreg,comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path, wfq_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands_k, ebands_kq
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

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_getgh1c=1,berryopt0=0, useylmgr1=0,master=0
 integer :: my_rank,nproc,mband,mband_kq,my_minb,my_maxb,nsppol,nkpt,nkpt_kq,idir,ipert
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor
 integer :: ib1,ib2,band,ik,ikq,timerev_q
 integer :: spin,istwf_k,istwf_kq,npw_k,npw_kq, comm_rpt
 integer :: mpw,mpw_k,mpw_kq,ierr,my_kstart,my_kstop,ncid
 integer :: n1,n2,n3,n4,n5,n6,nspden,ncerr
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1, interpolated
 real(dp) :: cpu,wall,gflops,ecut,eshift,eig0nk,dotr,doti
 logical :: i_am_master, gen_eigenpb
 type(wfd_t) :: wfd_k, wfd_kq
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(gkk_t) :: gkk2d
 character(len=500) :: msg, what
 character(len=fnlen) :: fname, gkkfilnam
!arrays
 integer :: g0_k(3),symq(4,2,cryst%nsym)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),nband(:,:),nband_kq(:,:),blkflg(:,:), wfd_istwfk(:)
 real(dp) :: kk(3),kq(3),qpt(3),phfrq(3*cryst%natom)
 real(dp),allocatable :: displ_cart(:,:,:),displ_red(:,:,:), eigens_kq(:,:,:)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: v1scf(:,:,:,:),gkk(:,:,:,:,:)
 real(dp),allocatable :: bras(:,:,:),kets(:,:,:),h1_kets(:,:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: dummy_vtrial(:,:),gvnlx1(:,:), gs1c(:,:), gkq_atm(:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),bks_mask_kq(:,:,:),keep_ur(:,:,:),keep_ur_kq(:,:,:)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:) !natom,nspinor*usecprj)

!************************************************************************

 what = "(GKK files)"; if (dtset%eph_task == -2) what = "GKQ file"
 write(msg, '(3a)') " Computation of electron-phonon coupling matrix elements ", trim(what), ch10
 call wrtout([std_out, ab_out], msg, do_flush=.True.)

 if (psps%usepaw == 1) then
   MSG_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm); i_am_master = my_rank == master

 ! Copy important dimensions
 natom = cryst%natom
 natom3 = 3 * natom
 nsppol = ebands_k%nsppol
 nspinor = ebands_k%nspinor
 nspden = dtset%nspden
 nkpt = ebands_k%nkpt
 mband = ebands_k%mband
 nkpt_kq = ebands_kq%nkpt
 mband_kq = ebands_kq%mband
 ecut = dtset%ecut
 !write(std_out, *)"ebands dims (b, k, s): ", ebands_k%mband, ebands_k%nkpt, ebands_k%nsppol
 !write(std_out, *)"ebands_kq dims (b, k, s): ", ebands_kq%mband, ebands_kq%nkpt, ebands_kq%nsppol

 qpt = dtset%qptn(:)

 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 n4=ngfft(4); n5=ngfft(5); n6=ngfft(6)

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 ! Initialize the wave function descriptors.
 ! For the time being, no memory distribution, each node has the full set of states.
 my_minb = 1; my_maxb = mband

 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask,(mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur,(mband, nkpt ,nsppol))
 nband=mband; bks_mask=.False.; keep_ur=.False.

 ABI_MALLOC(nband_kq, (nkpt_kq, nsppol))
 ABI_MALLOC(bks_mask_kq,(mband_kq, nkpt_kq, nsppol))
 ABI_MALLOC(keep_ur_kq,(mband_kq, nkpt_kq ,nsppol))
 nband_kq=mband_kq; bks_mask_kq=.False.; keep_ur_kq=.False.

 ! Distribute the k-points over the processors
 call xmpi_split_work(nkpt,comm,my_kstart,my_kstop)
 do ik=1,nkpt
   if (.not. (ik >= my_kstart .and. ik <= my_kstop)) cycle
   kk = ebands_k%kptns(:,ik)
   kq = kk + qpt
   ! Find the index of the k+q point
   call findqg0(ikq,g0_k,kq,nkpt_kq,ebands_kq%kptns(:,:), [1,1,1])
   bks_mask(:,ik,:) = .True.
   bks_mask_kq(:,ikq,:) = .True.
 end do

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1

 ! Initialize the wavefunction descriptors
 call wfd_init(wfd_k,cryst,pawtab,psps,keep_ur,mband,nband,nkpt,nsppol,bks_mask,&
   nspden,nspinor,ecut,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands_k%kptns,ngfft,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)
 ABI_FREE(wfd_istwfk)

 call wfd_k%print(header="Wavefunctions on the k-points grid",mode_paral='PERS')

 ABI_MALLOC(wfd_istwfk, (nkpt_kq))
 wfd_istwfk = 1

 call wfd_init(wfd_kq,cryst,pawtab,psps,keep_ur_kq,mband_kq,nband_kq,nkpt_kq,nsppol,bks_mask_kq,&
   nspden,nspinor,ecut,dtset%ecutsm,dtset%dilatmx,wfd_istwfk,ebands_kq%kptns,ngfft,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)
 ABI_FREE(wfd_istwfk)

 call wfd_kq%print(header="Wavefunctions on the q-shifted k-points grid",mode_paral='PERS')

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(nband_kq)
 ABI_FREE(bks_mask_kq)
 ABI_FREE(keep_ur_kq)

 ! Read wavefunctions on the k-points grid and q-shifted k-points grid.
 call wfd_k%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 call wfd_kq%read_wfk(wfq_path, iomode_from_fname(wfq_path))

 ! ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx,natom,n1,n2,n3,ph1d,cryst%xred)

 ! Find the appropriate value of mpw
 call find_mpw(mpw_k, ebands_k%kptns(:,:), nsppol, nkpt, cryst%gmet,ecut,comm)
 call find_mpw(mpw_kq, ebands_kq%kptns(:,:), nsppol, nkpt_kq, cryst%gmet,ecut,comm)
 mpw = max(mpw_k, mpw_kq)

 ! Allow PW-arrays dimensioned with mpw
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))

 ! Spherical Harmonics for useylm==1.
 ABI_MALLOC(ylm_k,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylm_kq,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylmgr_kq,(mpw, 3, psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

 ! TODO FOR PAW
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (natom, nspinor*usecprj))

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1  ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2     ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0 ! gvnlx1 is output
 ABI_MALLOC(gvnlx1, (2,usevnl))
 ABI_MALLOC(grad_berry, (2,nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 !1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 !2) Perform the setup needed for the non-local factors:
 !* Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 !* PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq,psps,pawtab,nspinor,nsppol,nspden,natom,&
   dtset%typat,cryst%xred,nfft,mgfft,ngfft,cryst%rprimd,dtset%nloalg,&
   usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda,&
   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab)

 ! Allocate vlocal. Note nvloc
 ABI_MALLOC(vlocal,(n4,n5,n6,gs_hamkq%nvloc))
 ! Allocate work space arrays.
 ABI_MALLOC(blkflg, (natom3,natom3))
 ABI_CALLOC(dummy_vtrial, (nfftf,nspden))

 call cwtime(cpu, wall, gflops, "start")

 interpolated = 0
 if (dtset%eph_use_ftinterp /= 0) then
   MSG_WARNING(sjoin("Enforcing FT interpolation for q-point", ktoa(qpt)))
   comm_rpt = xmpi_comm_self
   call dvdb%ftinterp_setup(dtset%ddb_ngqpt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)
   cplex = 2
   ABI_MALLOC(v1scf, (cplex, nfftf, nspden, dvdb%my_npert))
   call dvdb%ftinterp_qpt(qpt, nfftf, ngfftf, v1scf, dvdb%comm_rpt)
   interpolated = 1
 else
   ! Find the index of the q-point in the DVDB.
   db_iqpt = dvdb%findq(qpt)
   if (db_iqpt /= -1) then
     if (dtset%prtvol > 0) call wrtout(std_out, sjoin("Found: ",ktoa(qpt)," in DVDB with index ",itoa(db_iqpt)))
     ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
     ! This call allocates v1scf(cplex, nfftf, nspden, 3*natom))
     call dvdb%readsym_allv1(db_iqpt, cplex, nfftf, ngfftf, v1scf, comm)
   else
     MSG_WARNING(sjoin("Cannot find q-point:", ktoa(qpt), "in DVDB file"))
   end if
 end if

 ! Examine the symmetries of the q wavevector
 call littlegroup_q(cryst%nsym,qpt,symq,cryst%symrec,cryst%symafm,timerev_q,prtvol=dtset%prtvol)

 ! Allocate vlocal1 with correct cplex. Note nvloc
 ABI_MALLOC_OR_DIE(vlocal1,(cplex*n4,n5,n6,gs_hamkq%nvloc,natom3), ierr)

 ABI_MALLOC(displ_cart, (2,3*cryst%natom,3*cryst%natom))
 ABI_MALLOC(displ_red, (2,3*cryst%natom,3*cryst%natom))

 if (dtset%eph_task == 2) then
   ! Write GKK files (1 file for perturbation)
   ABI_MALLOC(gkk, (2*mband*nsppol,nkpt,1,1,mband_kq))

 else if (dtset%eph_task == -2) then
   ! Write GKQ file with all perturbations. gkq are given in the atom representation.
   ! TODO: mband_kq == mband
   ABI_MALLOC(gkq_atm, (2, mband_kq, mband, nkpt))
   if (i_am_master) then
     call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)
     fname = strcat(dtfil%filnam_ds(4), "_GKQ.nc")
#ifdef HAVE_NETCDF
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKQ file")
     NCF_CHECK(cryst%ncwrite(ncid))
     ! Write bands on k mesh.
     NCF_CHECK(ebands_ncwrite(ebands_k, ncid))
     ncerr = nctk_def_dims(ncid, [nctkdim_t('number_of_phonon_modes', natom3)], defmode=.True.)
     NCF_CHECK(ncerr)
     ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
       "symdynmat", "symv1scf", "dvdb_add_lr", "interpolated"])
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: "qdamp"]))

     ! Define EPH arrays
     ncerr = nctk_def_arrays(ncid, [ &
       nctkarr_t('qpoint', "dp" , 'number_of_reduced_dimensions'), &
       nctkarr_t('emacro_cart', "dp", 'number_of_cartesian_directions, number_of_cartesian_directions'), &
       nctkarr_t('becs_cart', "dp", "number_of_cartesian_directions, number_of_cartesian_directions, number_of_atoms"), &
       nctkarr_t("eigenvalues_kq", "dp", "max_number_of_states, number_of_kpoints, number_of_spins"), &
       nctkarr_t('phfreqs', "dp", 'number_of_phonon_modes'), &
       nctkarr_t('phdispl_cart', "dp", 'complex, number_of_phonon_modes, number_of_phonon_modes'), &
       !nctkarr_t('phdispl_cart_qvers', "dp", 'complex, number_of_phonon_modes, number_of_phonon_modes'), &
       nctkarr_t('phdispl_red', "dp", 'complex, number_of_phonon_modes, number_of_phonon_modes'), &
       nctkarr_t("gkq_representation", "char", "character_string_length"), &
       nctkarr_t('gkq', "dp", &
         'complex, max_number_of_states, max_number_of_states, number_of_phonon_modes, number_of_kpoints, number_of_spins') &
     ])
     NCF_CHECK(ncerr)
     ! Write data.
     NCF_CHECK(nctk_set_datamode(ncid))
     ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
       "symdynmat", "symv1scf", "dvdb_add_lr", "interpolated"], &
       [dtset%symdynmat, dtset%symv1scf, dtset%dvdb_add_lr, interpolated])
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: "qdamp"], [dvdb%qdamp]))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpoint"), qpt))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "emacro_cart"), dvdb%dielt))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "becs_cart"), dvdb%zeff))
     ABI_MALLOC(eigens_kq, (ebands_kq%mband, nkpt, nsppol))
     do ik=1,nkpt
       kk = ebands_k%kptns(:,ik)
       kq = kk + qpt
       ! Find the index of the k+q point
       call findqg0(ikq, g0_k, kq, nkpt_kq, ebands_kq%kptns, [1,1,1])
       eigens_kq(:, ik, :) = ebands_kq%eig(:, ikq, :)
     end do
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eigenvalues_kq"), eigens_kq))
     ABI_FREE(eigens_kq)
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phfreqs"), phfrq))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phdispl_cart'), displ_cart))
     ! Add phonon displacement for qvers
     !call ifc_fourq(ifc, cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red, nanaqdir="reduced")
     !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phdispl_cart_qvers'), displ_cart))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, 'phdispl_red'), displ_red))
     NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "gkq_representation"), "atom"))
#endif
   end if
 else
   MSG_ERROR(sjoin("Invalid value for eph_task:", itoa(dtset%eph_task)))
 end if

 ! Loop over all 3*natom perturbations.
 do ipc=1,natom3
   idir = mod(ipc-1, 3) + 1
   ipert = (ipc - idir) / 3 + 1
   write(msg, '(a,2i4)') " Treating ipert, idir = ", ipert, idir
   call wrtout(std_out, msg, do_flush=.True.)
   if (dtset%eph_task == 2) gkk = zero

   do spin=1,nsppol
     if (dtset%eph_task == -2) gkq_atm = zero

     ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
     call rf_transgrid_and_pack(spin,nspden,psps%usepaw,cplex,nfftf,nfft,ngfft,gs_hamkq%nvloc,&
               pawfgr,mpi_enreg,dummy_vtrial,v1scf(:,:,:,ipc),vlocal,vlocal1(:,:,:,:,ipc))

     ! Continue to initialize the Hamiltonian
     call gs_hamkq%load_spin(spin,vlocal=vlocal,with_nonlocal=.true.)

     ! Allocate workspace for wavefunctions. Make npw larger than expected.
     ABI_MALLOC(bras, (2, mpw*nspinor, mband))
     ABI_MALLOC(kets, (2, mpw*nspinor, mband))
     ABI_MALLOC(h1_kets, (2, mpw*nspinor, mband))

     ! GKA: This little block used to be right after the perturbation loop
     ! Prepare application of the NL part.
     call init_rf_hamiltonian(cplex,gs_hamkq,ipert,rf_hamkq,has_e1kbsc=.true.)
     call rf_hamkq%load_spin(spin,vlocal1=vlocal1(:,:,:,:,ipc),with_nonlocal=.true.)

     do ik=1,nkpt
       ! Only do a subset a k-points
       if (.not. (ik >= my_kstart .and. ik <= my_kstop)) cycle

       kk = ebands_k%kptns(:,ik)
       kq = kk + qpt
       ! Find the index of the k+q point
       call findqg0(ikq, g0_k, kq, nkpt_kq, ebands_kq%kptns, [1,1,1])

       ! Copy u_k(G)
       istwf_k = wfd_k%istwfk(ik); npw_k = wfd_k%npwarr(ik)
       ABI_CHECK(mpw >= npw_k, "mpw < npw_k")
       kg_k(:,1:npw_k) = wfd_k%kdata(ik)%kg_k
       do ib2=1,mband
         call wfd_k%copy_cg(ib2, ik, spin, kets(1,1,ib2))
       end do

       ! Copy u_kq(G)
       istwf_kq = wfd_kq%istwfk(ikq); npw_kq = wfd_kq%npwarr(ikq)
       ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
       kg_kq(:,1:npw_kq) = wfd_kq%kdata(ikq)%kg_k
       do ib1=1,mband_kq
         call wfd_kq%copy_cg(ib1, ikq, spin, bras(1,1,ib1))
       end do

       ! if PAW, one has to solve a generalized eigenproblem
       ! Be careful here because I will need sij_opt==-1
       gen_eigenpb = (psps%usepaw==1)
       sij_opt = 0; if (gen_eigenpb) sij_opt = 1
       ABI_MALLOC(gs1c, (2,npw_kq*nspinor*((sij_opt+1)/2)))

       ! GKA: Previous loop on 3*natom perturbations used to start here
       ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
       call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert,&    ! In
         cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k,&            ! In
         npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq,&           ! In
         dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)       ! Out

       ! Calculate dvscf * psi_k, results stored in h1_kets on the k+q sphere.
       ! Compute H(1) applied to GS wavefunction Psi(0)
       do ib2=1,mband
         eig0nk = ebands_k%eig(ib2,ik,spin)
         ! Use scissor shift on 0-order eigenvalue
         eshift = eig0nk - dtset%dfpt_sciss

         call getgh1c(berryopt0,kets(:,:,ib2),cwaveprj0,h1_kets(:,:,ib2),&
                      grad_berry,gs1c,gs_hamkq,gvnlx1,idir,ipert,eshift,mpi_enreg,optlocal,&
                      optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
       end do

       ABI_FREE(kinpw1)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)
       ABI_FREE(dkinpw)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(ph3d)
       ABI_FREE(gs1c)
       ABI_SFREE(ph3d1)

       ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
       ! The array eig1_k contains:
       !
       ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
       ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
       do ib2=1,mband
         do ib1=1,mband_kq
           call dotprod_g(dotr,doti,istwf_kq,npw_kq*nspinor,2,bras(1,1,ib1),h1_kets(1,1,ib2),&
             mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
           band = 2*ib2-1 + (spin-1) * 2 * mband

           if (dtset%eph_task == 2) then
             gkk(band,ik,1,1,ib1) = dotr
             gkk(band+1,ik,1,1,ib1) = doti
           else
             gkq_atm(:, ib1, ib2, ik) = [dotr, doti]
           end if

         end do
       end do

     end do ! ikpt

     ABI_FREE(bras)
     ABI_FREE(kets)
     ABI_FREE(h1_kets)
     call rf_hamkq%free()

     if (dtset%eph_task == -2) then
       ! Gather the k-points computed by all processes
       call xmpi_sum_master(gkq_atm, master, comm, ierr)
       if (i_am_master) then
         ! Write the netCDF file.
#ifdef HAVE_NETCDF
         ncerr = nf90_put_var(ncid, nctk_idname(ncid, "gkq"), gkq_atm, &
           start=[1, 1, 1, ipc, 1, 1], count=[2, mband, mband, 1, nkpt, spin])
         NCF_CHECK(ncerr)
#endif
       end if
     end if
   end do ! spin

   if (dtset%eph_task == 2) then
     ! Gather the k-points computed by all processes
     call xmpi_sum_master(gkk,master,comm,ierr)
     ! Init a gkk_t object
     call gkk_init(gkk,gkk2d,mband,nsppol,nkpt,1,1)
     ! Write the netCDF file.
     call appdig(ipc,dtfil%fnameabo_gkk,gkkfilnam)
     fname = strcat(gkkfilnam, ".nc")
#ifdef HAVE_NETCDF
     if (i_am_master) then
       NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKK file")
       NCF_CHECK(cryst%ncwrite(ncid))
       NCF_CHECK(ebands_ncwrite(ebands_k, ncid))
       call gkk_ncwrite(gkk2d, qpt, 1.0_dp,  ncid)
       NCF_CHECK(nf90_close(ncid))
     end if
#endif
     ! Free memory
     call gkk_free(gkk2d)
   end if
 end do ! ipc (loop over 3*natom atomic perturbations)

 call cwtime(cpu, wall, gflops, "stop")
 write(msg, '(2a)') " Computation of gkq matrix elements with ", trim(what)
 call wrtout([std_out, ab_out], msg, do_flush=.True.)
 call wrtout(std_out, sjoin("cpu-time:", sec2str(cpu), ",wall-time:", sec2str(wall)), do_flush=.True.)

 if (dtset%eph_task == -2 .and. i_am_master) then
#ifdef HAVE_NETCDF
   NCF_CHECK(nf90_close(ncid))
#endif
 end if

 ! ===========
 ! Free memory
 ! ===========
 ABI_SFREE(gkk)
 ABI_SFREE(gkq_atm)
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red)
 ABI_FREE(v1scf)
 ABI_FREE(vlocal1)
 ABI_FREE(gvnlx1)
 ABI_FREE(grad_berry)
 ABI_FREE(dummy_vtrial)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr_kq)
 ABI_FREE(blkflg)

 call gs_hamkq%free()
 call wfd_k%free()
 call wfd_kq%free()
 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)

end subroutine eph_gkk
!!***

!----------------------------------------------------------------------

!!****f* m_gkk/ncwrite_v1qnu
!! NAME
!!  ncwrite_v1qnu
!!
!! FUNCTION
!!  Compute \delta V_{q,nu)(r) and dump results to netcdf file.
!!  This routine should be called by a single processor.
!!
!! INPUT
!!  dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!!  dtset<dataset_type>= Input variables.
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!  out_ncpath=Name of the netcdf file.
!!
!! OUTPUT
!!  Only writing
!!
!! PARENTS
!!      m_eph_driver
!!
!! CHILDREN
!!      get_kg
!!
!! SOURCE

subroutine ncwrite_v1qnu(dvdb, dtset, ifc, out_ncpath)

 use m_bz_mesh, only : kpath_t, kpath_new

!Arguments ------------------------------------
!scalars
 type(dataset_type),target,intent(in) :: dtset
 type(dvdb_t),intent(inout) :: dvdb
 type(ifc_type),intent(in) :: ifc
 character(len=*),intent(in) :: out_ncpath

!Local variables-------------------------------
!scalars
 integer,parameter :: master = 0
 integer :: db_iqpt, cplex, nfft, comm, ip, idir, ipert, my_rank, interpolated, comm_rpt
 logical :: with_lr_model
#ifdef HAVE_NETCDF
 integer :: ncid, ncerr
#endif
!arrays
 integer :: ngfft(18)
 real(dp) :: phfreqs(dvdb%natom3),qpt(3)
 real(dp) :: displ_cart(2,3, dvdb%cryst%natom, dvdb%natom3), displ_red(2,dvdb%natom3,dvdb%natom3)
 real(dp),allocatable :: v1scf(:,:,:,:), v1_qnu(:,:,:,:)
 real(dp),allocatable :: v1lr_atm(:,:,:,:), v1lr_qnu(:,:,:,:)
#if 1
 integer :: iq, nu, iatom, ii, jj, kk
 real(dp) :: inv_qepsq, qtau, phre, phim, rtmp
 type(kpath_t) :: qpath
 real(dp) :: bounds(3,6), qpt_red(3), qpt_cart(3),  glr(3)
 real(dp) :: values(dvdb%natom3)

!************************************************************************

 ! +0.50000  +0.50000  +0.50000  # L
 ! +0.00000  +0.00000  +0.00000  # $\Gamma$
 ! +0.50000  +0.00000  +0.50000  # X
 ! +0.50000  +0.25000  +0.75000  # W
 ! +0.37500  +0.37500  +0.75000  # K
 ! +0.00000  +0.00000  +0.00000  # $\Gamma$
 ! +0.37500  +0.37500  +0.75000  # K

 ! +0.62500  +0.25000  +0.62500  # U
 ! +0.50000  +0.50000  +0.50000  # L
 ! +0.37500  +0.37500  +0.75000  # K
 ! +0.62500  +0.25000  +0.62500  # U
 ! +0.50000  +0.00000  +0.50000  # X

 bounds(:, 1) = tol3 * [+0.50000,  +0.50000, +0.50000] !  # L
 bounds(:, 2) = tol3 * [+0.00000,  +0.00000, +0.00000] !  # $\Gamma$
 bounds(:, 3) = tol3 * [+0.50000,  +0.00000, +0.50000] !  # X
 bounds(:, 4) = tol3 * [+0.37500,  +0.37500, +0.75000] !  # K
 bounds(:, 5) = tol3 * [+0.00000,  +0.00000, +0.00000] !  # $\Gamma$
 bounds(:, 6) = tol3 * [+0.50000,  +0.25000, +0.75000] !  # W

 qpath = kpath_new(bounds, dvdb%cryst%gprimd, dtset%ndivsm)

 do iq=1,qpath%npts
   qpt_red = qpath%points(:, iq)
   qpt_cart = two_pi * matmul(dvdb%cryst%gprimd, qpt_red)
   inv_qepsq = one / dot_product(qpt_cart, matmul(ifc%dielt, qpt_cart))
   call ifc%fourq(dvdb%cryst, qpt_red, phfreqs, displ_cart)
   do nu=1, dvdb%natom3
     glr = zero
     do iatom=1, dvdb%cryst%natom
       ! Phase factor exp(-i (q+G) . tau)
       qtau = - two_pi * dot_product(qpt_red, dvdb%cryst%xred(:,iatom))
       phre = cos(qtau); phim = sin(qtau)
       do jj=1,3
         do ii=1,3
           do kk=1,3
             rtmp = dvdb%qstar(ii, jj, kk, iatom) * qpt_cart(ii) * qpt_cart(jj)
             glr(1) = glr(1) + rtmp * (displ_cart(1, kk, iatom, nu) * phre - displ_cart(2, kk, iatom, nu) * phim)
             glr(2) = glr(2) + rtmp * (displ_cart(2, kk, iatom, nu) * phre + displ_cart(1, kk, iatom, nu) * phre)
           end do
         end do
       end do
     end do
     glr = half * (glr / inv_qepsq) * (four_pi / dvdb%cryst%ucvol)
     values(nu) = (glr(1) ** 2 + glr(2) ** 2) / (two *  phfreqs(nu))
   end do ! nu
   write(std_out, "(i0, 4(f9.6), /, (es18.6, 1x))") iq, qpt_red, phfreqs(nu), (values(nu), nu=1, 3*dvdb%natom)
 end do ! iqpt

 call qpath%free()
 return
#endif

 my_rank = xmpi_comm_rank(dvdb%comm)
 comm = dvdb%comm
 qpt = dtset%qptn

 call wrtout(std_out, sjoin(" Writing Delta V_{q,nu)(r) potentials to file:", out_ncpath), do_flush=.True.)
 call wrtout([std_out, ab_out], sjoin(ch10, "- Results stored in: ", out_ncpath))
 call wrtout(std_out, sjoin(" Using qpt:", ktoa(qpt)))
 !call wrtout([std_out, ab_out], " Use `abiopen.py out_V1QAVG.nc -e` to visualize results")
 call dvdb%print(unit=std_out)

 ! Define FFT mesh
 ngfft = dvdb%ngfft
 nfft = product(ngfft(1:3))

 if (dtset%eph_task == -16) then
   call wrtout([std_out, ab_out], " Assuming q-point already in the DVDB file. No interpolation.")
   interpolated = 0

 else if (dtset%eph_task == +16) then
   call wrtout([std_out, ab_out], " Using Fourier interpolation.")
    comm_rpt = xmpi_comm_self
    call dvdb%ftinterp_setup(dtset%ddb_ngqpt, 1, dtset%ddb_shiftq, nfft, ngfft, comm_rpt)
    interpolated = 1
 else
   MSG_ERROR(sjoin("Invalid value for eph_task:", itoa(dtset%eph_task)))
 end if

 with_lr_model = .True.

 ! Create netcdf file.
#ifdef HAVE_NETCDF
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, out_ncpath, comm))
   NCF_CHECK(dvdb%cryst%ncwrite(ncid))

   ! Add other dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nfft", nfft), nctkdim_t("nspden", dvdb%nspden), &
     nctkdim_t("natom3", 3 * dvdb%cryst%natom)], defmode=.True.)
   NCF_CHECK(ncerr)

   ! Define arrays
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("ngfft", "int", "three"), &
     nctkarr_t("qpt", "dp", "three"), &
     nctkarr_t("phfreqs", "dp", "natom3"), &
     nctkarr_t("displ_cart", "dp", "two, natom3, natom3"), &
     nctkarr_t("v1_qnu", "dp", "two, nfft, nspden, natom3")])
   NCF_CHECK(ncerr)

   if (with_lr_model) then
     NCF_CHECK(nctk_def_arrays(ncid, [nctkarr_t("v1lr_qnu", "dp", "two, nfft, nspden, natom3")]))
   end if

   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngfft"), ngfft(1:3)))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "qpt"), qpt))
 end if
#endif

 ABI_MALLOC(v1_qnu, (2, nfft, dvdb%nspden, dvdb%natom3))
 if (with_lr_model) then
   ABI_MALLOC(v1lr_atm, (2, nfft, dvdb%nspden, dvdb%natom3))
   ABI_MALLOC(v1lr_qnu, (2, nfft, dvdb%nspden, dvdb%natom3))
 end if

 ! Get phonon freqs and displacemented for this q-point.
 call ifc%fourq(dvdb%cryst, qpt, phfreqs, displ_cart, out_displ_red=displ_red)

 if (interpolated == 0) then
   ! Find the index of the q-point in the DVDB.
   db_iqpt = dvdb%findq(qpt)
   if (db_iqpt /= -1) then
     ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
     ! This call allocates v1scf(cplex, nfft, nspden, 3*natom))
     call dvdb%readsym_allv1(db_iqpt, cplex, nfft, ngfft, v1scf, comm)
   else
     MSG_ERROR(sjoin("Cannot find q-point:", ktoa(qpt), "in DVDB file"))
   end if
 else

   cplex = 2
   ABI_MALLOC(v1scf, (cplex, nfft, dvdb%nspden, dvdb%my_npert))
   call dvdb%ftinterp_qpt(qpt, nfft, ngfft, v1scf, dvdb%comm_rpt)
 end if

 ! Compute scattering potential in phonon representations instead ot atomic one.
 ! v1_qnu = \sum_{ka} phdispl{ka}(q,nu) D_{ka,q} V_scf(r)
 ! NOTE: prefactor 1/sqrt(2 w(q,nu)) is not included in the potentials saved to file.
 ! v1_qnu(2, nfft, nspden, natom3), v1scf(cplex, nfft, nspden, natom3)
 call v1atm_to_vqnu(cplex, nfft, dvdb%nspden, dvdb%natom3, v1scf, displ_red, v1_qnu)

 if (with_lr_model) then
   ! Compute LR model in the atomic representation then compute phonon representation in v1lr_qnu.
   v1lr_atm = zero
   do idir=1,3
     do ipert=1,dvdb%natom
       ip = (ipert - 1) * 3 + idir
       call dvdb%get_v1r_long_range(qpt, idir, ipert, nfft, ngfft, v1lr_atm(:,:,1,ip))
       if (dvdb%nspden == 2) v1lr_atm(:,:,2,ip) = v1lr_atm(:,:,1,ip)
     end do
   end do
   call v1atm_to_vqnu(2, nfft, dvdb%nspden, dvdb%natom3, v1lr_atm, displ_red, v1lr_qnu)
 end if

#ifdef HAVE_NETCDF
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "phfreqs"), phfreqs))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "displ_cart"), displ_cart))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "v1_qnu"), v1_qnu))
 if (with_lr_model) then
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "v1lr_qnu"), v1lr_qnu))
 end if
#endif

 ABI_FREE(v1scf)
 ABI_FREE(v1_qnu)
 ABI_SFREE(v1lr_atm)
 ABI_SFREE(v1lr_qnu)

#ifdef HAVE_NETCDF
 NCF_CHECK(nf90_close(ncid))
#endif
 call dvdb%close()

 call wrtout(std_out, "dvqnu file written", do_flush=.True.)

end subroutine ncwrite_v1qnu
!!***

!----------------------------------------------------------------------

!!****f* m_gkk/v1atm_to_vqnu
!! NAME
!!  v1atm_to_vqnu
!!
!! FUNCTION
!!  Receive potentials in atomic representation and return potential in phonon representation
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gkk
!!
!! CHILDREN
!!      get_kg
!!
!! SOURCE

subroutine v1atm_to_vqnu(cplex, nfft, nspden, natom3, v1_atm, displ_red, v1_qnu)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex, nfft, nspden, natom3
!arrays
 real(dp),intent(in) :: v1_atm(cplex, nfft, nspden, natom3)
 real(dp),intent(out) :: v1_qnu(2, nfft, nspden, natom3)
 real(dp),intent(in) :: displ_red(2, natom3, natom3)

!Local variables-------------------------------
!scalars
 integer :: nu, ip, ispden

!************************************************************************

 do nu=1,natom3
   ! v1_qnu = \sum_{ka} phdispl{ka}(q,nu) D_{ka,q} V_scf(r)
   ! NOTE: prefactor 1/sqrt(2 w(q,nu)) is not included in the potentials saved to file.
   ! v1_qnu(2, nfft, nspden, natom3), v1_atm(cplex, nfft, nspden, natom3)
   v1_qnu(:, :, :, nu) = zero
   do ip=1,natom3
     do ispden=1,nspden
       if (cplex == 2) then
         v1_qnu(1, :, ispden, nu) = v1_qnu(1, :, ispden, nu) + &
           displ_red(1,ip,nu) * v1_atm(1,:,ispden,ip) - displ_red(2,ip,nu) * v1_atm(2,:,ispden,ip)
         v1_qnu(2, :, ispden, nu) = v1_qnu(2, :, ispden, nu) + &
           displ_red(2,ip,nu) * v1_atm(1,:,ispden,ip) + displ_red(1,ip,nu) * v1_atm(2,:,ispden,ip)
       else
         ! Gamma point. d(q) = d(-q)* --> d is real.
         v1_qnu(1, :, ispden, nu) = v1_qnu(1, :, ispden, nu) + displ_red(1,ip,nu) * v1_atm(1,:,ispden,ip)
       end if
     end do
   end do
 end do

end subroutine v1atm_to_vqnu

!!***

!----------------------------------------------------------------------

!!****f* m_gkk/find_mpw
!! NAME
!!  find_mpw
!!
!! FUNCTION
!!  Look at all k-points and spins to find the maximum
!!  number of plane waves.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_gkk
!!
!! NOTES
!!
!! CHILDREN
!!      get_kg
!!
!! SOURCE

subroutine find_mpw(mpw, kpts, nsppol, nkpt, gmet, ecut, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: mpw
 integer,intent(in) :: nsppol, nkpt
 integer,intent(in) :: comm
 real(dp),intent(in) :: ecut
!arrays
 real(dp),intent(in) :: kpts(3,nkpt)
 real(dp),intent(in) :: gmet(3,3)

!Local variables ------------------------------
!scalars
 integer :: my_rank, cnt, nproc, ierr
 integer :: ispin, ikpt
 integer :: my_mpw, onpw
 integer,allocatable :: gtmp(:,:)
 real(dp) :: kpt(3)

 my_rank = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 mpw = 0; cnt=0
 do ispin=1,nsppol
   do ikpt=1,nkpt
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle
     kpt = kpts(:,ikpt)
     call get_kg(kpt,1,ecut,gmet,onpw,gtmp)
     ABI_FREE(gtmp)
     mpw = max(mpw, onpw)
   end do
 end do
 my_mpw = mpw; call xmpi_max(my_mpw, mpw, comm, ierr)

end subroutine find_mpw
!!***
subroutine gkq_atm_to_gkq_nu(gkq_atm, atm_mass,displ_red, natom3, ntypat, phfreqs, typat, gkq_nu)
!Arguments ------------------------------------
!scalar
 integer,intent(in) :: natom3, ntypat

!arrays
 integer,intent(in) :: typat(ntypat)
 real(dp),intent(in) :: gkq_atm(2,natom3),displ_red(2,natom3,natom3),phfreqs(natom3),atm_mass(ntypat)
 real(dp),intent(out) :: gkq_nu(2,natom3)
!Local variables
 integer :: ipc,nu,natom,iat,idir
 real(dp) :: freq_prefs(natom3),mass_pref
 
 freq_prefs = sqrt(one/two/phfreqs)
 nu = 1
 natom = natom3/3

 gkq_nu = zero

 do nu=1,natom3
   ipc = 0
   do iat = 1,natom
     mass_pref = sqrt(one/atm_mass(typat(iat)))
     do idir = 1,3
       ipc = ipc + 1
       gkq_nu(1,nu) = gkq_nu(1,nu) + mass_pref*(displ_red(1,ipc,nu)*gkq_atm(1,ipc) - displ_red(2,ipc,nu)*gkq_atm(2,ipc))
       gkq_nu(2,nu) = gkq_nu(2,nu) + mass_pref*(displ_red(1,ipc,nu)*gkq_atm(2,ipc) + displ_red(2,ipc,nu)*gkq_atm(1,ipc))
     end do !idir
     gkq_nu(:,nu) = freq_prefs(nu)*gkq_nu(:,nu)
   end do !iat
 end do !nu
 

end subroutine gkq_atm_to_gkq_nu

subroutine gh1c_atm2nu(h1_kets, atm_mass,displ_red, natom3, ntypat, phfreqs, typat, h1_kets_nu) 
!Arguments ------------------------------------
!scalar
 integer,intent(in) :: natom3, ntypat

!arrays
 integer,intent(in) :: typat(ntypat)
 real(dp),intent(in) :: h1_kets(2,natom3),displ_red(2,natom3,natom3),phfreqs(natom3),atm_mass(ntypat)
 real(dp),intent(out) :: h1_kets_nu(2,natom3)
!Local variables
 integer :: ipc,nu,natom,iat,idir
 real(dp) :: freq_prefs(natom3),mass_pref
 
 freq_prefs = sqrt(one/two/phfreqs)
 nu = 1
 natom = natom3/3

 h1_kets_nu = zero

 do nu=1,natom3
   ipc = 0
   do iat = 1,natom
     mass_pref = sqrt(one/atm_mass(typat(iat)))
     do idir = 1,3
       ipc = ipc + 1
       h1_kets_nu(1,nu) = h1_kets_nu(1,nu) + mass_pref*(displ_red(1,ipc,nu)*h1_kets(1,ipc) - displ_red(2,ipc,nu)*h1_kets(2,ipc))
       h1_kets_nu(2,nu) = h1_kets_nu(2,nu) + mass_pref*(displ_red(1,ipc,nu)*h1_kets(2,ipc) + displ_red(2,ipc,nu)*h1_kets(1,ipc))
     end do !idir
     h1_kets_nu(:,nu) = freq_prefs(nu)*h1_kets_nu(:,nu)
   end do !iat
 end do !nu

end subroutine gh1c_atm2nu

subroutine mag_squared(mag_sq, vec_in)

 real(dp),intent(in) :: vec_in(2)
 real(dp),intent(out) :: mag_sq

mag_sq = vec_in(1)**2 + vec_in(2)**2

end subroutine mag_squared

! Copied from m_optic
subroutine get_opt_matel(pmat, wfkfile, ddkfile_1, ddkfile_2, ddkfile_3, comm)

 character(len=fnlen),intent(in) :: ddkfile_1,ddkfile_2,ddkfile_3
 character(len=fnlen),intent(inout) :: wfkfile
 
 complex(dpc),intent(out),allocatable :: pmat(:,:,:,:,:)
 
integer,intent(in) :: comm

 real(dp), ABI_CONTIGUOUS pointer :: outeig(:)
 integer,parameter :: formeig0=0,formeig1=1,master=0
 integer :: ii, fform,iomode0, ierr,bdtot0_index,bdtot_index,ikpt,isppol,mband,nkpt,nband1,nsppol,my_rank,varid,nsym
#ifdef HAVE_NETCDF
 integer :: ncid
#endif
 character(len=fnlen) :: infiles(0:3)
 character(len=500) :: msg
 integer :: iomode_ddk(1:3)
 integer,allocatable :: symrel(:,:,:),nband(:)
 logical :: use_ncevk(1:3)
 real(dp) :: gprimd(3,3),gprimd_trans(3,3), rprimd(3,3)
 real(dp),target,allocatable :: eigen11(:),eigen12(:),eigen13(:)
 real(dp),allocatable :: eigen0(:),eig0tmp(:),eigtmp(:),symcart(:,:,:)
 type(wfk_t) :: wfks(1:3),wfk0
 type(hdr_type) :: hdr_ddk(3), hdr

 my_rank = xmpi_comm_rank(comm)
 infiles = [wfkfile, ddkfile_1, ddkfile_2, ddkfile_3]

 
   ! Open GS wavefunction file
   ! Note: Cannot use MPI-IO here because of prtwf=3.
   ! If prtwf==3, the DDK file does not contain the wavefunctions but
   ! this info is not reported in the header and the offsets in wfk_compute_offsets
   ! are always computed assuming the presence of the cg
   call nctk_fort_or_ncfile(wfkfile, iomode0, msg)
   if (len_trim(msg) /= 0) MSG_ERROR(msg)
   if (iomode0 == IO_MODE_MPI) iomode0 = IO_MODE_FORTRAN
   call wfk_open_read(wfk0,wfkfile,formeig0,iomode0,get_unit(),xmpi_comm_self)
   ! Get header from the gs file
   call hdr_copy(wfk0%hdr, hdr)

   nkpt = hdr%nkpt
   rprimd(:,:)=hdr%rprimd(:,:)
   nsppol=hdr%nsppol
   nsym=hdr%nsym

 ! Identify the type of RF Wavefunction files
 use_ncevk = .False.
 do ii=1,3
   use_ncevk(ii) = endswith(infiles(ii), "_EVK.nc")
 end do

 ! Read ddk here from WFK files or from EVK.nc (only the header in the latter case)
 do ii=1,3

   call nctk_fort_or_ncfile(infiles(ii), iomode_ddk(ii), msg)
   if (len_trim(msg) /= 0) MSG_ERROR(msg)
   if (iomode_ddk(ii) == IO_MODE_MPI) iomode_ddk(ii) = IO_MODE_FORTRAN

   if (.not. use_ncevk(ii)) then
     call wfk_open_read(wfks(ii), infiles(ii), formeig1, iomode_ddk(ii), get_unit(), xmpi_comm_self)
     call hdr_copy(wfks(ii)%hdr,hdr_ddk(ii))
   else

#ifdef HAVE_NETCDF
    NCF_CHECK(nctk_open_read(ncid, infiles(ii), xmpi_comm_self))
    call hdr_ncread(hdr_ddk(ii),ncid, fform)
    ABI_CHECK(fform /= 0, sjoin("Error while reading:", infiles(ii)))

    NCF_CHECK(nf90_close(ncid))
#else
    MSG_ERROR("Netcdf not available!")
#endif

   end if
 end do
  

!   if(any(iomode_ddk(:)/=iomode0))then
!     write(msg, "(5a)")&
!&      ' The ground-state and ddk files should have the same format,',ch10,&
!&      ' either FORTRAN binary or NetCDF, which is not the case.',ch10,&
!&      ' Action : see input variable iomode.'
!     MSG_ERROR(msg)
!   endif

   ! Perform basic consistency tests for the GS WFK and the DDK files, e.g.
   ! k-points and their order, spins, number of bands could differ in the four files.
   ! Note indeed that we must have the same quantities in all the files.

   if (.not. use_ncevk(1)) then

     write(msg, "(12a)")ch10,&
&      ' Check the consistency of the wavefunction files (esp. k point and number of bands). ',ch10,&
&      ' Will compare, pairwise ( 1/2, 2/3, 3/4 ), the four following files :',ch10,&
&      trim(wfkfile),ch10,trim(infiles(1)),ch10,trim(infiles(2)),ch10,trim(infiles(3))
     call wrtout(std_out,msg,'COLL')

!DEBUG
!  stop
!ENDDEBUG

     if (hdr%compare(hdr_ddk(1)) /= 0) then
       write(msg, "(3a)")" Ground-state wavefunction file and ddkfile ",trim(infiles(1))," are not consistent. See above messages."
       MSG_ERROR(msg)
     end if
     do ii=1,2
       if (wfks(ii)%compare(wfks(ii+1)) /= 0) then
         write(msg, "(2(a,i0,a))")" ddkfile", ii," and ddkfile ",ii+1, ", are not consistent. See above messages"
         MSG_ERROR(msg)
       end if
     enddo
   endif

!DEBUG
!  stop
!ENDDEBUG

 ABI_ALLOCATE(nband,(nkpt*nsppol))

 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)

!Get mband, as the maximum value of nband(nkpt)
 mband=maxval(nband(:))
 do ii=1,nkpt
   if (nband(ii) /= mband) then
     MSG_ERROR("nband must be constant across kpts")
   end if
 end do

 ! Read the eigenvalues of ground-state and ddk files
 ABI_ALLOCATE(eigen0,(mband*nkpt*nsppol))
 ! MG: Do not understand why not [...,3]
 ABI_ALLOCATE(eigen11,(2*mband*mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen12,(2*mband*mband*nkpt*nsppol))
 ABI_ALLOCATE(eigen13,(2*mband*mband*nkpt*nsppol))



 if (my_rank == master) then
   ABI_ALLOCATE(eigtmp,(2*mband*mband))
   ABI_ALLOCATE(eig0tmp,(mband))

   do ii=1,3
     if (.not. use_ncevk(ii)) cycle
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_open_read(ncid, infiles(ii), xmpi_comm_self))
     varid = nctk_idname(ncid, "h1_matrix_elements")
     outeig => eigen11
     if (ii == 2) outeig => eigen12
     if (ii == 3) outeig => eigen13
     NCF_CHECK(nf90_get_var(ncid, varid, outeig, count=[2, mband, mband, nkpt, nsppol]))
     NCF_CHECK(nf90_close(ncid))
#else
     MSG_ERROR("Netcdf not available!")
#endif
   end do

   bdtot0_index=0 ; bdtot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt
       nband1=nband(ikpt+(isppol-1)*nkpt)
       eigtmp = zero
       eig0tmp = zero

       call wfk0%read_eigk(ikpt,isppol,xmpio_single,eig0tmp)
       eigen0(1+bdtot0_index:nband1+bdtot0_index)=eig0tmp(1:nband1)

       ! Read DDK matrix elements from WFK
       do ii=1,3
         if (.not. use_ncevk(ii)) then
           call wfks(ii)%read_eigk(ikpt, isppol, xmpio_single, eigtmp)
           if (ii == 1) eigen11(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
           if (ii == 2) eigen12(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
           if (ii == 3) eigen13(1+bdtot_index:2*nband1**2+bdtot_index)=eigtmp(1:2*nband1**2)
           !ABI_CHECK(wfks(ii)%nband(ikpt,isppol) == nband1, "ddk1 nband1")
         end if
       end do
       bdtot0_index=bdtot0_index+nband1
       bdtot_index=bdtot_index+2*nband1**2
     end do
   end do

   call wfk0%close()
   do ii=1,3
     if (.not. use_ncevk(ii)) call wfks(ii)%close()
   end do

   ABI_FREE(eigtmp)
   ABI_FREE(eig0tmp)
 end if ! master

 call xmpi_bcast(eigen0,master,comm,ierr)
 call xmpi_bcast(eigen11,master,comm,ierr)
 call xmpi_bcast(eigen12,master,comm,ierr)
 call xmpi_bcast(eigen13,master,comm,ierr)

 ABI_ALLOCATE(symcart,(3,3,nsym))
 ABI_ALLOCATE(symrel,(3,3,nsym))
 !YG: we need to transpose gprimd since matrinv give the transpose of the inverse!
 symrel(:,:,:) = hdr%symrel(:,:,:)
 gprimd_trans = transpose(gprimd)
 call sym2cart(gprimd_trans,nsym,rprimd,symrel,symcart)

 ABI_ALLOCATE(pmat,(mband,mband,nkpt,3,nsppol))
 call pmat2cart(eigen11,eigen12,eigen13,mband,nkpt,nsppol,pmat,rprimd)
 !TODO: pmat renorm?
 !call pmat_renorm(fermie, eigen0, mband, nkpt, nsppol, pmat, scissor)

 ABI_FREE(eigen0)
 ABI_FREE(eigen11)
 ABI_FREE(eigen12)
 ABI_FREE(eigen13)
 ABI_FREE(symcart)
 ABI_FREE(symrel)
end subroutine get_opt_matel

subroutine getgh1c_nopaw(gs_hamkq,rf_hamkq,dtset,psps,ebands_k,ik,ikq,kk,kq,idir,ipert,cryst,wfd_k,wfd_kq,mpw,mband,mband_kq,nspinor,gkq_atm,mpi_enreg,cwaveprj0,spin)

 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 integer,intent(in) :: idir,ipert,mpw,mband,mband_kq,nspinor
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands_k 
 type(wfd_t),intent(in) :: wfd_k,wfd_kq
 real(dp),allocatable,intent(out) :: gkq_atm(:,:,:)
 type(mpi_type),intent(in) :: mpi_enreg
 integer,intent(in) :: spin
 

 integer,parameter :: useylmgr1=0,optlocal=1,optnl=2,opt_gvnlx1=0,usevnl=0,berryopt0=0,tim_getgh1c=1
 real(dp) :: eig0nk, eshift,dotr,doti
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: gs1c(:,:),gvnlx1(:,:),grad_berry(:,:),bras(:,:,:),kets(:,:,:),h1_kets(:,:,:)
 integer :: kg_k(3,mpw), kg_kq(3,mpw)
 real(dp),intent(in) :: kk(3), kq(3)
 integer :: nkpg, nkpg1,istwf_k,istwf_kq,npw_k,npw_kq,ib1,ib2,sij_opt,ik,ikq
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:) !natom,nspinor*usecprj)

 ABI_MALLOC(ylm_k,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylm_kq,(mpw, psps%mpsang*psps%mpsang*psps%useylm))
 ABI_MALLOC(ylmgr_kq,(mpw, 3, psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
 ABI_MALLOC(gvnlx1, (2,usevnl))
 ABI_MALLOC(grad_berry, (2,nspinor*(berryopt0/4)))

 ABI_MALLOC(gkq_atm, (2,mband,mband))
 ABI_MALLOC(cwaveprj0, (2,0))

 
         ABI_MALLOC(bras, (2, mpw*nspinor, mband))
         ABI_MALLOC(kets, (2, mpw*nspinor, mband))
         ABI_MALLOC(h1_kets, (2, mpw*nspinor, mband))
         !TODO: ik check here, and ik loop
         ! Only do a subset a k-points


         ! Copy u_k(G)
         istwf_k = wfd_k%istwfk(ik); npw_k = wfd_k%npwarr(ik)
         ABI_CHECK(mpw >= npw_k, "mpw < npw_k")
         kg_k(:,1:npw_k) = wfd_k%kdata(ik)%kg_k
         do ib2=1,mband
           call wfd_k%copy_cg(ib2, ik, spin, kets(1,1,ib2))
         end do

         ! Copy u_kq(G)
         istwf_kq = wfd_kq%istwfk(ikq); npw_kq = wfd_kq%npwarr(ikq)
         ABI_CHECK(mpw >= npw_kq, "mpw < npw_kq")
         kg_kq(:,1:npw_kq) = wfd_kq%kdata(ikq)%kg_k
         do ib1=1,mband_kq
           call wfd_kq%copy_cg(ib1, ikq, spin, bras(1,1,ib1))
         end do


         ! if PAW, one has to solve a generalized eigenproblem
         ! Be careful here because I will need sij_opt==-1
         sij_opt = 0; if (psps%usepaw==1) sij_opt = 1
         ABI_MALLOC(gs1c, (2,npw_kq*nspinor*((sij_opt+1)/2)))

         

         ! GKA: Previous loop on 3*natom perturbations used to start here
         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
     call getgh1c_setup(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert,&    ! In
       cryst%natom,cryst%rmet,cryst%gprimd,cryst%gmet,istwf_k,&            ! In
       npw_k,npw_kq,useylmgr1,kg_k,ylm_k,kg_kq,ylm_kq,ylmgr_kq,&           ! In
       dkinpw,nkpg,nkpg1,kpg_k,kpg1_k,kinpw1,ffnlk,ffnl1,ph3d,ph3d1)       ! Out
     


         ! Calculate dvscf * psi_k, results stored in h1_kets on the k+q sphere.
         ! Compute H(1) applied to GS wavefunction Psi(0)
         do ib2=1,mband
           eig0nk = ebands_k%eig(ib2,ik,spin)
           ! Use scissor shift on 0-order eigenvalue
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0,kets(:,:,ib2),cwaveprj0,h1_kets(:,:,ib2),&
                        grad_berry,gs1c,gs_hamkq,gvnlx1,idir,ipert,eshift,mpi_enreg,optlocal,&
                        optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
         
         end do
         ! Calculate elphmat(j,i) = <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation.
         ! The array eig1_k contains:
         !
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
         ! <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
         do ib2=1,mband
           do ib1=1,mband_kq
             call dotprod_g(dotr,doti,istwf_kq,npw_kq*nspinor,2,bras(1,1,ib1),h1_kets(1,1,ib2),&
               mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

             gkq_atm(:, ib1, ib2) = [dotr, doti]
           end do
         end do

        ! call getgh1c_nopaw(gs_hamkq,rf_hamkq,dtset,psps,kk,kq,idir,ipert,cryst,istwf_k,npw_k,npw_kq,kg_k,kg_kq,mpw)
       ABI_FREE(ylm_k)
       ABI_FREE(ylm_kq)
       ABI_FREE(ylmgr_kq)
       ABI_FREE(kinpw1)
       ABI_FREE(kpg1_k)
       ABI_FREE(kpg_k)
       ABI_FREE(dkinpw)
       ABI_FREE(ffnlk)
       ABI_FREE(ffnl1)
       ABI_FREE(ph3d)
       ABI_SFREE(ph3d1)
       call pawcprj_free(cwaveprj0)
       ABI_FREE(cwaveprj0)

         ABI_FREE(gs1c)
         ABI_FREE(bras)
         ABI_FREE(kets)
         ABI_FREE(h1_kets)

end subroutine getgh1c_nopaw

end module m_gkk
!!***
