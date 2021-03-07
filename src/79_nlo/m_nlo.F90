!!****m* ABINIT/m_nlo
!! NAME
!!  m_nlo
!!
!! FUNCTION
!! DFT calculations of linear and nonlinear optical properties
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2020 ABINIT group (MVeithen,MB,LB)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_nlo

 use defs_basis
 use m_errors
 use m_abicore
 use m_build_info
 use m_xmpi
 use m_tetrahedron
 use m_htetra
 use m_crystal
 use m_ebands
 use m_dtset
 use m_wfd
 use m_krank
 use m_ddk
 use m_hdr
 use m_pawcprj
 use m_bz_mesh
 use m_wfk

 use defs_datatypes, only : ebands_t, pseudopotential_type
 use defs_abitypes, only : mpi_type
 use m_fstrings, only : itoa, sjoin, strcat
 use m_io_tools, only : get_unit, iomode_from_fname
 use m_numeric_tools, only: arth, simpson
 use m_symtk,    only : matr3inv
 use m_time, only : timein
 use m_pawtab,         only : pawtab_type
 use m_fftcore, only : get_kg
 use m_mpinfo, only : initmpi_seq, destroy_mpi_enreg
 use m_kpts, only    : tetra_from_kptrlatt, kpts_timrev_from_kptopt
 use m_sort, only : sort_int
 implicit none

 private
!!***

 public :: trans_rate_1pa
 public :: trans_rate_d2pa
 public :: kk_linperm
 public :: linopt_coefs
 public :: d2pa_coefs
 public :: wfd_read_all
!!***

contains

subroutine trans_rate_1pa(bcorr, cryst, ebands, nw, pol, scissor, trans_rate, wmin, wmax, wmesh, comm, &
&                         dtset, pawtab, psps, wfd, wfk_path)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bcorr, comm, nw
 real(dp),intent(in) :: scissor, wmin, wmax
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wfd_t),intent(in) :: wfd
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 character(len=*),intent(in) :: wfk_path

!arrays
 complex(dpc),intent(in) :: pol(3)
 real(dp),intent(out) :: trans_rate(nw)
 real(dp),intent(out) :: wmesh(nw)

!Local variables ------------------------------
!scalars
 character(len=500) :: wfk_path_fbz
 integer :: ierr, ibztet, cnt 
 integer :: nkbz, mpw, mpw_fbz, npw_k, istwf_k
 integer :: my_rank, nproc, tetra_num, num_sband=0, num_fband=0, mband
 integer :: ikpt, isppol, fband, sband, ihash, jtetra, iband, ifband, isband
 integer :: isym, itim, timrev, sym, nkpt_fullbz
 integer, parameter :: master=0
 logical :: compare_matrix_elements
 real(dp) :: step, ediff, population_diff, renorm_factor, max_occ, difference, total_diff=0
 type(htetra_t) :: htetra, htetra_fbz
 type(t_tetrahedron) :: tetra
 type(ddkop_t) :: ddkop, ddkop_fbz
 type(ebands_t) :: ebands_fbz
 type(kmesh_t) :: kmesh
 integer :: usecprj
 !WFK
 type(wfk_t) :: wfk
 type(MPI_type) :: mpi_enreg
 type(hdr_type) :: hdr
 type(krank_t) :: krank, krank_fbz
 type(hdr_type) :: hdr_fbz
 integer :: wfk_unit, sym_index, krank_prev
 integer,parameter :: formeig0 = 0
 type(wfd_t) :: wfd_fbz
 
!arrays
 character(len=500) :: msg, errmsg
 complex(dpc) :: prfs_tet(4)
 real(dp) :: p_fs(2), p_fs_fbz(2), vk(2,3,ebands%mband,ebands%mband)
 integer :: ik, ik_fbz, my_start, my_stop
 integer,allocatable :: sbands(:), fbands(:) 
 real(dp) :: klatt(3,3), rlatt(3,3), dweight(4,nw), tweight(4,nw), kpt(3), kpt_sym(3), weights(nw,2), weights_fbz(nw,2)
 real(dp),allocatable :: energy_fs(:), energy_fs_fbz(:), dweights(:,:,:,:),tweights(:,:), v_fbz(:,:,:,:,:,:)
 real(dp) :: integrand, v_bks(2,3)
 real(dp),allocatable :: kbz(:,:), cg_ket(:,:), cg_bra(:,:), ket_wfd(:,:), ket_wfk(:,:)
 real(dp),pointer :: eigen_fbz(:,:,:)
 integer,allocatable :: kg_k(:,:), iperm(:), rank(:), rank_fbz(:)
 complex(dpc),allocatable :: pmat(:,:,:,:,:), pol_rotated(:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:)
!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 
 
 ABI_CALLOC(fbands, (ebands%mband))
 ABI_CALLOC(sbands, (ebands%mband))

 mpw = maxval(wfd%npwarr)
 !WFK
 if (my_rank == 0) then
 wfk_unit = get_unit()
 call wfk_open_read(wfk, "copy_WFK", formeig0, iomode_from_fname(wfk_path), wfk_unit, xmpi_comm_self)
 call hdr_copy(wfk%hdr, hdr)
 end if

 ! Apply scissor shift to states entirely above the Fermi level
 !TODO: Should not update ebands with scissor shift because subsequent calls will be affected by the shift
 print *, "final bands"
   do iband = 1,ebands%mband
     if (all(ebands%eig(iband,:,:) > ebands%fermie)) then 
       num_fband = num_fband + 1
       fbands(num_fband) = iband
       !ebands%eig(iband,:,:) = ebands%eig(iband,:,:) + scissor 
     end if
   end do
 print *, fbands

 print *, "Starting bands"
 do iband = 1,ebands%mband
   if (any(fbands == iband)) cycle
   num_sband = num_sband + 1
   sbands(num_sband) = iband
 end do
 print *, sbands

 step = (wmax - wmin)/(nw - 1)
 wmesh = arth(wmin, step, nw)
 print *, "wmesh"
 print *, wmesh

 rlatt = ebands%kptrlatt
 call matr3inv(rlatt, klatt)

 ABI_ALLOCATE(energy_fs, (ebands%nkpt))

 
 print *, "nkpt"
 print *, ebands%nkpt
 ABI_ALLOCATE(dweights, (ebands%mband, ebands%mband,nw,ebands%nkpt))
 ABI_ALLOCATE(tweights, (nw,ebands%nkpt))
 ABI_MALLOC(cg_ket, (2, mpw*ebands%nspinor))
 ABI_MALLOC(cg_bra, (2, mpw*ebands%nspinor))
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (cryst%natom, ebands%nspinor*usecprj))

 weights = zero
 weights_fbz = zero

 trans_rate = zero
 call xmpi_split_work(ebands%nkpt,comm,my_start,my_stop)

 !WFD
 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)
 !WFK
 !call initmpi_seq(mpi_enreg)
 !ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, dtset%ngfft)
 call ebands_apply_scissors(ebands, scissor)

 timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 ABI_MALLOC(pol_rotated, (3, cryst%nsym*(timrev + 1)))
 ABI_MALLOC(iperm, (cryst%nsym*(timrev + 1)))
 iperm = [(isym, isym=1,cryst%nsym*(timrev + 1))]
 do itim = 0,timrev
   do isym=1,cryst%nsym
     pol_rotated(:,cryst%nsym*itim + isym) = (1-2*itim)*matmul(transpose(cryst%symrel_cart(:,:,isym)), pol)
   end do !isym
 end do !itim
 ABI_MALLOC(rank, (cryst%nsym*(timrev + 1)))


 htetra = tetra_from_kptrlatt(cryst, ebands%kptopt, ebands%kptrlatt, ebands%nshiftk, ebands%shiftk, ebands%nkpt, ebands%kptns, comm, msg, ierr)
 mband = ebands%mband
 call make_mesh(kmesh, cryst, 3, ebands%kptrlatt, ebands%nshiftk, ebands%shiftk)
 krank_fbz = krank_from_kptrlatt(kmesh%nbz, kmesh%bz, ebands%kptrlatt)
 compare_matrix_elements = .false.
 if (compare_matrix_elements) then

   print *, "compare = true"

   wfk_path_fbz = "GaAs_otfo_DS5_WFK"
   call wfk_read_eigenvalues(wfk_path_fbz, eigen_fbz, hdr_fbz, comm)
   ebands_fbz = ebands_from_hdr(hdr_fbz, maxval(hdr_fbz%nband), eigen_fbz)
   krank_fbz = krank_from_kptrlatt(ebands_fbz%nkpt, ebands_fbz%kptns, ebands_fbz%kptrlatt)
   call wfd_read_all(cryst, dtset, ebands_fbz, dtset%ngfft, pawtab, psps, wfd_fbz, wfk_path_fbz, comm)
   ABI_FREE(eigen_fbz)
   mpw_fbz = maxval(wfd_fbz%npwarr)
   ddkop_fbz = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw_fbz, wfd%ngfft)

   call ebands_apply_scissors(ebands_fbz, scissor)
   htetra_fbz = tetra_from_kptrlatt(cryst, ebands_fbz%kptopt, ebands_fbz%kptrlatt, ebands_fbz%nshiftk, ebands_fbz%shiftk, ebands_fbz%nkpt, ebands_fbz%kptns, comm, msg, ierr)

   ABI_MALLOC(v_fbz, (2,3,ebands_fbz%mband, ebands_fbz%mband, ebands_fbz%nkpt, ebands_fbz%nsppol))
   ABI_MALLOC(rank_fbz, (ebands_fbz%nkpt))
   ABI_MALLOC(energy_fs_fbz,(ebands_fbz%nkpt))

   do ik=1,ebands_fbz%nkpt
     kpt = ebands_fbz%kptns(:,ik)
     rank_fbz(ik) = krank_fbz%get_index(kpt)
     npw_k = wfd_fbz%npwarr(ik); istwf_k = wfd_fbz%istwfk(ik)
     call ddkop%setup_spin_kpoint(dtset, cryst, psps, isppol, ebands_fbz%kptns(:,ik), istwf_k, npw_k, wfd_fbz%kdata(ik)%kg_k)
     do isppol=1,ebands%nsppol
       do sband = 1,ebands%mband
         energy_fs_fbz = ebands_fbz%eig(fband,ik,isppol) - ebands_fbz%eig(sband,ik,isppol)
         call wfd_fbz%copy_cg(sband,ik,isppol,cg_ket)
         do fband = sband,ebands%mband
            call wfd_fbz%copy_cg(fband,ik,isppol,cg_bra)
            v_bks = ddkop%get_vnondiag(ebands_fbz%eig(fband,ik,isppol), istwf_k, npw_k, ebands_fbz%nspinor, cg_bra, cg_ket, cwaveprj0)
            v_fbz(:,:,sband,fband,ik,isppol) = v_bks
         end do !fband
       end do !sband
     end do ! isppol
   end do !ik
 end if

 krank = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt)

 max_occ = two/ebands%nspinor
 nkpt_fullbz = ebands%nshiftk*ebands%kptrlatt(1,1)*ebands%kptrlatt(2,2)*ebands%kptrlatt(3,3)
 call identk(ebands%kptns, ebands%nkpt, nkpt_fullbz, crystal%nsym, timrev + 1, crystal%symrec, crystal%symafm, kbz, ktab, ktabi, ktabo, nkbz)

 if (ebands%kptopt == 3) then !all k points are inequivalent
    do isppol = 1,ebands%nsppol   
      do ik = my_start, my_stop 
        npw_k = wfd%npwarr(ik); istwf_k = wfd%istwfk(ik)
        call ddkop%setup_spin_kpoint(dtset, cryst, psps, isppol, ebands%kptns(:,ik), istwf_k, npw_k, wfd%kdata(ik)%kg_k)
        do sband = 1,mband
          call wfd%copy_cg(sband,ik,isppol,cg_ket)
          do fband = sband,mband
            population_diff = ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)
            if (abs(population_diff) < 1.0d-12) cycle
            energy_fs = ebands%eig(fband,:,isppol) - ebands%eig(sband,:,isppol)
            call htetra%get_onewk_wvals(ik, 0, nw, wmesh, max_occ, ebands%nkpt, energy_fs, weights)
            if (all(abs(weights(:,1)) .lt. 1.0d-12)) cycle
            call wfd%copy_cg(fband,ik,isppol,cg_bra)
            v_bks = ddkop%get_vnondiag(ebands%eig(sband,ik,isppol), istwf_k, npw_k, ebands%nspinor, cg_bra, cg_ket, cwaveprj0)
            p_fs = pol(1)*v_bks(:,1) + pol(2)*v_bks(:,2) &
&              + pol(3)*v_bks(:,3)
            integrand =(p_fs(1)**2 + p_fs(2)**2)*population_diff
            trans_rate = trans_rate + integrand*weights(:,1)
          end do !fband
        end do !sband
      end do !ik
    end do !isppol
    call xmpi_sum(trans_rate, comm, ierr)
    if (ierr /= 0) MSG_ERROR("Error in xmpi sum for transition rate")
    trans_rate = trans_rate/ebands%nkpt

 else if (ebands%kptopt == 1) then
   if (my_rank == master) then
     call wrtout(std_out, "Generating tetrahedra in the irreducible Brillouin zone", "COLL")  
   end if
     do isppol = 1,ebands%nsppol   
       do ik = 1,ebands%nkpt
         !vk(2,3,mband,mband) are matrix elements between bands at given ik (in irreducible wedge) 
         vk = zero
         npw_k = wfd%npwarr(ik); istwf_k = wfd%istwfk(ik)
         ! ik loop outside so this is only called once per k point
         call ddkop%setup_spin_kpoint(dtset, cryst, psps, isppol, ebands%kptns(:,ik), istwf_k, npw_k, wfd%kdata(ik)%kg_k)
         ! Precompute all needed matrix elements in IBZ
         do sband = 1,mband
           call wfd%copy_cg(sband,ik,isppol,cg_ket)
           call ddkop%apply(ebands%eig(sband,ik,isppol), npw_k, ebands%nspinor, cg_ket, cwaveprj0)
           do fband = sband,mband
             ! Don't bother computing if they have the same population
             if (abs(ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)) .lt. 1.0d-12) cycle
             call wfd%copy_cg(fband,ik,isppol,cg_bra)
             v_bks = ddkop%get_braket(ebands%eig(sband,ik,isppol),istwf_k, npw_k, ebands%nspinor, cg_bra, mode="cart")
             ! renorm? Just seems to make errors larger.
             renorm_factor = (ebands%eig(fband,ik,isppol) - ebands%eig(sband,ik,isppol))/(ebands%eig(fband,ik,isppol) - ebands%eig(sband,ik,isppol) - scissor)
             renorm_factor = one
             vk(:,:,sband,fband) = v_bks*(renorm_factor)
           end do !fband
         end do !sband
         
         ! Loop to get ranks of all k points in star (with duplicates)
         sym_index = 1
         do itim = 0,timrev
           do isym = 1,cryst%nsym
             kpt_sym = matmul(cryst%symrec(:,:,isym), ebands%kptns(:,ik))
             rank(sym_index) = krank_fbz%get_index(kpt_sym)
             sym_index = sym_index + 1
           end do !isym
         end do !itim
         ! Sort the array containing ranks so that we can easily pick the first instance of each
         call sort_int((timrev + 1)*cryst%nsym, rank, iperm)


         ! Loop over bands on outside allows us to only compute tetrahedron weights once
         do sband = 1,mband-2
           ! fband loop starts as sband to avoid double counting transitions
           do fband = sband,mband-2
             population_diff = ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)
             ! Population filter again
             if (abs(population_diff) .lt. 1.0d-12) cycle
             energy_fs = ebands%eig(fband,:,isppol) - ebands%eig(sband,:,isppol)
             ! get weights
             call htetra%get_onewk_wvals(ik, 0, nw, wmesh, max_occ, ebands%nkpt, energy_fs, weights)
             ! filter based on weights
             if (all(abs(weights(:,1)) .lt. 1.0d-12)) cycle
             ! initialize so that loop runs the first time
             krank_prev = -2
             ! Loop over all isym,itim
             do isym = 1,(timrev + 1)*cryst%nsym
               ! multiple sym ops can give the same k point, skip if this is repeat (remember the array is sorted by krank)
               if (rank(isym) == krank_prev) cycle
               krank_prev = rank(isym)
               ! iperm sorted so that iperm(isym) is the index of the symmetry operation
               sym = iperm(isym)
               p_fs = pol_rotated(1, isym)*vk(:,1,sband,fband) + pol_rotated(2,isym)*vk(:,2,sband,fband) &
&                 + pol_rotated(3,isym)*vk(:,3,sband,fband)
               if (compare_matrix_elements) then
                 print *, "sym applied"
                 print *, sym
                 print *, "ik_fbz"
                 print *, krank_prev
                 print *, rank
                 print *, "ibz k"
                 print *, ebands%kptns(:,ik)
                 print *, "fbz k"
                 print *, ebands_fbz%kptns(:,krank_prev)
                 !print *, "eig"
                 energy_fs_fbz =  ebands_fbz%eig(fband, :, isppol) - ebands_fbz%eig(sband,:,isppol)
                 print *," ediff fbz"
                 print *,ebands_fbz%eig(fband, krank_prev, isppol) - ebands_fbz%eig(sband,krank_prev,isppol)

                 print *," ediff ibz"
                 print *,ebands%eig(fband, ik, isppol) - ebands%eig(sband,ik,isppol)
                 call htetra_fbz%get_onewk_wvals(krank_prev, 0, nw, wmesh, max_occ, ebands_fbz%nkpt, energy_fs_fbz, weights_fbz)

                 p_fs_fbz = pol(1)*v_fbz(:,1,sband, fband, krank_prev, isppol) + pol(2)*v_fbz(:,2,sband, fband, krank_prev, isppol) &
                         & + pol(3)*v_fbz(:,3,sband, fband, krank_prev, isppol)
                 difference = (p_fs(1)**2 + p_fs(2)**2) - (p_fs_fbz(1)**2 + p_fs_fbz(2)**2)
                 total_diff = total_diff + difference
                 print *, "ibz weights"
                 print *, weights(:,1)
                 print *, "fbz weights"
                 print *, weights_fbz(:,1)
                 if (abs(difference) .gt. 1.0d-4) then
                   print *, "ik ibz"
                   print *, ik
                   print *, "sband"
                   print *, sband
                   print *, "fband"
                   print *, fband
                   print *, "kpt"
                   print *, "diff"
                   print *, difference
                 end if
               end if
               integrand =(p_fs(1)**2 + p_fs(2)**2)*population_diff
               trans_rate = trans_rate + integrand*weights(:,1)
             end do !isym
           end do !fband
         end do !sband
       end do !ik
     end do !isppol
     call xmpi_sum(trans_rate, comm, ierr)
     if (ierr /= 0) MSG_ERROR("Error in xmpi sum for transition rate")
     trans_rate = trans_rate/nkpt_fullbz
 end if 
 print *, "Total diff"
 print *, total_diff

 call htetra%free
 ABI_DEALLOCATE(energy_fs)
 ABI_DEALLOCATE(sbands)
 ABI_DEALLOCATE(fbands)
 ABI_DEALLOCATE(dweights)
 ABI_DEALLOCATE(tweights)
 ABI_DEALLOCATE(cg_ket)
 ABI_DEALLOCATE(cg_bra)
 ABI_FREE(pol_rotated)
 ABI_FREE(iperm)
 ABI_FREE(rank)
 if (compare_matrix_elements) then 
   ABI_FREE(v_fbz)
   ABI_FREE(rank_fbz)
   ABI_FREE(energy_fs_fbz)
 end if

end subroutine trans_rate_1pa

!TODO: Fix this function, currently doesn't work
subroutine trans_rate_d2pa(bcorr, cryst, ebands, nw, pmat, pol_1, pol_2, scissor, total_trans_rate, wmin, wmax, wmesh, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm, bcorr
 real(dp),intent(in) :: scissor, wmin, wmax
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands

!arrays
 complex(dpc),intent(in) :: pmat(:,:,:,:,:), pol_1(3), pol_2(3)
 integer,intent(in) :: nw
 real(dp),intent(out) :: total_trans_rate(nw)
 real(dp),intent(out) :: wmesh(nw)

!Local variables ------------------------------
!scalars
 integer :: ierr
 integer :: my_rank, nproc, tetra_num, num_fband=0, num_sband=0
 integer :: fband, sband, ik
 integer :: ikpt, isppol, iband, ifband, isband, iw
 integer, parameter :: master=0
 real(dp) :: step, detuning, trans_rate
 type(t_tetrahedron) :: tetra
 
!arrays
 character(len=500) :: msg, errmsg
 complex(dpc) :: integrand_fs, int_summand
 integer :: fbands(20), sbands(20)
 integer,allocatable :: bz2ibz_indexes(:) 
 real(dp) :: klatt(3,3), rlatt(3,3)
 real(dp) :: energy_is
 complex(dpc) :: p_is_1, p_is_2, p_fi_1, p_fi_2
 real(dp) :: integrand, detuning_k
 real(dp),allocatable :: energy_fs(:), dweights(:,:), tweights(:,:)

!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 
 
 step = (wmax - wmin)/(nw - 1)
 wmesh = arth(wmin, step, nw)


!TODO: Make the band selection match that in the 1PA function
print *, "Starting bands"
do sband = 1,ebands%mband
  if (all(ebands%occ(sband,:,:) == 0)) cycle
  num_sband = num_sband + 1
  print *, sband
  sbands(num_sband) = sband
end do

print *, "Final bands"
do fband = 1,ebands%mband
  if (all(ebands%occ(fband,:,:) /= 0)) cycle
  num_fband = num_fband + 1
  print *, fband
  fbands(num_fband) = fband
end do

 do isppol = 1,ebands%nsppol
   do iband = 1,ebands%mband
     if (all(ebands%eig(iband,:,isppol) > ebands%fermie)) then 
       write(std_out, '(a,i2)') "Scissor shift applied to band", iband
       ebands%eig(iband,:,isppol) = ebands%eig(iband,:,isppol) + scissor 
     end if
   end do
 end do

 if (ebands%kptopt == 3) then !all k points are inequivalent
   if (my_rank == master) then
     call wrtout(std_out, "Generating inequivalent tetrahedra across the full Brillouin zone", "COLL")  
   end if
     
     rlatt = ebands%kptrlatt
     call matr3inv(rlatt, klatt)
     ABI_ALLOCATE(bz2ibz_indexes, (ebands%nkpt))
     ABI_ALLOCATE(dweights, (nw, ebands%nkpt))
     ABI_ALLOCATE(tweights, (nw, ebands%nkpt))
     ABI_ALLOCATE(energy_fs, (ebands%nkpt))
     bz2ibz_indexes = [(ikpt, ikpt=1,ebands%nkpt)]
     
     
     
     call init_tetra(bz2ibz_indexes, cryst%gprimd, klatt, ebands%kptns, ebands%nkpt, tetra, ierr, errmsg, comm)
     if (ierr /= 0) MSG_ERROR(errmsg)

     total_trans_rate = zero
     do isppol = 1,ebands%nsppol
       do isband = 1,num_sband
         sband = sbands(isband)
         do ifband = 1,num_fband
           fband = fbands(ifband)
           energy_fs = ebands%eig(fband, :, isppol) - ebands%eig(sband, :, isppol)
           call tetra_blochl_weights(tetra, energy_fs, two*wmin, two*wmax, one, nw, ebands%nkpt, bcorr, tweights, dweights, comm) 
           do ik=1,ebands%nkpt
               do iw = 1,nw
                 if (dweights(iw,ik) == zero) cycle
               	 int_summand = zero 
                 do iband = 1,ebands%mband
                   p_is_1 = pol_1(1)*pmat(iband,sband,ik,1,isppol) + pol_1(2)*pmat(iband,sband,ik,2,isppol) + pol_1(3)*pmat(iband,sband,ik,3,isppol)
                   p_is_2 = pol_2(1)*pmat(iband,sband,ik,1,isppol) + pol_2(2)*pmat(iband,sband,ik,2,isppol) + pol_2(3)*pmat(iband,sband,ik,3,isppol)
                   p_fi_1 = pol_1(1)*pmat(fband,iband,ik,1,isppol) + pol_1(2)*pmat(fband,iband,ik,2,isppol) + pol_1(3)*pmat(fband,iband,ik,3,isppol)
                   p_fi_2 = pol_2(1)*pmat(fband,iband,ik,1,isppol) + pol_2(2)*pmat(fband,iband,ik,2,isppol) + pol_2(3)*pmat(fband,iband,ik,3,isppol)

                   detuning_k = ebands%eig(iband, ik, isppol) - ebands%eig(sband, ik, isppol) &
&                               - wmesh(iw)
                   int_summand = int_summand + (p_is_1*p_fi_2 + p_is_2*p_fi_1)/(detuning_k + (0, 0.01))
                 end do !iband
                 integrand = abs(int_summand)**2*(ebands%occ(sband,ik,isppol) &
&                                                          - ebands%occ(fband,ik,isppol))
                 total_trans_rate(iw) = total_trans_rate(iw) + integrand*dweights(iw,ik)
               end do !iw
           end do !ik
         end do !ifband
       end do !isband
     end do !isppol

     call xmpi_sum(total_trans_rate, comm, ierr)
     total_trans_rate = total_trans_rate/four
     if (ierr /= 0) MSG_ERROR("Error in xmpi sum for transition rate")
 call destroy_tetra(tetra)
 ABI_DEALLOCATE(bz2ibz_indexes)
 ABI_DEALLOCATE(energy_fs)
 ABI_DEALLOCATE(dweights)
 ABI_DEALLOCATE(tweights)
 end if


end subroutine trans_rate_d2pa

subroutine kk_linperm(eps_imag, len_abs, eps_real, wmesh)

 integer,intent(in) :: len_abs

 real(dp),intent(in) :: eps_imag(len_abs)
 real(dp),intent(out) :: eps_real(len_abs)

 real(dp),intent(in) :: wmesh(len_abs)

 integer :: iw, iwp
 real(dp) :: wstep, w, wp
 real(dp),parameter :: c = 137_dp
 real(dp) :: integrand(len_abs)
 
 wstep = wmesh(2) - wmesh(1)
 
 do iw = 1,len_abs
   w = wmesh(iw)
   do iwp = 1,len_abs
     if (iwp == iw) cycle
     integrand(iwp) = wmesh(iwp)*eps_imag(iwp)/(wmesh(iwp)**2 - w**2)
   end do
   eps_real(iw) = simpson(wstep,integrand)
 end do

 eps_real = one + two/(pi)*eps_real

 end subroutine kk_linperm


 subroutine linopt_coefs(alpha_units_per_cm, bcorr, cryst, ebands, fname_root, ngfft, nw, pawtab, pol, scissor, wmax, comm, dtset, psps, wfk0_path)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bcorr, comm, nw
 real(dp),intent(in) :: scissor, wmax
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands
 integer,intent(in) :: alpha_units_per_cm
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 character(len=*),intent(in) :: wfk0_path
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!arrays
 complex(dpc),intent(in) :: pol(3)
 character(len=*),intent(in) :: fname_root
 integer,intent(in) :: ngfft(18)

!Local variables ------------------------------
!scalars
 integer :: my_rank, nproc 
 integer :: iw, funit
 integer, parameter :: master=0
 real(dp),parameter :: wminzero=zero
 !WFD
 type(wfd_t) :: wfd
 
!arrays
 character(len=500) :: msg, errmsg, fname
 real(dp) :: trans_rate(nw),eps_mag(nw)
 real(dp) :: index_alpha(2,nw), eps(2,nw), wmesh(nw)
 complex(dpc),allocatable :: pmat(:,:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 integer,allocatable :: wfd_istwfk(:),nband(:,:)

!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 

 ABI_MALLOC(nband, (ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, ebands%nkpt, ebands%nsppol))
 !WFD
 ! Read in all wave functions (memory intensive)
 bks_mask = .True.; keep_ur = .False.
 !WFK
 !bks_mask = .False.; keep_ur = .False.
 ABI_MALLOC(wfd_istwfk, (ebands%nkpt))
 wfd_istwfk = 1
 nband = dtset%mband

 !WFD
 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, dtset%mband, nband, ebands%nkpt, ebands%nsppol, bks_mask,&
               dtset%nspden, ebands%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))
 !WFK
 ! Either init with bks mas False (don't allocate any memory but still use wfd info) or don't initialize at all

 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(nband)

 call trans_rate_1pa(bcorr, cryst, ebands, nw, pol, scissor, trans_rate, wminzero, wmax, wmesh, comm, dtset, pawtab, psps, wfd, wfk0_path)
 eps(2,:) = trans_rate*(1_dp/wmesh**2)*one/(two_pi)
 eps(2,:) = trans_rate*(1_dp/wmesh**2)*one/(four*pi)
 eps(2,1) = 0
 call kk_linperm(eps(2,:), nw, eps(1,:), wmesh)
 eps_mag = sqrt(eps(1,:)**2 + eps(2,:)**2)

 index_alpha(1,:) = sqrt(one/two*(eps_mag + eps(1,:)))
 index_alpha(2,:) = sqrt(one/two*(eps_mag - eps(1,:)))*two*wmesh/137_dp

 if (alpha_units_per_cm==1) index_alpha(2,:) = index_alpha(2,:)/5.29d-9

 if (my_rank == master) then
   funit = get_unit()
   fname = trim(fname_root)//'_linopt.out'
   open(funit, file = fname)
   do iw = 1,nw
     write(funit, *) wmesh(iw), index_alpha(2,iw), index_alpha(1,iw), eps(1,iw), eps(2,iw)
   end do
   close(funit)
 end if

 !WFD (WFK?)
 !call wfd%free()

 end subroutine linopt_coefs

 subroutine d2pa_coefs(alpha2_d, bcorr, cryst, ebands, fname_root, nw, pmat, pol_1, pol_2, scissor, wmin, wmax, wmesh, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm, bcorr
 real(dp),intent(in) :: scissor, wmin, wmax
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands

!arrays
 complex(dpc),intent(in) :: pmat(:,:,:,:,:), pol_1(3), pol_2(3)
 integer,intent(in) :: nw
 real(dp),intent(out) :: alpha2_d(nw)
 real(dp),intent(out) :: wmesh(nw)
 character(len=*),intent(in) :: fname_root

!Local variables ------------------------------
!scalars
 integer :: ierr, funit, iw
 integer :: my_rank, nproc 
 integer, parameter :: master=0
 real(dp) :: step
 
!arrays
 character(len=500) :: msg, errmsg, fname
 real(dp) :: trans_rate(nw)

!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 
 
 step = (wmax - wmin)/(nw - 1)
 wmesh = arth(wmin, step, nw)

 call trans_rate_d2pa(bcorr, cryst, ebands, nw, pmat, pol_1, pol_2, scissor, trans_rate, wmin, wmax, wmesh, comm)
 alpha2_d = 4_dp/(1**2*137**2*wmesh**3)*trans_rate
 alpha2_d(1) = 0
 print *, "alpha 2"
 print *, alpha2_d

 if (my_rank == master) then
   fname = trim(fname_root)//'_d2pa.out'
   funit = get_unit()
   open(funit, file = fname)
   do iw = 1,nw
     write(funit, *) wmesh(iw), alpha2_d(iw)
   end do
   close(funit)
 end if

 end subroutine d2pa_coefs

 subroutine wfd_read_all(cryst, dtset, ebands, ngfft, pawtab, psps, wfd, wfk_path, comm)
!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wfd_t),intent(out) :: wfd
 character(len=*),intent(in) :: wfk_path
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 integer,intent(in) :: comm

!arrays
 integer,intent(in) :: ngfft(18)

!Local variables ------------------------------
!scalars
 
!arrays
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 integer,allocatable :: wfd_istwfk(:),nband(:,:)


 ABI_MALLOC(nband, (ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(bks_mask, (ebands%mband, ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(keep_ur, (ebands%mband, ebands%nkpt, ebands%nsppol))
 ! Read in all wave functions (memory intensive)
 bks_mask = .True.; keep_ur = .False.
 ABI_MALLOC(wfd_istwfk, (ebands%nkpt))
 wfd_istwfk = 1
 nband = ebands%mband

 !WFD
 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, ebands%mband, nband, ebands%nkpt, ebands%nsppol, bks_mask,&
               dtset%nspden, ebands%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)
 call wfd%read_wfk(wfk_path, iomode_from_fname(wfk_path))
 !WFK
 ! Either init with bks mas False (don't allocate any memory but still use wfd info) or don't initialize at all

 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(nband)

 end subroutine wfd_read_all


end module m_nlo


