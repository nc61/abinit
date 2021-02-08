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
 use m_ddk
 use m_hdr
 use m_pawcprj
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
 implicit none

 private
!!***

 public :: trans_rate_1pa
 public :: trans_rate_d2pa
 public :: kk_linperm
 public :: linopt_coefs
 public :: d2pa_coefs
!!***

contains

subroutine trans_rate_1pa(bcorr, cryst, ebands, nw, polarization, scissor, trans_rate, wmin, wmax, wmesh, comm, &
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
 complex(dpc),intent(in) :: polarization(3)
 real(dp),intent(out) :: trans_rate(nw)
 real(dp),intent(out) :: wmesh(nw)

!Local variables ------------------------------
!scalars
 integer :: ierr, ibztet, cnt 
 integer :: nkbz, mpw, npw_k, istwf_k
 integer :: my_rank, nproc, tetra_num, num_sband=0, num_fband=0
 integer :: ikpt, isppol, fband, sband, ihash, jtetra, iband, ifband, isband
 integer, parameter :: master=0
 real(dp) :: step, ediff
 type(htetra_t) :: htetra
 type(t_tetrahedron) :: tetra
 type(ddkop_t) :: ddkop
 integer :: usecprj
 !WFK
 type(wfk_t) :: wfk
 type(MPI_type) :: mpi_enreg
 type(hdr_type) :: hdr
 integer :: wfk_unit
 integer,parameter :: formeig0 = 0
 
!arrays
 character(len=500) :: msg, errmsg
 complex(dpc) :: prfs_tet(4)
 real(dp) :: p_fs(2)
 integer :: isym(4), itim(4)
 integer :: ik, my_start, my_stop
 integer,allocatable :: bz2ibz_indexes(:),sbands(:), fbands(:)
 real(dp) :: klatt(3,3), rlatt(3,3), dweight(4,nw), tweight(4,nw), kk(3), weights(nw,2)
 real(dp),allocatable :: energy_fs(:),dweights(:,:,:,:),tweights(:,:)
 real(dp) :: integrand
 real(dp),allocatable :: kbz(:,:), cg_ket(:,:), cg_bra(:,:), ket_wfd(:,:), ket_wfk(:,:)
 integer,allocatable :: kg_k(:,:)
 complex(dpc),allocatable :: pmat(:,:,:,:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:)
 real(dp) :: vk(2,3)
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
 !mpw = maxval(hdr%npwarr)

 ABI_ALLOCATE(ket_wfd,(2, mpw))
 ABI_ALLOCATE(ket_wfk,(2, mpw))

 print *, "mpw"
 print *, mpw

 print *, "wfd ket"
 call wfd%copy_cg(4,1,1,ket_wfd)
 print *, ket_wfd
 
 if (my_rank == 0) then
 print *, "wfk ket"
 call wfk%read_bks(4, 1, 1, xmpio_single, cg_bks=ket_wfk)
 print *, ket_wfk
 end if

 print *, "subtraction"
 print *, ket_wfd - ket_wfk

 print *, "npwarr"
 print *, wfd%npwarr(1)

 ABI_FREE(ket_wfd)
 ABI_FREE(ket_wfk)

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

 ABI_ALLOCATE(bz2ibz_indexes, (ebands%nkpt))
 ABI_ALLOCATE(energy_fs, (ebands%nkpt))

 bz2ibz_indexes = [(ikpt, ikpt=1,ebands%nkpt)]
 
 ABI_ALLOCATE(dweights, (ebands%mband, ebands%mband,nw,ebands%nkpt))
 ABI_ALLOCATE(tweights, (nw,ebands%nkpt))
 ABI_MALLOC(cg_ket, (2, mpw*ebands%nspinor))
 ABI_MALLOC(cg_bra, (2, mpw*ebands%nspinor))
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (cryst%natom, ebands%nspinor*usecprj))

 trans_rate = zero
 call xmpi_split_work(ebands%nkpt,comm,my_start,my_stop)

 !WFD
 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)
 !WFK
 !call initmpi_seq(mpi_enreg)
 !ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, dtset%ngfft)


 if (ebands%kptopt == 3) then !all k points are inequivalent
   if (my_rank == master) then
     call wrtout(std_out, "Generating inequivalent tetrahedra across the full Brillouin zone", "COLL")  
   end if
     call htetra_init(htetra, bz2ibz_indexes, cryst%gprimd, klatt, ebands%kptns, & 
&                     ebands%nkpt, ebands%kptns, ebands%nkpt, ierr, errmsg, comm, 2)
     ! precompute weights so that k can be on outer loop so that ddkop%setup_spin_point is only run once per k
    !do isppol = 1,ebands%nsppol   
    !  do sband = 1,num_sband
    !    isband = sbands(isband)
    !    do ifband = 1,num_fband
    !      fband = fbands(ifband)
    !      !energy_fs = ebands%eig(fband,:,isppol) - ebands%eig(sband,:,isppol) + scissor
    !      !TODO: Figure out why this is running more slowly with MPI than single core.
    !      !call htetra%blochl_weights(energy_fs, wmin, wmax, one, nw, ebands%nkpt, bcorr, tweights(:,:), dweights(sband,fband,:,:), comm) 
    !    end do !fband
    !  end do !sband
    !end do !isppol
     if (ierr /= 0) MSG_ERROR(errmsg)
     do isppol = 1,ebands%nsppol   
       do ik = my_start, my_stop 
         !wfk%read_band_block([1,10], ik, isppol, xmpio_single, kg_k=kg_k, cg_k)
         !WFD
         npw_k = wfd%npwarr(ik); istwf_k = wfd%istwfk(ik)
         call ddkop%setup_spin_kpoint(dtset, cryst, psps, isppol, ebands%kptns(:,ik), istwf_k, npw_k, wfd%kdata(ik)%kg_k)
         !WFK
         !istwf_k = 1
         !call get_kg(ebands%kptns(:,ik), istwf_k, dtset%ecut, cryst%gmet, npw_k, kg_k)
         !call ddkop%setup_spin_kpoint(dtset, cryst, psps, isppol, ebands%kptns(:,ik), istwf_k, npw_k, kg_k)
         do sband = 1,num_sband
           isband = sbands(isband)
           call wfd%copy_cg(sband,ik,isppol,cg_ket)
           do ifband = 1,num_fband
             fband = fbands(ifband)

             !WFD
             call wfd%copy_cg(fband,ik,isppol,cg_bra)
             !WFK
             !call wfk%read_bks(sband, ik, isppol, xmpio_single, cg_bks=cg_ket)
             !call wfk%read_bks(fband, ik, isppol, xmpio_single, cg_bks=cg_bra)

             vk = ddkop%get_vnondiag(ebands%eig(fband,ik,isppol), istwf_k, npw_k, ebands%nspinor, cg_bra, cg_ket, cwaveprj0)
             ediff = ebands%eig(fband,ik,isppol) - ebands%eig(sband,ik,isppol) 
             !vk = vk*(ediff+scissor)/ediff 

             p_fs = polarization(1)*vk(:,1) + polarization(2)*vk(:,2) &
&               + polarization(3)*vk(:,3)
             integrand =(p_fs(1)**2 + p_fs(2)**2)*(ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol))
             energy_fs = ebands%eig(fband,:,isppol) - ebands%eig(sband,:,isppol) + scissor
             call htetra%get_onewk_wvals(ik, 1, nw, wmesh, one, ebands%nkpt, energy_fs, weights)
             !print *, "weights"
             !print *, weights(:,1)
             trans_rate = trans_rate + integrand*weights(:,1)/(ebands%nkpt)
           end do !fband
         end do !sband
       end do !ik
     end do !isppol
     call xmpi_sum(trans_rate, comm, ierr)
     if (ierr /= 0) MSG_ERROR("Error in xmpi sum for transition rate")

 !TODO: Fix this or delete it. This section tries to calculate properties using only the irreducible wedge and symmetry properties (rotating matrix elements). Just need to account for phase correctly
 else if (ebands%kptopt == 1) then
   if (my_rank == master) then
     call wrtout(std_out, "Generating tetrahedra in the irreducible Brillouin zone", "COLL")  
   end if
 end if 

 call htetra%free
 ABI_DEALLOCATE(bz2ibz_indexes)
 ABI_DEALLOCATE(energy_fs)
 ABI_DEALLOCATE(sbands)
 ABI_DEALLOCATE(fbands)
 ABI_DEALLOCATE(dweights)
 ABI_DEALLOCATE(tweights)
 ABI_DEALLOCATE(cg_ket)
 ABI_DEALLOCATE(cg_bra)
 !WFK
 !ABI_DEALLOCATE(kg_k)
 !call destroy_mpi_enreg(mpi_enreg)
 !close(wfk_unit)

end subroutine trans_rate_1pa

!TODO: Fix this function, currently doesn't work
subroutine trans_rate_d2pa(bcorr, cryst, ebands, nw, pmat, polarization_1, polarization_2, scissor, total_trans_rate, wmin, wmax, wmesh, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm, bcorr
 real(dp),intent(in) :: scissor, wmin, wmax
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands

!arrays
 complex(dpc),intent(in) :: pmat(:,:,:,:,:), polarization_1(3), polarization_2(3)
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
                   p_is_1 = polarization_1(1)*pmat(iband,sband,ik,1,isppol) + polarization_1(2)*pmat(iband,sband,ik,2,isppol) + polarization_1(3)*pmat(iband,sband,ik,3,isppol)
                   p_is_2 = polarization_2(1)*pmat(iband,sband,ik,1,isppol) + polarization_2(2)*pmat(iband,sband,ik,2,isppol) + polarization_2(3)*pmat(iband,sband,ik,3,isppol)
                   p_fi_1 = polarization_1(1)*pmat(fband,iband,ik,1,isppol) + polarization_1(2)*pmat(fband,iband,ik,2,isppol) + polarization_1(3)*pmat(fband,iband,ik,3,isppol)
                   p_fi_2 = polarization_2(1)*pmat(fband,iband,ik,1,isppol) + polarization_2(2)*pmat(fband,iband,ik,2,isppol) + polarization_2(3)*pmat(fband,iband,ik,3,isppol)

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


 subroutine linopt_coefs(alpha_units_per_cm, bcorr, cryst, ebands, fname_root, ngfft, nw, pawtab, polarization, scissor, wmax, comm, dtset, psps, wfk0_path)

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
 complex(dpc),intent(in) :: polarization(3)
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

 call trans_rate_1pa(bcorr, cryst, ebands, nw, polarization, scissor, trans_rate, wminzero, wmax, wmesh, comm, dtset, pawtab, psps, wfd, wfk0_path)
 eps(2,:) = trans_rate*(1_dp/wmesh**2)*one/(two_pi)
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

 subroutine d2pa_coefs(alpha2_d, bcorr, cryst, ebands, fname_root, nw, pmat, polarization_1, polarization_2, scissor, wmin, wmax, wmesh, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm, bcorr
 real(dp),intent(in) :: scissor, wmin, wmax
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands

!arrays
 complex(dpc),intent(in) :: pmat(:,:,:,:,:), polarization_1(3), polarization_2(3)
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

 call trans_rate_d2pa(bcorr, cryst, ebands, nw, pmat, polarization_1, polarization_2, scissor, trans_rate, wmin, wmax, wmesh, comm)
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


end module m_nlo


