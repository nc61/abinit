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
&                         dtset, pawtab, psps, wfd)

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

!arrays
 complex(dpc),intent(in) :: pol(3)
 real(dp),intent(out) :: trans_rate(nw)
 real(dp),intent(out) :: wmesh(nw)

!Local variables ------------------------------
!scalars
 integer :: ierr 
 integer :: mpw, npw_k, istwf_k
 integer :: my_rank, nproc, mband
 integer :: isppol, sband, fband, iband 
 integer :: isym, itim, timrev, idir, jdir
 integer, parameter :: master=0
 real(dp) :: step, population_diff, renorm_factor, max_occ 
 type(htetra_t) :: htetra
 type(ddkop_t) :: ddkop 
 integer :: usecprj
 
!arrays
 character(len=500) :: msg 
 complex(dpc) :: vk(3,ebands%mband,ebands%mband), vk_sf(3), transrate_tensor(nw,3,3), pol_rot(3), transrate_times_polrot(3),integrand(nw,3,3)
 integer :: ik, my_start, my_stop, iw
 real(dp) :: klatt(3,3), rlatt(3,3), weights(nw,2), energy_fs(ebands%nkpt)
 real(dp) ::  v_bks(2,3)
 real(dp),allocatable :: cg_ket(:,:), cg_bra(:,:)
 complex(dpc),allocatable :: pol_rotated(:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:)
!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 

 mpw = maxval(wfd%npwarr)
 print *, "mpw"
 print *, mpw

 step = (wmax - wmin)/(nw - 1)
 wmesh = arth(wmin, step, nw)

 rlatt = ebands%kptrlatt
 call matr3inv(rlatt, klatt)

 ABI_CALLOC(cg_ket, (2, mpw*ebands%nspinor))
 ABI_CALLOC(cg_bra, (2, mpw*ebands%nspinor))
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (cryst%natom, ebands%nspinor*usecprj))

 call xmpi_split_work(ebands%nkpt,comm,my_start,my_stop)

!!WFD
 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)
 call ebands_apply_scissors(ebands, scissor)

 timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 print *, "Timrev"
 print *, timrev
 ABI_MALLOC(pol_rotated, (3, cryst%nsym*(timrev + 1)))

 do itim = 0,timrev
   do isym=1,cryst%nsym
     pol_rotated(:,cryst%nsym*itim + isym) = (1-2*itim)*matmul(transpose(cryst%symrel_cart(:,:,isym)), pol)
   end do !isym
 end do !itim
 print *, "pol rotated"
 print *, pol_rotated

 htetra = tetra_from_kptrlatt(cryst, ebands%kptopt, ebands%kptrlatt, ebands%nshiftk, ebands%shiftk, ebands%nkpt, ebands%kptns, comm, msg, ierr)
 if (ierr /= 0) MSG_ERROR(msg)
 mband = ebands%mband
 max_occ = two/ebands%nspinor

 transrate_tensor = zero
 if (ebands%kptopt == 3) then !all k points are inequivalent
   ! do nothing
 else if (ebands%kptopt == 1) then
    do isppol = 1,ebands%nsppol   
      do ik = my_start,my_stop
        !vk(3,mband,mband) are matrix elements between bands at given ik (in irreducible wedge) 
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
            vk(:,sband,fband) = cmplx(v_bks(1,:), v_bks(2,:),kind=dp)*(renorm_factor)
          end do !fband
        end do !sband
        
        ! Loop over bands on outside allows us to only compute tetrahedron weights once
        do sband = 1,mband
          ! fband loop starts as sband to avoid double counting transitions
          do fband = sband,mband
            population_diff = ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)
            ! Population filter again
            if (abs(population_diff) .lt. 1.0d-12) cycle
            energy_fs = ebands%eig(fband,:,isppol) - ebands%eig(sband,:,isppol)
            ! get weights
            call htetra%get_onewk_wvals(ik, 0, nw, wmesh, max_occ, ebands%nkpt, energy_fs, weights)
            ! filter based on weights
            if (all(abs(weights(:,1)) .lt. 1.0d-12)) cycle
            vk_sf = vk(:,sband,fband)
            do idir = 1,3
              do jdir = idir,3
                integrand(:,idir,jdir) = weights(:,1)*conjg(vk_sf(idir))*vk_sf(jdir)
              end do !idir
            end do !jdir

            transrate_tensor = transrate_tensor + integrand*population_diff*ebands%wtk(ik)
          end do !fband
        end do !sband
      end do !ik
    end do !isppol
    call xmpi_sum(transrate_tensor, comm, ierr)
    if (ierr /= 0) MSG_ERROR("Error in xmpi sum for transition rate")
 end if 
 do idir=1,3
   do jdir=1,idir-1
     transrate_tensor(:,idir,jdir) = conjg(transrate_tensor(:,jdir,idir))
   end do !jdir
 end do !idir
 
 trans_rate = zero 
 do itim=0,timrev
   do isym=1,cryst%nsym
     pol_rot = pol_rotated(:,cryst%nsym*itim + isym)
     do iw=1,nw
       ! zhemv: double precision complex matrix*vector
       call zhemv('u',3,one,transrate_tensor(iw,:,:),3,pol_rot,1,zero,transrate_times_polrot,1)
       trans_rate(iw) = dot_product(pol_rot, transrate_times_polrot)
     end do ! iw
   end do !isym
 end do !itim

  call htetra%free
  call pawcprj_free(cwaveprj0)
  ABI_FREE(pol_rotated)
  ABI_FREE(cg_ket)
  ABI_FREE(cg_bra)

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

 call trans_rate_1pa(bcorr, cryst, ebands, nw, pol, scissor, trans_rate, wminzero, wmax, wmesh, comm, dtset, pawtab, psps, wfd)
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


