!!****m* ABINIT/m_nlo
!! NAME
!!  m_nlo
!!
!! FUNCTION
!! DFT calculations of linear and nonlinear optical properties
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2020 ABINIT group (NC)
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
 use m_special_funcs,  only : gaussian, fermi_dirac
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
 public :: trans_rate_2pa
 public :: kk_linperm
 public :: linopt_coefs
 public :: coefs_2pa
 public :: print_matrix_elements
!!***

contains

subroutine trans_rate_1pa(bcorr, cryst, ebands, ngfft, nw, pol, scissor, trans_rate, wmin, wmax, wmesh, comm, &
&                         dtset, pawtab, psps, wfk0_path)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: bcorr,comm,nw
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: scissor,wmin,wmax
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!arrays
 complex(dpc),intent(in) :: pol(3)
 real(dp),intent(out) :: trans_rate(nw)
 real(dp),intent(out) :: wmesh(nw)

!Local variables ------------------------------
!scalars
 integer :: ierr,intmeth 
 integer :: mpw,npw_k,istwf_k
 integer :: my_rank,nproc,mband
 integer :: isppol,sband,fband,iband 
 integer :: isym,itim,timrev,idir,jdir
 integer, parameter :: master=0,intmeth_tetra=2,intmeth_gauss=1
 real(dp) :: step,population_diff,renorm_factor,max_occ
 real(dp) :: fermi_electrons,fermi_holes,carrier_temp
 type(htetra_t) :: htetra
 type(ddkop_t) :: ddkop 
 type(wfd_t) :: wfd
 integer :: usecprj
 
!arrays
 character(len=500) :: msg 
 complex(dpc) :: vk(3,ebands%mband,ebands%mband),vk_sf(3),transrate_tensor(nw,3,3),pol_rot(3),transrate_times_polrot(3),integrand(nw,3,3)
 integer :: ik,my_start,my_stop,iw
 real(dp) :: weights(nw,2),energy_fs(ebands%nkpt)
 real(dp) :: v_bks(2,3),sym_inv_trans(3,3)
 real(dp),allocatable :: cg_ket(:,:),cg_bra(:,:)
 complex(dpc),allocatable :: pol_rotated(:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 integer,allocatable :: wfd_istwfk(:),nband(:,:)
 complex(dpc) :: anan = (zero, zero)

!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 

 ! Trim buffer from maximum band number
 mband = dtset%mband - dtset%nbdbuf

 ! Implement (band,kpt,spin) mask to read relevant wave functions from wfd
 ABI_MALLOC(nband, (ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, ebands%nkpt, ebands%nsppol))
 keep_ur = .False.

 ! Assign each processor its k points so each can read only the wave functions it will use
 call xmpi_split_work(ebands%nkpt,comm,my_start,my_stop)
 ! Read in wave functions
 bks_mask = .false.
 bks_mask(1:mband,my_start:my_stop,:) = .true.

 ABI_MALLOC(wfd_istwfk, (ebands%nkpt))
 wfd_istwfk = 1
 nband = dtset%mband

 ! Initialize wfd and read wave functions
 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, dtset%mband, nband, ebands%nkpt, ebands%nsppol, bks_mask,&
               dtset%nspden, ebands%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(nband)

 mpw = maxval(wfd%npwarr)

 ! Define the optical frequency grid
 step = (wmax - wmin)/(nw - 1)
 wmesh = arth(wmin, step, nw)

 ! Initialize vectors for fourier coefficients (c_G) of wave functions
 ABI_CALLOC(cg_ket, (2, mpw*ebands%nspinor))
 ABI_CALLOC(cg_bra, (2, mpw*ebands%nspinor))
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (cryst%natom, ebands%nspinor*usecprj))

 ! Create ddk object
 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 ! Apply scissor shift to band structure
 call ebands_apply_scissors(ebands, scissor)

 timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 ABI_MALLOC(pol_rotated, (3, cryst%nsym*(timrev + 1)))

 ! Take the polarization vector and rotate by all crystal symmetries
 do itim = 0,timrev
   do isym=1,cryst%nsym
     call matr3inv(cryst%symrel_cart(:,:,isym), sym_inv_trans)
     pol_rotated(:,cryst%nsym*itim + isym) = (1-2*itim)*matmul(transpose(sym_inv_trans), pol)
   end do !isym
 end do !itim

 max_occ = two/ebands%nspinor

 ! TODO: This loop is an unfinished feature to calculate the optical properties at finite temperature
 !       or in the presence of excited carriers. Ideally, it should be possible for electrons and holes
 !       to have a different temperature to study optical response of a thermalizing material
 if (.false.) then
   fermi_holes = dtset%dfield(1)
   fermi_electrons = dtset%dfield(2)
   carrier_temp = dtset%dfield(3)
   do isppol = 1,ebands%nsppol
     do ik = my_start,my_stop
       do sband = 1,mband
         if (sband .le. 4) then
           ebands%occ(sband,ik,isppol) = max_occ*fermi_dirac(ebands%eig(sband,ik,isppol), fermi_holes, carrier_temp)
         else
           ebands%occ(sband,ik,isppol) = max_occ*fermi_dirac(ebands%eig(sband,ik,isppol), fermi_electrons, carrier_temp)
         end if
       end do !sband
     end do !ik
   end do !isppol
 end if

 ! Set integration method (tetrahedral or gaussian)
 intmeth = dtset%eph_intmeth
 ! Create tetrahedra object for integration
 if (intmeth == intmeth_tetra) then
   htetra = tetra_from_kptrlatt(cryst, ebands%kptopt, ebands%kptrlatt, ebands%nshiftk, &
&                               ebands%shiftk, ebands%nkpt, ebands%kptns, comm, msg, ierr)
    if (ierr /= 0) ABI_ERROR(msg)
 end if
 transrate_tensor = zero
 do isppol = 1,ebands%nsppol   
   do ik = my_start,my_stop
     ! vk(3,mband,mband) are matrix elements between bands at given ik (in irreducible wedge) 
     vk = zero
     npw_k = wfd%npwarr(ik); istwf_k = wfd%istwfk(ik)
     ! The ik loop is outside the band look so this ddk setup only needs to be called once (or twice with spin) per point
     call ddkop%setup_spin_kpoint(dtset, cryst, psps, isppol, ebands%kptns(:,ik), istwf_k, npw_k, wfd%kdata(ik)%kg_k)
     ! Precompute all needed matrix elements in the irreducible BZ
     do sband = 1,mband
       ! Get the wave function Fourier coefficients of the initial band and apply the ddk operator
       call wfd%copy_cg(sband,ik,isppol,cg_ket)
       call ddkop%apply(ebands%eig(sband,ik,isppol), npw_k, ebands%nspinor, cg_ket, cwaveprj0)
       do fband = sband,mband 
         ! Don't bother continuing if bands have the same population
         if (abs(ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)) .lt. 1.0d-12) cycle
         ! Read in the final band and compute the bracket <f|ddk|s>
         call wfd%copy_cg(fband,ik,isppol,cg_bra)
         v_bks = ddkop%get_braket(ebands%eig(sband,ik,isppol),istwf_k, npw_k, ebands%nspinor, cg_bra, mode="cart")
         ! Adjust the optical matrix elements to be consistent with the scissor shift
         ! TODO: Determine the appropriate adjustment function
         renorm_factor = one + (scissor)/(ebands%eig(fband,ik,isppol) - ebands%eig(sband,ik,isppol) - scissor)
         ! Convert to complex data type and apply renormalization
         vk(:,sband,fband) = cmplx(v_bks(1,:), v_bks(2,:),kind=dp)*(renorm_factor)
       end do !fband
     end do !sband

     ! Loop over bands on outside allows us to only compute tetrahedron weights once
      do sband = 1,mband
       do fband = sband,mband
         population_diff = ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)
         ! Population filter again
         if (abs(population_diff) .lt. 1.0d-12) cycle
         energy_fs = ebands%eig(fband,:,isppol) - ebands%eig(sband,:,isppol)

         ! Calculate weights for either tetrahedral or gaussian integration
         if (intmeth == intmeth_tetra) then
           call htetra%get_onewk_wvals(ik, 0, nw, wmesh, max_occ, ebands%nkpt, energy_fs, weights)
         else if (intmeth == intmeth_gauss) then
           weights(:,1) = max_occ*gaussian(wmesh - energy_fs(ik), 0.0005_dp)
         end if
         ! Skip if weight is negligible
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
 if (ierr /= 0) ABI_ERROR("Error in xmpi sum for transition rate")

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
  ABI_FREE(cwaveprj0)
  ABI_FREE(pol_rotated)
  ABI_FREE(cg_ket)
  ABI_FREE(cg_bra)

end subroutine trans_rate_1pa

subroutine trans_rate_2pa(bcorr, cryst, ebands, ngfft, nw, pol1, pol2, scissor, trans_rate, wmesh1, wmesh2, comm, &
&                         dtset, pawtab, psps, wfk0_path)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: bcorr,comm,nw
 real(dp),intent(in) :: scissor 
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!arrays
 integer,intent(in) :: ngfft(18)
 complex(dpc),intent(in) :: pol1(3),pol2(3)
 real(dp),intent(out) :: trans_rate(nw)
 real(dp),intent(in) :: wmesh1(nw),wmesh2(nw)

!Local variables ------------------------------
!scalars
 integer :: ierr 
 integer :: mpw,npw_k,istwf_k
 integer :: my_rank,nproc,mband
 integer :: isppol,sband,fband,iband 
 integer :: isym,itim,timrev,idir,jdir,kdir,ldir
 integer :: intmeth
 integer,parameter :: master=0,intmeth_tetra=2,intmeth_gauss=1
 real(dp) :: step,population_diff,renorm_factor,max_occ,energy_is
 type(htetra_t) :: htetra
 type(ddkop_t) :: ddkop 
 type(wfd_t) :: wfd
 integer :: usecprj
 
!arrays
 character(len=500) :: msg 
 complex(dpc) :: vk(3,ebands%mband,ebands%mband),vk_sf(3),transrate_tensor(nw,3,3,3,3),pol1_rot(3),pol2_rot(3),transrate_times_polrot(3),summand(nw,3,3)
 complex(dpc) :: path_sum(nw,3,3),integrand(nw,3,3,3,3),tens_1c(nw,3,3,3),tens_2c(nw,3,3),tens_3c(nw,3),vk_is(3),vk_fi(3)
 complex(dpc) :: detuning_1(nw),detuning_2(nw)
 integer :: ik,my_start,my_stop,iw,nsym
 real(dp) :: klatt(3,3),rlatt(3,3),weights(nw,2),energy_fs(ebands%nkpt)
 real(dp) ::  v_bks(2,3),sym_inv_trans(3,3)
 real(dp),allocatable :: cg_ket(:,:),cg_bra(:,:)
 complex(dpc),allocatable :: pol1_rotated(:,:),pol2_rotated(:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 integer,allocatable :: wfd_istwfk(:),nband(:,:)
!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 

 call xmpi_split_work(ebands%nkpt,comm,my_start,my_stop)

 ABI_MALLOC(nband, (ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, ebands%nkpt, ebands%nsppol))

 bks_mask = .False.; keep_ur = .False.
 bks_mask(1:dtset%mband,my_start:my_stop,:) = .true.
 ABI_MALLOC(wfd_istwfk, (ebands%nkpt))
 wfd_istwfk = 1
 nband = dtset%mband

 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, dtset%mband, nband, ebands%nkpt, ebands%nsppol, bks_mask,&
               dtset%nspden, ebands%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(nband)

 mpw = maxval(wfd%npwarr)

 ABI_CALLOC(cg_ket, (2, mpw*ebands%nspinor))
 ABI_CALLOC(cg_bra, (2, mpw*ebands%nspinor))
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (cryst%natom, ebands%nspinor*usecprj))

 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)
 call ebands_apply_scissors(ebands, scissor)

 mband = ebands%mband - dtset%nbdbuf
 max_occ = two/ebands%nspinor

 intmeth = dtset%eph_intmeth
 transrate_tensor = zero

 timrev = kpts_timrev_from_kptopt(ebands%kptopt)
 ABI_MALLOC(pol1_rotated, (3, cryst%nsym*(timrev + 1)))
 ABI_MALLOC(pol2_rotated, (3, cryst%nsym*(timrev + 1)))

 do itim = 0,timrev
   do isym=1,cryst%nsym
     call matr3inv(cryst%symrel_cart(:,:,isym), sym_inv_trans)
     pol1_rotated(:,cryst%nsym*itim + isym) = (1-2*itim)*matmul(cryst%symrel_cart(:,:,isym), pol1)
     pol2_rotated(:,cryst%nsym*itim + isym) = (1-2*itim)*matmul(cryst%symrel_cart(:,:,isym), pol2)
   end do !isym
 end do !itim

 ! Create tetrahedra object for integration
 if (intmeth == intmeth_tetra) then
   htetra = tetra_from_kptrlatt(cryst, ebands%kptopt, ebands%kptrlatt, ebands%nshiftk, &
&                               ebands%shiftk, ebands%nkpt, ebands%kptns, comm, msg, ierr)
   if (ierr /= 0) ABI_ERROR(msg)
 end if

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
         !if (abs(ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)) .lt. 1.0d-12) cycle
         call wfd%copy_cg(fband,ik,isppol,cg_bra)
         v_bks = ddkop%get_braket(ebands%eig(sband,ik,isppol),istwf_k, npw_k, ebands%nspinor, cg_bra, mode="cart")
         renorm_factor = (ebands%eig(fband,ik,isppol) - ebands%eig(sband,ik,isppol))/(ebands%eig(fband,ik,isppol) - ebands%eig(sband,ik,isppol) - scissor)
         vk(:,sband,fband) = cmplx(v_bks(1,:), v_bks(2,:),kind=dp)*(renorm_factor)
         vk(:,fband,sband) = conjg(vk(:,sband,fband))
       end do !fband
     end do !sband
     
     ! Loop over bands on outside allows us to only compute tetrahedron weights once
     !do sband = 1,mband
     do sband = 1,mband 
       ! fband loop starts as sband to avoid double counting transitions
       do fband = sband,mband
         population_diff = ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)
         ! Population filter again
         if (abs(population_diff) .lt. 1.0d-12) cycle
         energy_fs = ebands%eig(fband,:,isppol) - ebands%eig(sband,:,isppol)
         ! get weights
         weights = zero
         if (intmeth == 1) then
           call htetra%get_onewk_wvals(ik, 0, nw, wmesh1 + wmesh2, max_occ, ebands%nkpt, energy_fs, weights)
         else if (intmeth == 2) then
           weights(:,1) = max_occ*gaussian(wmesh1 + wmesh2 - energy_fs(ik), 0.0005_dp)
         end if
         ! filter based on weights
         if (all(abs(weights(:,1)) .lt. 1.0d-12)) cycle
         path_sum = zero
         do iband=1,mband
           energy_is = ebands%eig(sband,ik,isppol) - ebands%eig(iband,ik,isppol)
           vk_is = vk(:,iband,sband)
           vk_fi = vk(:,fband,iband)
           !TODO: Dynamic imaginary regularization factor in denominator
           do idir = 1,3
             do jdir = 1,3
               summand(:,idir,jdir) = vk_fi(idir)*vk_is(jdir)/(energy_is - wmesh1 + (0,0.005))  &
&                                     + vk_fi(jdir)*vk_is(idir)/(energy_is - wmesh2 + (0,0.005))
             end do !idir
           end do !jdir
           path_sum = path_sum + summand
         end do !iband
         integrand = zero
         do idir=1,3
           do jdir=1,3
             do kdir=1,3
               do ldir=1,3
                 integrand(:,idir,jdir,kdir,ldir) = conjg(path_sum(:,idir,jdir))*path_sum(:,kdir,ldir)*weights(:,1)
               end do ! ldir
             end do !kdir
           end do !jdir
         end do !idir
         transrate_tensor = transrate_tensor + integrand*population_diff*ebands%wtk(ik)
       end do !fband
     end do !sband
   end do !ik
 end do !isppol
 call xmpi_sum(transrate_tensor, comm, ierr)
 if (ierr /= 0) ABI_ERROR("Error in xmpi sum for transition rate")
 
 trans_rate = zero

 nsym = cryst%nsym

 if (ebands%kptopt == 3) then
   nsym = 1
   timrev = 0
 end if
 ! Contract the tensor for the total two-photon transition rate
 do itim=0,timrev
   do isym=1,nsym
     pol1_rot = pol1_rotated(:,nsym*itim + isym)
     pol2_rot = pol2_rotated(:,nsym*itim + isym)

     tens_1c = zero
     do idir=1,3
       tens_1c = tens_1c + transrate_tensor(:,:,:,:,idir)*pol1_rot(idir)
     end do
     tens_2c = zero
     do idir=1,3
       tens_2c = tens_2c + tens_1c(:,:,:,idir)*pol1_rot(idir)
     end do
     tens_3c = zero
     do idir=1,3
       tens_3c = tens_3c + tens_2c(:,:,idir)*pol2_rot(idir)
     end do
     do idir=1,3
       trans_rate = trans_rate + tens_3c(:,idir)*pol2_rot(idir)
     end do
   end do !isym
 end do !itim
 trans_rate = trans_rate/(nsym*(timrev + 1))

 call htetra%free
 call pawcprj_free(cwaveprj0)
 ABI_FREE(pol1_rotated)
 ABI_FREE(pol2_rotated)
 ABI_FREE(cg_ket)
 ABI_FREE(cg_bra)

end subroutine trans_rate_2pa

subroutine kk_linperm(eps_imag, len_abs, eps_real, wmesh)

 integer,intent(in) :: len_abs
 real(dp),intent(in) :: eps_imag(len_abs)
 real(dp),intent(out) :: eps_real(len_abs)
 real(dp),intent(in) :: wmesh(len_abs)

 integer :: iw, iwp
 real(dp) :: wstep, w, wp
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
 integer,intent(in) :: bcorr,comm,nw
 real(dp),intent(in) :: scissor,wmax
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
 integer :: my_rank,nproc 
 integer :: iw,funit
 integer,parameter :: master=0
 real(dp),parameter :: wminzero=zero,c_au = 137_dp
 
!arrays
 character(len=500) :: msg,errmsg,fname
 real(dp) :: trans_rate(nw),eps_mag(nw)
 real(dp) :: index_alpha(2,nw),eps(2,nw),wmesh(nw)

!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 

 call trans_rate_1pa(bcorr, cryst, ebands, ngfft, nw, pol, scissor, trans_rate, wminzero, wmax, wmesh, comm, dtset, pawtab, psps, wfk0_path)
 eps(2,:) = trans_rate*(1_dp/wmesh**2)*one/(two_pi)
 eps(2,1) = 0
 call kk_linperm(eps(2,:), nw, eps(1,:), wmesh)
 eps_mag = sqrt(eps(1,:)**2 + eps(2,:)**2)

 index_alpha(1,:) = sqrt(one/two*(eps_mag + eps(1,:)))
 index_alpha(2,:) = sqrt(one/two*(eps_mag - eps(1,:)))*two*wmesh/c_au

 if (alpha_units_per_cm==1) index_alpha(2,:) = index_alpha(2,:)/5.29d-9

 if (my_rank == master) then
   funit = get_unit()
   fname = 'linopt_OUT' 
   open(funit, file = fname)
   do iw = 1,nw
     write(funit, *) wmesh(iw), index_alpha(2,iw), index_alpha(1,iw), eps(1,iw), eps(2,iw), trans_rate(iw)
   end do
   close(funit)
 end if

end subroutine linopt_coefs

subroutine coefs_2pa(bcorr, cryst, ebands, fname_root, ngfft, nw, pawtab, pol1, pol2, scissor, w1min, w1max, comm, dtset, psps, wfk0_path, w2)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm, bcorr
 real(dp),intent(in) :: scissor, w1min, w1max
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 character(len=*),intent(in) :: wfk0_path
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 real(dp),intent(in),optional :: w2

!arrays
 complex(dpc),intent(in) :: pol1(3), pol2(3)
 integer,intent(in) :: nw
 character(len=*),intent(in) :: fname_root
 integer,intent(in) :: ngfft(18)

!Local variables ------------------------------
!scalars
 integer :: ierr, funit, iw
 integer :: my_rank, nproc, pump_index 
 integer, parameter :: master=0
 real(dp) :: step, wmesh1(nw), wmesh2(nw)
 real(dp),parameter :: c_au=137_dp,au_to_cmGW=28.28
 type(wfd_t) :: wfd
 
!arrays
 character(len=500) :: msg, errmsg, fname
 real(dp) :: trans_rate(nw), alpha2(nw)

!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 
 
 step = (w1max - w1min)/(nw - 1)
 wmesh1 = arth(w1min, step, nw)

 ! If second frequency argument not present, degenerate 2PA is calculated
 if (present(w2)) then
   wmesh2 = w2
 else
   wmesh2 = wmesh1
 end if

   call trans_rate_2pa(bcorr, cryst, ebands, ngfft, nw, pol1, pol2, scissor, trans_rate, wmesh1, wmesh2, comm, &
&                      dtset, pawtab, psps, wfk0_path)

 !TODO: Refractive index is currently hard coded.
 !      Need two options: compute the index using linopt_coefs, use constant value,
 !                        or allow user to input dispersion
 alpha2 = 4_dp/(3.2**2*c_au**2*wmesh1**3)*trans_rate
 alpha2(1) = zero
 alpha2 = alpha2*au_to_cmGW

 !TODO: Output NetCDF file
 if (my_rank == master) then
   fname = 'nd2pa_OUT' 
   funit = get_unit()
   open(funit, file = fname)
   do iw = 1,nw
     write(funit, *) wmesh1(iw),  wmesh2(iw), alpha2(iw), trans_rate(iw)
   end do
   close(funit)
 end if

end subroutine coefs_2pa

!TODO: Put reused initialization steps (reading wfd, initializing ddk) here
subroutine init_optics()

end subroutine init_optics

subroutine print_matrix_elements(cryst, ebands, ngfft, comm, dtset, pawtab, psps, wfk0_path)
                         

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
 integer,intent(in) :: ngfft(18)
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(inout) :: ebands
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)


!Local variables ------------------------------
!scalars
 integer :: ierr 
 integer :: mpw,npw_k,istwf_k
 integer :: my_rank,nproc,mband
 integer :: isppol,sband,fband,iband 
 integer :: isym,itim,timrev,idir,jdir
 integer,parameter :: master=0
 real(dp) :: step,population_diff,renorm_factor,max_occ
 type(htetra_t) :: htetra
 type(ddkop_t) :: ddkop 
 type(wfd_t) :: wfd
 integer :: usecprj,intmeth
 
!arrays
 character(len=500) :: msg,fname 
 complex(dpc),allocatable :: vk(:,:,:),v(:,:,:,:,:)
 integer :: ik,my_start,my_stop,funit
 real(dp) ::  v_bks(2,3),sym_inv_trans(3,3)
 real(dp),allocatable :: cg_ket(:,:),cg_bra(:,:)
 complex(dpc),allocatable :: pol_rotated(:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 integer,allocatable :: wfd_istwfk(:),nband(:,:)
 complex(dpc) :: anan = (zero,zero)
!************************************************************************
 my_rank = xmpi_comm_rank(comm)
 nproc   = xmpi_comm_size(comm) 

 mband = dtset%mband - dtset%nbdbuf
 ABI_MALLOC(nband, (ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, ebands%nkpt, ebands%nsppol))
 keep_ur = .False.

 call xmpi_split_work(ebands%nkpt,comm,my_start,my_stop)

 bks_mask = .false.
 bks_mask(1:mband,my_start:my_stop,:) = .true.

 ABI_MALLOC(wfd_istwfk, (ebands%nkpt))
 wfd_istwfk = 1
 nband = dtset%mband

 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, dtset%mband, nband, ebands%nkpt, ebands%nsppol, bks_mask,&
               dtset%nspden, ebands%nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)
 ABI_FREE(nband)

 mpw = maxval(wfd%npwarr)

 ABI_CALLOC(cg_ket, (2, mpw*ebands%nspinor))
 ABI_CALLOC(cg_bra, (2, mpw*ebands%nspinor))
 usecprj = 0
 ABI_MALLOC(cwaveprj0, (cryst%natom, ebands%nspinor*usecprj))

 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 ABI_MALLOC(vk, (mband, mband, 3))
 ABI_MALLOC(v, (mband, mband, ebands%nkpt, 3, ebands%nsppol))

 do isppol = 1,ebands%nsppol   
   !For debugging purposes, set the vk array to NaN to make sure we're actually assigning
   vk = anan/anan
   do ik = 1,ebands%nkpt
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
         !if (abs(ebands%occ(sband,ik,isppol) - ebands%occ(fband,ik,isppol)) .lt. 1.0d-12) cycle
         call wfd%copy_cg(fband,ik,isppol,cg_bra)
         v_bks = ddkop%get_braket(ebands%eig(sband,ik,isppol),istwf_k, npw_k, ebands%nspinor, cg_bra, mode="cart")
         print *, sband, fband, ik
         print *, v_bks
         ! renorm? Just seems to make errors larger.
         !renorm_factor = (ebands%eig(fband,ik,isppol) - ebands%eig(sband,ik,isppol))/(ebands%eig(fband,ik,isppol) - ebands%eig(sband,ik,isppol) - scissor)
         renorm_factor = one
         vk(fband,sband,:) = cmplx(v_bks(1,:), v_bks(2,:),kind=dp)*(renorm_factor)
         vk(sband,fband,:) = conjg(vk(fband,sband,:))
       end do !fband
     end do !sband
   v(:,:,ik,:,isppol) = vk
   end do !ik
 end do !isppol
 
 call pawcprj_free(cwaveprj0)
 ABI_FREE(cg_ket)
 ABI_FREE(cg_bra)

 if (my_rank == master) then
   fname = "matrix_elements.out"
   funit = get_unit()
   open(funit, file = fname)
   write(funit, *) v
   close(funit)
 end if

  ABI_FREE(v)
  ABI_FREE(vk)
end subroutine print_matrix_elements

end module m_nlo


