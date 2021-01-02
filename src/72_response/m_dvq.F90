!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_dvq

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_dtset
 use m_crystal
 use m_hamiltonian
 use m_ebands
 use m_getgh1c
 use m_dvdb
 use m_mpinfo
 use m_initylmg
 use m_pawcprj
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use defs_abitypes,   only : MPI_type
 use defs_datatypes,  only : ebands_t, pseudopotential_type
 use m_pawtab,        only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
 use m_fstrings,      only : ktoa, sjoin
 implicit none

 private

 type, private :: ham_targets_t
   real(dp),allocatable :: ffnlk(:,:,:,:), ffnl1(:,:,:,:)
   real(dp),allocatable :: kpg_k(:,:), kpg1_k(:,:)
   real(dp),allocatable :: ph3d(:,:,:), ph3d1(:,:,:)
   real(dp),allocatable :: dkinpw(:), kinpw1(:)
   contains
     procedure :: free => ham_targets_free   ! Free memory.
 end type ham_targets_t

!!****t* m_dvq/dvqop_t
!! NAME
!!  dvqop_t
!!
!! FUNCTION
!!  This object provides a simplified interface to compute matrix elements of the
!!  velocity operator with the DFPT routines.
!!
!! SOURCE

 type,public :: dvqop_t

  integer :: inclvkb
  ! Option for calculating the matrix elements of [Vnl,r].
  ! 0 to exclude commutator, 2 to include it

  integer :: usepaw
  ! 0 for NC, 1 for PAW

  integer :: mpw
  ! Maximum number of plane-waves over k-points (used to dimension arrays)

  integer :: cplex

  integer :: nspden
 
  integer :: natom3

  integer :: ngfft(18)
 
  integer :: nfft

  integer :: nfftf

  integer :: ngfftf(18)

  integer :: ddb_ngqpt(3)

  real(dp) :: ddb_shiftq(3)

  integer :: use_ftinterp

  type(dvdb_t),pointer :: dvdb

  real(dp) :: kpoint(3)
  ! K-point (set in setup_spin_kpoint)

  real(dp) :: eig0nk

  real(dp) :: dfpt_sciss = zero

  real(dp) :: rprimd(3,3)

  type(MPI_type),pointer :: mpi_enreg => null()

  type(gs_hamiltonian_type),allocatable :: gs_hamkq(:)

  type(rf_hamiltonian_type),allocatable :: rf_hamkq(:)

  type(ham_targets_t), private, allocatable :: htg(:)
  ! Store arrays targetted by the hamiltonians.

  real(dp), private, allocatable :: gh1c(:,:,:)
   !gh1c, (2, mpw*nspinor, 3*natom))

  real(dp), private, allocatable :: gs1c(:,:,:)
   ! gs1c, (2, mpw*nspinor, 3*psps%usepaw))

 contains

  procedure :: load_vlocal1 => dvqop_load_vlocal1 
  !Load in the potential for each atomic perturbation

  procedure :: setup_spin_kpoint => dvqop_setup_spin_kpoint

  procedure :: apply => dvqop_apply


 end type dvqop_t

 public :: dvqop_new

 CONTAINS

!----------------------------------------------------------------------

!!****f* m_dvq/dvqop_new
!! NAME
!!  ddkop_new
!!
!! FUNCTION
!!  Build new object. Use dtset%inclvkb to determine whether non-local part should be included.
!!
!! INPUTS
!! dtset<dataset_type>=All input variables for this dataset.
!! cryst<crystal_t>=Crystal structure.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! mpi_enreg=information about MPI parallelization
!! mpw=Maximum number of plane-waves over k-points.
!! ngfft(18)=contain all needed information about 3D FFT
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(dvqop_t) function dvqop_new(dtset, dvdb, cryst, pawtab, psps, mpi_enreg, mpw, ngfft, ngfftf) result(new)

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),target,intent(in) :: mpi_enreg
 integer,intent(in) :: mpw
 type(dvdb_t),intent(inout),target :: dvdb
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: cplex1 = 1
 integer :: nfft, mgfft, ipert, natom3

! *************************************************************************

 ABI_CHECK(dtset%usepaw == 0, "PAW not tested/implemented!")

 new%inclvkb = dtset%inclvkb
 new%usepaw = dtset%usepaw
 new%dfpt_sciss = dtset%dfpt_sciss
 new%mpw = mpw

 new%rprimd = cryst%rprimd
 new%mpi_enreg => mpi_enreg
 new%nspden = dtset%nspden

 natom3 = 3*cryst%natom
 new%natom3 = natom3

 nfft = product(ngfft(1:3))
 mgfft = maxval(ngfft(1:3))

 new%nfft = nfft
 new%nfftf = product(ngfftf(1:3))

 new%ngfft = ngfft
 new%ngfftf = ngfftf
 new%cplex = cplex1
 new%use_ftinterp = dtset%eph_use_ftinterp
 new%ddb_ngqpt = dtset%ddb_ngqpt
 new%ddb_shiftq = dtset%ddb_shiftq

 new%dvdb => dvdb
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 ABI_MALLOC(new%gh1c, (2, new%mpw*dtset%nspinor, natom3))
 ABI_MALLOC(new%gs1c, (2, new%mpw*dtset%nspinor, natom3))

 ABI_MALLOC(new%gs_hamkq, (natom3))
 ABI_MALLOC(new%rf_hamkq, (natom3))

 ABI_MALLOC(new%htg, (natom3))
 do ipert=1,natom3
   ! ==== Initialize most of the Hamiltonian (and derivative) ====
   ! 1) Allocate all arrays and initialize quantities that do not depend on k and spin.
   ! 2) Perform the setup needed for the non-local factors:
   ! * Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
   ! * PAW: Initialize the overlap coefficients and allocate the Dij coefficients.
   call init_hamiltonian(new%gs_hamkq(ipert), psps, pawtab, dtset%nspinor, dtset%nsppol, dtset%nspden, cryst%natom,&
     cryst%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg)
     !paw_ij=paw_ij,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
     !usecprj=usecprj,ph1d=ph1d,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

   ! Prepare application of the NL part.
   call init_rf_hamiltonian(cplex1, new%gs_hamkq(ipert), ipert, new%rf_hamkq(ipert), has_e1kbsc=.true.)
     !&paw_ij1=paw_ij1,comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
     !&mpi_spintab=mpi_enreg%my_isppoltab)
 end do


end function dvqop_new

subroutine dvqop_load_vlocal1(self,qpt,spin,pawfgr,comm)

 class(dvqop_t),intent(inout) :: self
 integer,intent(in) :: spin
 integer :: nfftf
 real(dp),intent(in) :: qpt(3)
 type(pawfgr_type),intent(in) :: pawfgr
 integer,intent(in) :: comm

 integer :: ipc,cplex,interpolated,comm_rpt,db_iqpt
 real(dp) :: vdummy(self%nfftf,self%nspden), vdummy_reshaped(self%ngfft(4),self%ngfft(5),self%ngfft(6),self%gs_hamkq(1)%nvloc)
 real(dp) :: vlocal1_reshaped(self%cplex*self%ngfft(4),self%ngfft(5),self%ngfft(6),self%gs_hamkq(1)%nvloc)
 real(dp),allocatable :: vlocal1(:,:,:,:)
 
 vdummy = zero
 
 if (self%use_ftinterp /= 0) then
   MSG_WARNING(sjoin("Enforcing FT interpolation for q-point", ktoa(qpt)))
   comm_rpt = xmpi_comm_self
   call self%dvdb%ftinterp_setup(self%ddb_ngqpt, 1, self%ddb_shiftq, self%nfftf, self%ngfftf, comm_rpt)
   cplex = 2
   ABI_MALLOC(vlocal1, (cplex, self%nfftf, self%nspden, self%dvdb%my_npert))
   call self%dvdb%ftinterp_qpt(qpt, self%nfftf, self%ngfftf, vlocal1, self%dvdb%comm_rpt)
   interpolated = 1
 else
   ! Find the index of the q-point in the DVDB.
   db_iqpt = self%dvdb%findq(qpt)
   if (db_iqpt /= -1) then
     ! Read or reconstruct the dvscf potentials for all 3*natom perturbations.
     ! This call allocates vlocal1(cplex, nfftf, nspden, 3*natom))
     call self%dvdb%readsym_allv1(db_iqpt, self%cplex, self%nfftf, self%ngfftf, vlocal1, comm)
   else
     MSG_WARNING(sjoin("Cannot find q-point:", ktoa(qpt), "in DVDB file"))
   end if
 end if

 do ipc=1,self%natom3
   call rf_transgrid_and_pack(spin,self%nspden,self%usepaw,self%cplex,self%nfftf,self%nfft, self%ngfft,&
                              self%gs_hamkq(ipc)%nvloc, pawfgr, self%mpi_enreg, vdummy, vlocal1, vdummy_reshaped, vlocal1_reshaped)
   call self%rf_hamkq(ipc)%load_spin(spin,vlocal1=vlocal1_reshaped,with_nonlocal=.true.)
   call self%gs_hamkq(ipc)%load_spin(spin,vlocal=vdummy_reshaped,with_nonlocal=.true.)
 end do
 ABI_FREE(vlocal1)
 

end subroutine dvqop_load_vlocal1
!!***
!----------------------------------------------------------------------

!!****f* m_dvq/dvqop_setup_spin_kpoint
!! NAME
!!  dvqop_setup_spin_kpoint
!!
!! FUNCTION
!!  Prepare internal tables that depend on k-point/spin
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  cryst<crystal_t>=Crystal structure.
!!  psps<pseudopotential_type>=Variables related to pseudopotentials.
!!  spin: spin index
!!  kpoint(3): K-point in reduced coordinates.
!!  istwkf_k: defines storage of wavefunctions for this k-point
!!  npw_k: Number of planewaves.
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvqop_setup_spin_kpoint(self, dtset, cryst, psps, spin, kpoint, kqpoint, istwf_k, istwf_kq, npw_k, npw_kq, kg_k, kg_kq)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: spin, npw_k, npw_kq, istwf_k, istwf_kq
 class(dvqop_t),intent(inout) :: self
 type(crystal_t) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg_k(3,npw_k),kg_kq(3,npw_k)
 real(dp),intent(in) :: kpoint(3), kqpoint(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: nkpt1=1, nsppol1=1
 type(mpi_type) :: mpienreg_seq
!arrays
 integer :: npwarr(nkpt1), dummy_nband(nkpt1*nsppol1)
 integer :: nkpg, nkpg1, useylmgr1, optder !, nylmgr1
 integer :: ipc, ipert, idir
 real(dp),allocatable :: ylm_k(:,:),ylm_kq(:,:),ylmgr1_k(:,:,:),ylmgr1_kq(:,:,:)


!************************************************************************

 ABI_CHECK(npw_k <= self%mpw, "npw_k > mpw!")
 ABI_CHECK(npw_kq <= self%mpw, "npw_kq > mpw!")
 self%kpoint = kpoint

 ! Set up the spherical harmonics (Ylm) at k+q if useylm = 1
 useylmgr1 = 0; optder = 0
 if (psps%useylm == 1) then
   useylmgr1 = 1; optder = 1
 end if

 ABI_MALLOC(ylm_k, (npw_k, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylm_kq, (npw_kq, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylmgr1_k, (npw_k, 3+6*(optder/2), psps%mpsang**2*psps%useylm*useylmgr1))
 ABI_MALLOC(ylmgr1_kq, (npw_kq, 3+6*(optder/2), psps%mpsang**2*psps%useylm*useylmgr1))

 if (psps%useylm == 1) then
   ! Fake MPI_type for sequential part. dummy_nband and nsppol1 are not used in sequential mode.
   call initmpi_seq(mpienreg_seq)
   dummy_nband = 0; npwarr = npw_k
   call initylmg(cryst%gprimd, kg_k, kpoint, nkpt1, mpienreg_seq, psps%mpsang, npw_k, dummy_nband, nkpt1, &
      npwarr, nsppol1, optder, cryst%rprimd, ylm_k, ylmgr1_k)
   call destroy_mpi_enreg(mpienreg_seq)
 end if

 do ipc=1,self%natom3
   idir = mod(ipc-1, 3) + 1
   ipert = (ipc - idir) / 3 + 1
   call self%htg(ipc)%free()

   ! Continue to initialize the Hamiltonian
   call self%gs_hamkq(ipc)%load_spin(spin, with_nonlocal=.true.)
   call self%rf_hamkq(ipc)%load_spin(spin, with_nonlocal=.true.)

   ! We need ffnl1 and dkinpw for all dir,atom. Note that the Hamiltonian objects use pointers to keep a reference
   ! to the output results of this routine.
   ! This is the reason why we need to store the targets in self%htg
 call getgh1c_setup(self%gs_hamkq(ipc), self%rf_hamkq(ipc), dtset, psps, kpoint, kqpoint, idir, ipert, & ! In
   cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, npw_k, npw_kq, &            ! In
   useylmgr1, kg_k, ylm_k, kg_kq, ylm_kq, ylmgr1_k, &                                       ! In
   self%htg(ipc)%dkinpw, nkpg, nkpg1, self%htg(ipc)%kpg_k, self%htg(ipc)%kpg1_k, &     ! Out
   self%htg(ipc)%kinpw1, self%htg(ipc)%ffnlk, self%htg(ipc)%ffnl1, &                   ! Out
   self%htg(ipc)%ph3d, self%htg(ipc)%ph3d1)                                             ! Out
 end do

 ABI_FREE(ylm_k)
 ABI_FREE(ylmgr1_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr1_kq)

end subroutine dvqop_setup_spin_kpoint
!!***
!----------------------------------------------------------------------

!!****f* m_dvq/dvqop_apply
!! NAME
!!  dvqop_apply
!!
!! FUNCTION
!!
!! INPUTS
!!  eig0nk: Eigenvalue associated to the wavefunction.
!!  npw_k: Number of planewaves.
!!  nspinor: Number of spinor components.
!!  cwave(2,npw_k*nspinor)=input wavefunction in reciprocal space
!!  cwaveprj(natom,nspinor*usecprj)=<p_lmn|C> coefficients for wavefunction |C> (and 1st derivatives)
!!     if not allocated or size=0, they are locally computed (and not sorted)!!
!!
!! SIDE EFFECTS
!! Stores:
!!  gh1c(2,npw1*nspinor)= <G|H^(1)|C> or  <G|H^(1)-lambda.S^(1)|C> on the k+q sphere
!!                        (only kinetic+non-local parts if optlocal=0)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dvqop_apply(self, eig0nk, npw_k, nspinor, cwave, cwaveprj)

!Arguments ------------------------------------
!scalars
 class(dvqop_t),intent(inout) :: self
 integer,intent(in) :: npw_k, nspinor
 real(dp),intent(in) :: eig0nk
!arrays
 real(dp),intent(inout) :: cwave(2,npw_k*nspinor)
 type(pawcprj_type),intent(inout) :: cwaveprj(:,:)

!Local variables-------------------------------
!scalars
 integer,parameter :: berryopt0 = 0, optlocal0 = 0, tim_getgh1c = 1, usevnl0 = 0, opt_gvnlx1 = 0
 integer :: sij_opt, ispinor, ipws, ipw, optnl
 integer :: ipc, idir, ipert
 real(dp) :: eshift
!arrays
 real(dp) :: grad_berry(2,(berryopt0/4)), gvnlx1(2,usevnl0)
 real(dp),pointer :: dkinpw(:),kinpw1(:)

!************************************************************************

 self%eig0nk = eig0nk

 if (self%inclvkb /= 0) then
   ! optlocal0 = 0: local part of H^(1) is not computed in gh1c=<G|H^(1)|C>
   ! optnl = 2: non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
   ! opt_gvnlx1 = option controlling the use of gvnlx1 array:
   optnl = 2 !; if (self%inclvkb == 0) optnl = 0

   eshift = self%eig0nk - self%dfpt_sciss
   do ipc=1,self%natom3
     idir = mod(ipc-1, 3) + 1
     ipert = (ipc - idir) / 3 + 1
     sij_opt = self%gs_hamkq(ipc)%usepaw
     call getgh1c(berryopt0, cwave, cwaveprj, self%gh1c(:,:,ipc), &
       grad_berry, self%gs1c(:,:,ipc), self%gs_hamkq(ipc), gvnlx1, idir, ipert, eshift, self%mpi_enreg, optlocal0, &
       optnl, opt_gvnlx1, self%rf_hamkq(ipc), sij_opt, tim_getgh1c, usevnl0)
   end do

 else
   MSG_ERROR("inclvkb = 0 not implemented")
 end if

end subroutine dvqop_apply
!!***
subroutine ham_targets_free(self)

!Arguments ------------------------------------
!scalars
 class(ham_targets_t),intent(inout) :: self

!************************************************************************

 ABI_SFREE(self%ffnlk)
 ABI_SFREE(self%ffnl1)
 ABI_SFREE(self%kpg_k)
 ABI_SFREE(self%kpg1_k)
 ABI_SFREE(self%dkinpw)
 ABI_SFREE(self%kinpw1)
 ABI_SFREE(self%ph3d)
 ABI_SFREE(self%ph3d1)

end subroutine ham_targets_free
!!***
 end module m_dvq
