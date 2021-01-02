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
 end do
 ABI_FREE(vlocal1)
 

end subroutine dvqop_load_vlocal1
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
