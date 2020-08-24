! copyright (c) 2013,  los alamos national security, llc (lans)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_benthos
!
!> \brief MPAS ocean benthos
!> \author Nicole Jeffery
!> \date   07/15/2020
!> \details
!> This module contains routines for computing the benthic tracer tendencies
!> and the bottom boundary flux on ecosys tracers due to the benthos
!
!-----------------------------------------------------------------------

module ocn_tracer_benthos

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timer
   use mpas_timekeeping
   use mpas_forcing
   use mpas_stream_manager
   use mpas_log
   use ocn_constants
!   use ocn_framework_forcing

   use bgc_mod
   use bgc_parms
   use benthos_mod
   use benthos_parms

   implicit none

   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: benthos_tracer_compute, &
!             ocn_tracer_ecosys_benthos_flux_compute,  &
             ocn_get_benthos_data,  &
             ocn_benthos_forcing_write_restart,  &
             ocn_benthos_tracer_init

   integer, public:: &
      numBenthicColumnsMax

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!-----------------------------------------------------------------------
!  name the necessary benthos derived types
!  all of these are defined in benthos_mod
!-----------------------------------------------------------------------

  type(benthos_indices_type)    , public :: benthos_indices
  type(benthos_element_indices_type) , public :: benthos_element_indices
  type(benthos_input_type)      , public :: benthos_input
  type(benthos_forcing_type)    , public :: benthos_forcing
  type(benthos_output_type)     , public :: benthos_output
  type(benthos_diagnostics_type), public :: benthos_diagnostic_fields
  type(benthos_flux_diagnostics_type), public :: benthos_flux_diagnostic_fields

! hold indices in tracer pool corresponding to each eco tracer array
  type(benthos_indices_type), public :: benthosIndices
  type(benthos_element_indices_type), public :: elementIndices
  type(BGC_indices_type) :: ecosysIndices


!***********************************************************************

contains

!***********************************************************************
!
!  routine benthos_tracer_compute
!
!> \brief   computes a seafloor tracer tendency due to benthos biogeochemistry
!> \author  Nicole Jeffery
!> \date    07/29/2020
!> \details
!>  This routine computes a seafloor tracer tendency due to the benthos biogeochemistry
!
!-----------------------------------------------------------------------

   subroutine benthos_tracer_compute(benthosTracers, benthosPhPrev3d,benthosTracersStorage, &
      activeTracers, nBenthicTracers, ecosysTracers, &
      nEcosysTracers, forcingPool, nCellsSolve, nBenthicVertLevels, bottomDepthObserved, maxLevelCell, &
      nVertLevels, layerThickness, indexTemperature, indexSalinity, benthosTracersTend,dt,err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      ! 1 dimensional arrays
      real (kind=RKIND), dimension(:), intent(in) :: &
         bottomDepthObserved

      ! two dimensional arrays
      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness   ! zMid?

      ! three dimensional arrays
      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         benthosTracers, & ! benthos BGC
         ecosysTracers,  & ! ocean BGC
         activeTracers     ! T and S

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         benthosPhPrev3d, &    ! pH
         benthosTracersStorage ! Deep storage

      integer, dimension(:), intent(in) :: &
         maxLevelCell      ! last active level

      ! scalars
      integer, intent(in) :: &
         nBenthicTracers,    &
         nEcosysTracers,     &
         nCellsSolve,        &
         nBenthicVertLevels, &
         nVertLevels,        &
         indexTemperature,   &
         indexSalinity

      real (KIND=RKIND), intent(in) :: dt

      type (mpas_pool_type), intent(inout) :: forcingPool

      !
      ! two dimensional pointers
      !
      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
        benthosTracersTend

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), pointer :: benthosAuxiliaryPool  ! additional forcing fields
      type (mpas_pool_type), pointer :: seafloorStoragePool   ! additional forcing fields
      type (mpas_pool_type), pointer :: benthosDiagFieldsLevel1Pool !,  &  ! diagnostics
!                                        benthosDiagFieldsLevel2,  &
!                                        benthosDiagFieldsLevel3,  &
!                                        benthosDiagFieldsLevel4,  &
!                                        benthosDiagFieldsLevel5

      logical, pointer ::  &
         config_set_benthos_flags, &
         config_use_carbonate_saturation, &
         config_use_secondary_reactions, &
         config_use_depth_dependent_porosity, &
         config_use_non_zero_diffusivity, &
         config_use_sedimentation, &
         config_use_fast_sedimentation, &
         config_use_positive_sedimentation, &
         config_use_bgc_sinking_flux, &
         config_use_step_initial_profiles, &
         config_use_benthic_reactions, &
         config_use_flux_correction, &
         config_use_ocean_conc, &
         config_use_deep_source, &
         config_use_constant_diffusivity, &
         config_benthosTracers_diagnostic_fields_level1 !,  &
!         config_benthosTracers_diagnostic_fields_level2,  &
!         config_benthosTracers_diagnostic_fields_level3,  &
!         config_benthosTracers_diagnostic_fields_level4,  &
!         config_benthosTracers_diagnostic_fields_level5

      ! seafloorStorage

!      real (kind=RKIND), dimension(:,:), pointer:: &
!         benthosTracersStorage

      !
      ! level 1 diagnostics
      !

      real (kind=RKIND), dimension(:,:), pointer :: &
         netBenthosReactionTend,   &
         netBenthosTransportTend,  &
         netBenthosTend,           &
         primarySourceTend,        &
         primarySinkTend,          &
         secondarySourceTend,      &
         secondarySinkTend,        &
         carbonateSourceTend,      &
         carbonateSinkTend,        &
         netBenthosTracersCell,    &
         netElements,              &
         netElementsReactionTend,  &
         netElementsTransportTend, &
         netElementsTend

      real (kind=RKIND) :: convertLengthMKStoCGS = 100.0_RKIND
      real (kind=RKIND) :: convertLengthCGStoMKS = 0.01_RKIND
      real (kind=RKIND) :: zTop, zBot

      integer :: iCell, iLevel, iTracer, numColumns, column, nElements

      err = 0

      call mpas_timer_start("benthos source-sink")

      call mpas_pool_get_subpool(forcingPool, 'benthosAuxiliary', benthosAuxiliaryPool)
      call mpas_pool_get_subpool(forcingPool, 'seafloorStorage', seafloorStoragePool)
      call mpas_pool_get_subpool(forcingPool, 'benthosDiagFieldsLevel1', benthosDiagFieldsLevel1Pool)

 !     call mpas_pool_get_array(seafloorStorage, 'benthosTracersStorage', benthosTracersStorage)
 !     call mpas_pool_get_array(seafloorStorage, 'benthosPhPrev3d', benthosPhPrev3d)

      call mpas_pool_get_config(ocnConfigs, 'config_set_benthos_flags',  config_set_benthos_flags)
      call mpas_pool_get_config(ocnConfigs, 'config_use_carbonate_saturation',  config_use_carbonate_saturation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_secondary_reactions',  config_use_secondary_reactions)
      call mpas_pool_get_config(ocnConfigs, 'config_use_depth_dependent_porosity',  config_use_depth_dependent_porosity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_non_zero_diffusivity',  config_use_non_zero_diffusivity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_sedimentation',  config_use_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_fast_sedimentation',  config_use_fast_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_positive_sedimentation',  config_use_positive_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_bgc_sinking_flux',  config_use_bgc_sinking_flux)
      call mpas_pool_get_config(ocnConfigs, 'config_use_step_initial_profiles',  config_use_step_initial_profiles)
      call mpas_pool_get_config(ocnConfigs, 'config_use_benthic_reactions',  config_use_benthic_reactions)
      call mpas_pool_get_config(ocnConfigs, 'config_use_flux_correction',  config_use_flux_correction)
      call mpas_pool_get_config(ocnConfigs, 'config_use_ocean_conc',  config_use_ocean_conc)
      call mpas_pool_get_config(ocnConfigs, 'config_use_deep_source',  config_use_deep_source)
      call mpas_pool_get_config(ocnConfigs, 'config_use_constant_diffusivity',  config_use_constant_diffusivity)
      call mpas_pool_get_config(ocnConfigs, 'config_benthosTracers_diagnostic_fields_level1', &
              config_benthosTracers_diagnostic_fields_level1)

      if (config_benthosTracers_diagnostic_fields_level1) then

         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netBenthosReactionTend', netBenthosReactionTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netBenthosTransportTend', netBenthosTransportTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netBenthosTend', netBenthosTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'primarySourceTend', primarySourceTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'primarySinkTend', primarySinkTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'secondarySourceTend', secondarySourceTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'secondarySinkTend', secondarySinkTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'carbonateSourceTend', carbonateSourceTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'carbonateSinkTend', carbonateSinkTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netBenthosTracersCell', netBenthosTracersCell)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netElements', netElements)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netElementsReactionTend', netElementsReactionTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netElementsTransportTend', netElementsTransportTend)
         call mpas_pool_get_array(benthosDiagFieldsLevel1Pool, 'netElementsTend', netElementsTend)

         if (associated(netElements))  nElements = size(netElements, dim=1)
      else
         nElements = 7  ! default 
      end if

      numColumns = 1
      column = 1
      !DWJ 08/05/2016: This loop needs OpenMP added to it.
      do iCell=1,nCellsSolve
          benthos_input%bottom_level(column) = maxLevelCell(iCell)
          benthos_input%oceanBottomTemperature(column) = activeTracers(indexTemperature, maxLevelCell(iCell), iCell)
          benthos_input%oceanBottomSalinity(column) = activeTracers(indexSalinity, maxLevelCell(iCell), iCell)
          benthos_input%oceanBottomSilicate(column) = ecosysTracers(ecosysIndices%sio3_ind, maxLevelCell(iCell), iCell)
          benthos_input%deepStorage(column,benthos_indices%poca_ind) = benthosTracersStorage(benthosIndices%poca_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%pocb_ind) = benthosTracersStorage(benthosIndices%pocb_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%pocc_ind) = benthosTracersStorage(benthosIndices%pocc_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%pona_ind) = benthosTracersStorage(benthosIndices%pona_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%ponb_ind) = benthosTracersStorage(benthosIndices%ponb_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%ponc_ind) = benthosTracersStorage(benthosIndices%ponc_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%popa_ind) = benthosTracersStorage(benthosIndices%popa_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%popb_ind) = benthosTracersStorage(benthosIndices%popb_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%popc_ind) = benthosTracersStorage(benthosIndices%popc_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%o2_ind) = benthosTracersStorage(benthosIndices%o2_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%nh4_ind) = benthosTracersStorage(benthosIndices%nh4_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%h2po4_ind) = benthosTracersStorage(benthosIndices%h2po4_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%co2_ind) = benthosTracersStorage(benthosIndices%co2_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%no3_ind) = benthosTracersStorage(benthosIndices%no3_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%mno2a_ind) = benthosTracersStorage(benthosIndices%mno2a_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%mno2b_ind) = benthosTracersStorage(benthosIndices%mno2b_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%mn_ind) = benthosTracersStorage(benthosIndices%mn_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%feoh3a_ind) = benthosTracersStorage(benthosIndices%feoh3a_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%feoh3b_ind) = benthosTracersStorage(benthosIndices%feoh3b_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%fe_ind) = benthosTracersStorage(benthosIndices%fe_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%fepa_ind) = benthosTracersStorage(benthosIndices%fepa_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%fepb_ind) = benthosTracersStorage(benthosIndices%fepb_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%so4_ind) = benthosTracersStorage(benthosIndices%so4_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%h2s_ind) = benthosTracersStorage(benthosIndices%h2s_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%ch4_ind) = benthosTracersStorage(benthosIndices%ch4_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%dic_ind) = benthosTracersStorage(benthosIndices%dic_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%alk_ind) = benthosTracersStorage(benthosIndices%alk_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%hco3_ind) = benthosTracersStorage(benthosIndices%hco3_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%fes_ind) = benthosTracersStorage(benthosIndices%fes_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%fes2_ind) = benthosTracersStorage(benthosIndices%fes2_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%s_ind) = benthosTracersStorage(benthosIndices%s_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%caco3a_ind) = benthosTracersStorage(benthosIndices%caco3a_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%caco3b_ind) = benthosTracersStorage(benthosIndices%caco3b_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%co3_ind) = benthosTracersStorage(benthosIndices%co3_ind,iCell)
          benthos_input%deepStorage(column,benthos_indices%camgco3_ind) = benthosTracersStorage(benthosIndices%camgco3_ind,iCell)

          do iLevel=1,nBenthicVertLevels
             benthos_input%PH_PREV_3D(iLevel,column) = benthosPhPrev3d(iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%poca_ind) = benthosTracers(benthosIndices%poca_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%pocb_ind) = benthosTracers(benthosIndices%pocb_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%pocc_ind) = benthosTracers(benthosIndices%pocc_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%pona_ind) = benthosTracers(benthosIndices%pona_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%ponb_ind) = benthosTracers(benthosIndices%ponb_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%ponc_ind) = benthosTracers(benthosIndices%ponc_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%popa_ind) = benthosTracers(benthosIndices%popa_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%popb_ind) = benthosTracers(benthosIndices%popb_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%popc_ind) = benthosTracers(benthosIndices%popc_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%o2_ind) = benthosTracers(benthosIndices%o2_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%nh4_ind) = benthosTracers(benthosIndices%nh4_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%h2po4_ind) = benthosTracers(benthosIndices%h2po4_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%co2_ind) = benthosTracers(benthosIndices%co2_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%no3_ind) = benthosTracers(benthosIndices%no3_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%mno2a_ind) = benthosTracers(benthosIndices%mno2a_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%mno2b_ind) = benthosTracers(benthosIndices%mno2b_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%mn_ind) = benthosTracers(benthosIndices%mn_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%feoh3a_ind) = benthosTracers(benthosIndices%feoh3a_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%feoh3b_ind) = benthosTracers(benthosIndices%feoh3b_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fe_ind) = benthosTracers(benthosIndices%fe_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fepa_ind) = benthosTracers(benthosIndices%fepa_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fepb_ind) = benthosTracers(benthosIndices%fepb_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%so4_ind) = benthosTracers(benthosIndices%so4_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%h2s_ind) = benthosTracers(benthosIndices%h2s_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%ch4_ind) = benthosTracers(benthosIndices%ch4_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%dic_ind) = benthosTracers(benthosIndices%dic_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%alk_ind) = benthosTracers(benthosIndices%alk_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%hco3_ind) = benthosTracers(benthosIndices%hco3_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fes_ind) = benthosTracers(benthosIndices%fes_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%fes2_ind) = benthosTracers(benthosIndices%fes2_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%s_ind) = benthosTracers(benthosIndices%s_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%caco3a_ind) = benthosTracers(benthosIndices%caco3a_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%caco3b_ind) = benthosTracers(benthosIndices%caco3b_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%co3_ind) = benthosTracers(benthosIndices%co3_ind,iLevel,iCell)
             benthos_input%benthosTracerBulk(iLevel,column,benthos_indices%camgco3_ind) = benthosTracers(benthosIndices%camgco3_ind,iLevel,iCell)
          end do ! iLevel

         call benthos_SourceSink(benthos_input, benthos_forcing,bottomDepthObserved(iCell), &
                             benthos_output, benthos_diagnostic_fields, nBenthicVertLevels, &
                             numBenthicColumnsMax,dt,err)

!         do iLevel=1,maxLevelCell(iCell)

!            ecosysTracersTend(ecosysIndices%po4_ind,iLevel,iCell) =  &
!               ecosysTracersTend(ecosysIndices%po4_ind,iLevel,iCell)   &
!               + benthos_output%BGC_tendencies(iLevel,column,BGC_indices%po4_ind)*layerThickness(iLevel,iCell)
!            ecosysTracersTend(ecosysIndices%no3_ind,iLevel,iCell) =  &
!               ecosysTracersTend(ecosysIndices%no3_ind,iLevel,iCell)   &
!               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%no3_ind)*layerThickness(iLevel,iCell)
!            ecosysTracersTend(ecosysIndices%sio3_ind,iLevel,iCell) =  &
!               ecosysTracersTend(ecosysIndices%sio3_ind,iLevel,iCell)   &
!               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%sio3_ind)*layerThickness(iLevel,iCell)
!            ecosysTracersTend(ecosysIndices%nh4_ind,iLevel,iCell) =  &
!               ecosysTracersTend(ecosysIndices%nh4_ind,iLevel,iCell)   &
!               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%nh4_ind)*layerThickness(iLevel,iCell)
!            ecosysTracersTend(ecosysIndices%fe_ind,iLevel,iCell) =  &
!               ecosysTracersTend(ecosysIndices%fe_ind,iLevel,iCell)   &
!               + BGC_output%BGC_tendencies(iLevel,column,BGC_indices%fe_ind)*layerThickness(iLevel,iCell)
 
            !
            ! level 1 diagnostics
            !

!            if (config_ecosysTracers_diagnostic_fields_level1) then
!               ecosys_diag_PAR_avg(iLevel,iCell) =  &
!                           BGC_diagnostic_fields%diag_PAR_avg(iLevel,column)
!               ecosys_diag_POC_FLUX_IN(iLevel,iCell) =  &
!                           BGC_diagnostic_fields%diag_POC_FLUX_IN(iLevel,column)*convertLengthCGStoMKS
!               ecosys_diag_CaCO3_FLUX_IN(iLevel,iCell) =  &
!                           BGC_diagnostic_fields%diag_CaCO3_FLUX_IN(iLevel,column)*convertLengthCGStoMKS
! 
!            endif  !  config_ecosysTracers_diagnostic_fields_level1

            !
            ! level 2 diagnostics
            !
 !...
         !
         ! level 1 diagnostics with no depth dependence
         !
 !        if (config_ecosysTracers_diagnostic_fields_level1) then
 !           ecosys_diag_photoC_TOT_zint(iCell) =  &
 !                       BGC_diagnostic_fields%diag_photoC_TOT_zint(column)*convertLengthCGStoMKS
 !           ecosys_diag_photoC_NO3_TOT_zint(iCell) =  &
 !                       BGC_diagnostic_fields%diag_photoC_NO3_TOT_zint(column)*convertLengthCGStoMKS
 !           ecosys_diag_O2_ZMIN(iCell) =  &
 !                       BGC_diagnostic_fields%diag_O2_ZMIN(column)
!         endif  !  config_ecosysTracers_diagnostic_fields_level1

         !
         ! level 3 diagnostics with no depth dependence
         !

       end do ! iCell

      call mpas_timer_stop("benthos source-sink")

   !--------------------------------------------------------------------

   end subroutine benthos_tracer_compute!}}}

!***********************************************************************
!
!  routine benthos_tracer_init
!
!> \brief   Initializes
!> \author  Nicole Jeffery
!> \date    07/16/2020
!> \details
!>  This routine initializes indices for benthos tracers
!
!-----------------------------------------------------------------------

   subroutine ocn_benthos_tracer_init(domain,err)!{{{


!NOTE:  called from mpas_ocn_forward_mode.F

      type (domain_type), intent(inout) :: domain !< Input/Output: domain information

      integer, intent(out) :: err !< Output: error flag

      type (mpas_pool_type), pointer :: forcingPool
      type (mpas_pool_type), pointer :: tracersPool
      type (mpas_pool_type), pointer :: seafloorTracersPool
      type (mpas_pool_type), pointer :: seafloorStoragePool
      type (mpas_pool_type), pointer :: benthosTracersPool
      type (mpas_pool_type), pointer :: benthosAuxiliaryPool
      type (mpas_pool_type), pointer :: benthosDiagFieldsLevel1Pool
!      type (mpas_pool_type), pointer :: ecosysMonthlyForcing


      real (kind=RKIND), dimension(:,:,:), pointer :: &
          benthosTracers

      ! seafloorStorage

      real (kind=RKIND), dimension(:,:), pointer:: &
          benthosTracersStorage, &
          benthosPhPrev3d

! input flux components in benthosAuxiliary
      real (kind=RKIND), dimension(:), pointer :: &
          riverFluxSed,        &
          riverFluxPOC

      ! input flux components in benthosMonthlyForcing
      real (kind=RKIND), dimension(:), pointer :: &
          riverFluzSed,    &
          riverFluzPOC

      ! scalars
      integer :: &
         nBenthicTracers,      &
         nElements,            &
         numBenthicColumnsMax, &
         iBenthicTracers,      &
         iVertLevels

      ! scalar pointers
      integer, pointer :: &
          nBenthicVertLevels,       &
          index_dummy,              &
          config_set_benthos_flags

!      character(len=strKIND) :: &
!         forcingIntervalMonthly,  &
!         forcingReferenceTimeMonthly

      logical, pointer :: &
          config_do_restart, &
          config_use_carbonate_saturation, &
          config_use_secondary_reactions, &
          config_use_depth_dependent_porosity, &
          config_use_non_zero_diffusivity, &
          config_use_sedimentation, &
          config_use_fast_sedimentation, &
          config_use_positive_sedimentation, &
          config_use_bgc_sinking_flux, &
          config_use_step_initial_profiles, &
          config_use_benthic_reactions, &
          config_use_flux_correction, &
          config_use_ocean_conc, &
          config_use_deep_source, &
          config_use_constant_diffusivity

      ! local variables

      logical :: &
          test_benthos_flags

      !
      ! get tracers pools
      !

      err = 0

      !
      ! Get tracer group so we can get the number of tracers in it
      !

      call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
      call mpas_pool_get_subpool(domain % blocklist % structs, 'tracers', tracersPool)
      call mpas_pool_get_subpool(forcingPool, 'seafloorTracers', seafloorTracersPool)
      call mpas_pool_get_subpool(forcingPool, 'seafloorStorage', seafloorStoragePool)
      call mpas_pool_get_subpool(forcingPool, 'benthosDiagFieldsLevel1', benthosDiagFieldsLevel1Pool)

      call mpas_pool_get_array(seafloorTracersPool, 'benthosTracers', benthosTracers)
      call mpas_pool_get_array(seafloorTracersPool, 'benthosPhPrev3d', benthosPhPrev3d)
      call mpas_pool_get_array(seafloorStoragePool, 'benthosTracersStorage', benthosTracersStorage)

!      call mpas_log_write('Error test warning, mpas_ocn_tracer_benthos.F',MPAS_LOG_WARN)
!      call mpas_log_write('Check benthos_tracer_cnt: $i',MPAS_LOG_ERR,intArgs=(/benthos_tracer_cnt/))

      if (associated(benthosTracers)) then   ! associated for benthosGRP

      nBenthicTracers = size(benthosTracers, dim=1)
      if (benthos_tracer_cnt /= nBenthicTracers) then
         call mpas_log_write('benthos_tracer_cnt /= nBenthicTracers: $i',MPAS_LOG_CRIT,intArgs=(/nBenthicTracers/))
         err = 1
         return
      end if

      !
      ! pull nBenthicVertLevels out of the mesh structure
      !

      call mpas_pool_get_dimension(domain % blocklist % dimensions, 'nBenthicVertLevels', nBenthicVertLevels)

!-----------------------------------------------------------------------
!  initialize benthos parameters
!-----------------------------------------------------------------------

      allocate( benthos_indices%short_name(benthos_tracer_cnt) )
      allocate( benthos_indices%long_name(benthos_tracer_cnt) )
      allocate( benthos_indices%units(benthos_tracer_cnt) )

! no need to allocate the above fields for ecosysIndices (only use the index)

!-----------------------------------------------------------------------
!  sets most of benthos parameters
!  sets namelist defaults
!-----------------------------------------------------------------------

      call benthos_parms_init   ! initialize namelist defaults
      !
      ! for now only do 1 column at a time
      !
      T0_Kelvin_benthos = T0_Kelvin

      numBenthicColumnsMax = 1

      benthos_indices%poca_ind    =  1
      benthos_indices%pocb_ind    =  2
      benthos_indices%pocc_ind    =  3
      benthos_indices%pona_ind    =  4
      benthos_indices%ponb_ind    =  5
      benthos_indices%ponc_ind    =  6
      benthos_indices%popa_ind    =  7
      benthos_indices%popb_ind    =  8
      benthos_indices%popc_ind    =  9
      benthos_indices%o2_ind      = 10
      benthos_indices%nh4_ind     = 11
      benthos_indices%h2po4_ind   = 12
      benthos_indices%co2_ind     = 13
      benthos_indices%no3_ind     = 14
      benthos_indices%mno2a_ind   = 15
      benthos_indices%mno2b_ind   = 16
      benthos_indices%mn_ind      = 17
      benthos_indices%feoh3a_ind  = 18
      benthos_indices%feoh3b_ind  = 19
      benthos_indices%fe_ind      = 20
      benthos_indices%fepa_ind    = 21
      benthos_indices%fepb_ind    = 22
      benthos_indices%so4_ind     = 23
      benthos_indices%h2s_ind     = 24
      benthos_indices%ch4_ind     = 25
      benthos_indices%dic_ind     = 26
      benthos_indices%alk_ind     = 27
      benthos_indices%hco3_ind    = 28
      benthos_indices%fes_ind     = 29
      benthos_indices%fes2_ind    = 30
      benthos_indices%s_ind       = 31
      benthos_indices%caco3a_ind  = 32
      benthos_indices%caco3b_ind  = 33
      benthos_indices%co3_ind     = 34
      benthos_indices%camgco3_ind = 35

      benthos_element_indices%carbon_ind     = 1
      benthos_element_indices%oxygen_ind     = 2
      benthos_element_indices%nitrogen_ind   = 3
      benthos_element_indices%phosphorus_ind = 4
      benthos_element_indices%sulfur_ind     = 5
      benthos_element_indices%manganese_ind  = 6
      benthos_element_indices%iron_ind       = 7

      nElements = 7
!      call mpas_log_write('test assign indices: benthos_indices%caco3a_ind = $i',MPAS_LOG_ERR,intArgs=(/benthos_indices%caco3a_ind/))
!      call mpas_log_write('total number of benthos tracers = $i',MPAS_LOG_ERR,intArgs=(/benthos_tracer_cnt/))
!      err = 1
!      return

      call mpas_pool_get_dimension(seafloorTracersPool, 'index_pocaBenthos',  index_dummy)
      benthosIndices%poca_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_pocbBenthos',  index_dummy)
      benthosIndices%pocb_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_poccBenthos',  index_dummy)
      benthosIndices%pocc_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_ponaBenthos',  index_dummy)
      benthosIndices%pona_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_ponbBenthos',  index_dummy)
      benthosIndices%ponb_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_poncBenthos',  index_dummy)
      benthosIndices%ponc_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_popaBenthos',  index_dummy)
      benthosIndices%popa_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_popbBenthos',  index_dummy)
      benthosIndices%popb_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_popcBenthos',  index_dummy)
      benthosIndices%popc_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_o2Benthos',  index_dummy)
      benthosIndices%o2_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_nh4Benthos',  index_dummy)
      benthosIndices%nh4_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_h2po4Benthos',  index_dummy)
      benthosIndices%h2po4_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_co2Benthos',  index_dummy)
      benthosIndices%co2_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_no3Benthos',  index_dummy)
      benthosIndices%no3_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_mno2aBenthos',  index_dummy)
      benthosIndices%mno2a_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_mno2bBenthos',  index_dummy)
      benthosIndices%mno2b_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_mnBenthos',  index_dummy)
      benthosIndices%mn_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_feoh3aBenthos',  index_dummy)
      benthosIndices%feoh3a_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_feoh3bBenthos',  index_dummy)
      benthosIndices%feoh3b_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_feBenthos',  index_dummy)
      benthosIndices%fe_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_fepaBenthos',  index_dummy)
      benthosIndices%fepa_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_fepbBenthos',  index_dummy)
      benthosIndices%fepb_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_so4Benthos',  index_dummy)
      benthosIndices%so4_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_h2sBenthos',  index_dummy)
      benthosIndices%h2s_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_ch4Benthos',  index_dummy)
      benthosIndices%ch4_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_dicBenthos',  index_dummy)
      benthosIndices%dic_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_alkBenthos',  index_dummy)
      benthosIndices%alk_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_hco3Benthos',  index_dummy)
      benthosIndices%hco3_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_fesBenthos',  index_dummy)
      benthosIndices%fes_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_fes2Benthos',  index_dummy)
      benthosIndices%fes2_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_sBenthos',  index_dummy)
      benthosIndices%s_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_caco3aBenthos',  index_dummy)
      benthosIndices%caco3a_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_caco3bBenthos',  index_dummy)
      benthosIndices%caco3b_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_co3Benthos',  index_dummy)
      benthosIndices%co3_ind  = index_dummy
      call mpas_pool_get_dimension(seafloorTracersPool, 'index_camgco3Benthos',  index_dummy)
      benthosIndices%camgco3_ind  = index_dummy

      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netCarbon',  index_dummy)
      elementIndices%carbon_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netOxygen',  index_dummy)
      elementIndices%oxygen_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netNitrogen',  index_dummy)
      elementIndices%nitrogen_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netPhosphorus',  index_dummy)
      elementIndices%phosphorus_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netSulfur',  index_dummy)
      elementIndices%sulfur_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netManganese',  index_dummy)
      elementIndices%manganese_ind  = index_dummy
      call mpas_pool_get_dimension(benthosDiagFieldsLevel1Pool, 'index_netIron',  index_dummy)
      elementIndices%iron_ind  = index_dummy

!      call mpas_pool_get_dimension(tracersPool, 'index_DMS',  index_dummy)
!      dmsIndices%dms_ind  = index_dummy
!      call mpas_pool_get_dimension(tracersPool, 'index_DMSP', index_dummy)
!      dmsIndices%dmsp_ind = index_dummy

      call mpas_pool_get_dimension(tracersPool, 'index_PO4',      index_dummy)
      ecosysIndices%po4_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_NO3',      index_dummy)
      ecosysIndices%no3_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_SiO3',     index_dummy)
      ecosysIndices%sio3_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_NH4',      index_dummy)
      ecosysIndices%nh4_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_Fe',       index_dummy)
      ecosysIndices%fe_ind              = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_O2',       index_dummy)
      ecosysIndices%o2_ind              = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DIC',      index_dummy)
      ecosysIndices%dic_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DIC_ALT_CO2', index_dummy)
      ecosysIndices%dic_alt_co2_ind     = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_ALK',      index_dummy)
      ecosysIndices%alk_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DOC',      index_dummy)
      ecosysIndices%doc_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DON',      index_dummy)
      ecosysIndices%don_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DOFe',     index_dummy)
      ecosysIndices%dofe_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DOP',      index_dummy)
      ecosysIndices%dop_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DOPr',     index_dummy)
      ecosysIndices%dopr_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_DONr',     index_dummy)
      ecosysIndices%donr_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_zooC',     index_dummy)
      ecosysIndices%zooC_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spChl',    index_dummy)
      ecosysIndices%spChl_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spC',      index_dummy)
      ecosysIndices%spC_ind             = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spFe',     index_dummy)
      ecosysIndices%spFe_ind            = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_spCaCO3',  index_dummy)
      ecosysIndices%spCaCO3_ind         = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatChl',  index_dummy)
      ecosysIndices%diatChl_ind         = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatC',    index_dummy)
      ecosysIndices%diatC_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatFe',   index_dummy)
      ecosysIndices%diatFe_ind          = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diatSi',   index_dummy)
      ecosysIndices%diatSi_ind          = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diazChl',  index_dummy)
      ecosysIndices%diazChl_ind         = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diazC',    index_dummy)
      ecosysIndices%diazC_ind           = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_diazFe',   index_dummy)
      ecosysIndices%diazFe_ind          = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_phaeoChl', index_dummy)
      ecosysIndices%phaeoChl_ind        = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_phaeoC',   index_dummy)
      ecosysIndices%phaeoC_ind          = index_dummy
      call mpas_pool_get_dimension(tracersPool, 'index_phaeoFe',  index_dummy)
      ecosysIndices%phaeoFe_ind         = index_dummy

      call mpas_pool_get_config(ocnConfigs, 'config_set_benthos_flags',  config_set_benthos_flags)
      call mpas_pool_get_config(ocnConfigs, 'config_use_carbonate_saturation',  config_use_carbonate_saturation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_secondary_reactions',  config_use_secondary_reactions)
      call mpas_pool_get_config(ocnConfigs, 'config_use_depth_dependent_porosity',  config_use_depth_dependent_porosity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_non_zero_diffusivity',  config_use_non_zero_diffusivity)
      call mpas_pool_get_config(ocnConfigs, 'config_use_sedimentation',  config_use_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_fast_sedimentation',  config_use_fast_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_positive_sedimentation',  config_use_positive_sedimentation)
      call mpas_pool_get_config(ocnConfigs, 'config_use_bgc_sinking_flux',  config_use_bgc_sinking_flux)
      call mpas_pool_get_config(ocnConfigs, 'config_use_step_initial_profiles',  config_use_step_initial_profiles)
      call mpas_pool_get_config(ocnConfigs, 'config_use_benthic_reactions',  config_use_benthic_reactions)
      call mpas_pool_get_config(ocnConfigs, 'config_use_flux_correction',  config_use_flux_correction)
      call mpas_pool_get_config(ocnConfigs, 'config_use_ocean_conc',  config_use_ocean_conc)
      call mpas_pool_get_config(ocnConfigs, 'config_use_deep_source',  config_use_deep_source)
      call mpas_pool_get_config(ocnConfigs, 'config_use_constant_diffusivity',  config_use_constant_diffusivity)

    ! Set configuration flags

      call benthos_flags_init(config_set_benthos_flags, config_use_carbonate_saturation, config_use_secondary_reactions, &
           config_use_depth_dependent_porosity, config_use_non_zero_diffusivity, config_use_sedimentation, &
           config_use_fast_sedimentation, config_use_positive_sedimentation, config_use_bgc_sinking_flux, config_use_step_initial_profiles, &
           config_use_benthic_reactions, config_use_flux_correction, config_use_ocean_conc, &
           config_use_deep_source, config_use_constant_diffusivity)

!NJ-TEST
test_benthos_flags = .true.
      if (test_benthos_flags) then
         call mpas_log_write('test benthos_flags_init')
         call mpas_log_write('config_set_benthos_flags: $i',MPAS_LOG_CRIT,intArgs=(/config_set_benthos_flags/))
         if (config_use_carbonate_saturation) then
            call mpas_log_write('config_use_carbonate_saturation = true')
         else
            call mpas_log_write('config_use_carbonate_saturation = false')
         end if
         if (config_use_secondary_reactions) then
            call mpas_log_write('config_use_secondary_reactions = true')
         else
            call mpas_log_write('config_use_secondary_reactions = false')
         end if
         if (config_use_depth_dependent_porosity) then
            call mpas_log_write('config_use_depth_dependent_porosity = true')
         else
            call mpas_log_write('config_use_depth_dependent_porosity = false')
         end if
         if (config_use_non_zero_diffusivity) then
            call mpas_log_write('config_use_non_zero_diffusivity = true')
         else
            call mpas_log_write('config_use_non_zero_diffusivity = false')
         end if
         if (config_use_sedimentation) then
            call mpas_log_write('config_use_sedimentation = true')
         else
            call mpas_log_write('config_use_sedimentation = false')
         end if
         if (config_use_fast_sedimentation) then
            call mpas_log_write('config_use_fast_sedimentation = true')
         else
            call mpas_log_write('config_use_fast_sedimentation = false')
         end if
         if (config_use_positive_sedimentation) then
            call mpas_log_write('config_use_positive_sedimentation = true')
         else
            call mpas_log_write('config_use_positive_sedimentation = false')
         end if
         if (config_use_bgc_sinking_flux) then
            call mpas_log_write('config_use_bgc_sinking_flux = true')
         else
            call mpas_log_write('config_use_bgc_sinking_flux = false')
         end if
         if (config_use_step_initial_profiles) then
            call mpas_log_write('config_use_step_initial_profiles = true')
         else
            call mpas_log_write('config_use_step_initial_profiles = false')
         end if
         if (config_use_benthic_reactions) then
            call mpas_log_write('config_use_benthic_reactions = true')
         else
            call mpas_log_write('config_use_benthic_reactions = false')
         end if
         if (config_use_flux_correction) then
            call mpas_log_write('config_use_flux_correction = true')
         else
            call mpas_log_write('config_use_flux_correction = false')
         end if
         if (config_use_ocean_conc) then
            call mpas_log_write('config_use_ocean_conc = true')
         else
            call mpas_log_write('config_use_ocean_conc = false')
         end if
         if (config_use_deep_source) then
            call mpas_log_write('config_use_deep_source = true')
         else
            call mpas_log_write('config_use_deep_source = false')
         end if
         if (config_use_constant_diffusivity) then
            call mpas_log_write('config_use_constant_diffusivity = true')
         else
            call mpas_log_write('config_use_constant_diffusivity = false')
         end if
         err = 1
         return
     end if ! test_benthos_flags

     ! benthos_init sets short and long names, units in benthos_indices
     ! also sets benthos model parameter fields

      call benthos_init(benthos_indices, benthos_element_indices, nBenthicVertLevels, nBenthicTracers, nElements)

!NOTES:

!also check short_name with mpas variable name

!-----------------------------------------------------------------------
!  allocate input, forcing, diagnostic arrays
!-----------------------------------------------------------------------

      allocate ( benthos_input%benthosTracerBulk(nBenthicVertLevels+1, numBenthicColumnsMax, benthos_tracer_cnt) )
!      allocate ( benthos_input%tracerConc(nBenthicVertLevels+1, numBenthicColumnsMax, benthos_tracer_cnt))
      allocate ( benthos_input%deepStorage(numBenthicColumnsMax, benthos_tracer_cnt))
      allocate ( benthos_input%oceanBottomDepth(numBenthicColumnsMax))
      allocate ( benthos_input%tracerType(benthos_tracer_cnt) )
      allocate ( benthos_input%PH_PREV_3D(nBenthicVertLevels+1, numBenthicColumnsMax))
      allocate ( benthos_input%bottom_level(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomTemperature(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomSalinity(numBenthicColumnsMax))
      allocate ( benthos_input%oceanBottomSilicate(numBenthicColumnsMax))
      allocate ( benthos_output%benthosTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate ( benthos_output%benthosReactionTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate ( benthos_output%benthosTransportTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, benthos_tracer_cnt) )

      allocate ( benthos_output%elementTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, element_tracer_cnt) )
      allocate ( benthos_output%elementReactionTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, element_tracer_cnt) )
      allocate ( benthos_output%elementTransportTendencies(nBenthicVertLevels+1, numBenthicColumnsMax, element_tracer_cnt) )
      allocate ( benthos_output%elementDeepStorage(numBenthicColumnsMax, element_tracer_cnt) )

      allocate ( benthos_forcing%deepStoragePrescribedFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate ( benthos_forcing%benthosOceanSedimentFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
        !---------------------------------------------------------------------------
        !   allocate flux diagnostic output fields
        !   all fluxes are positive down
        !---------------------------------------------------------------------------

      allocate (benthos_flux_diagnostic_fields%deepStorageDiffVelFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%deepStorageFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%benthosOceanDiffVelFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%benthosOceanFlux(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%elementSedimentation(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%elementBurialExchange(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_flux_diagnostic_fields%elementOceanExchange(numBenthicColumnsMax, element_tracer_cnt) )

      !---------------------------------------------------------------------------
      !   allocate diagnostic output fields
      !---------------------------------------------------------------------------


      ! 2D column integrated diagnostics

      allocate (benthos_diagnostic_fields%diag_primarySourceTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_primarySinkTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_secondarySourceTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_secondarySinkTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_carbonateSourceTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_carbonateSinkTend(numBenthicColumnsMax, benthos_tracer_cnt) )

      allocate (benthos_diagnostic_fields%diag_netElements(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netBenthosTracersCell(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netBenthosReactionTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netBenthosTransportTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netBenthosTend(numBenthicColumnsMax, benthos_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netElementsReactionTend(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netElementsTransportTend(numBenthicColumnsMax, element_tracer_cnt) )
      allocate (benthos_diagnostic_fields%diag_netElementsTend(numBenthicColumnsMax, element_tracer_cnt) )

! initialize monthly forcing to be read from file (sediment fluxes)

!      forcingIntervalMonthly = "0000-01-00_00:00:00"
!      forcingReferenceTimeMonthly = "0000-01-15_00:00:00"

      call MPAS_pool_get_config(domain % configs, 'config_do_restart', config_do_restart)

!      call MPAS_forcing_init_group( forcingGroupHead,  &
!                "benthosMonthlyClimatology", &
!                domain, &
!                '0000-01-01_00:00:00', &
!                '0000-01-01_00:00:00', &
!                '0001-00-00_00:00:00', &
!                config_do_restart)

      call mpas_pool_get_subpool(forcingPool, 'benthosAuxiliary', benthosAuxiliaryPool)
      call mpas_pool_get_array(benthosAuxiliaryPool,'riverFluxSed', riverFluxSed)
      call mpas_pool_get_array(benthosAuxiliaryPool,'riverFluxPOC', riverFluxPOC)

!      call mpas_pool_get_subpool(domain % blocklist % structs, 'benthosMonthlyForcing', benthosMonthlyForcing)

!      call mpas_pool_get_array(benthosMonthlyForcing, 'riverFluzSed', riverFluzSed)
!      call mpas_pool_get_array(benthosMonthlyForcing, 'riverFluzPOC', riverFluzPOC)

!      call MPAS_forcing_init_field( domain % streamManager, &
!                forcingGroupHead, &
!                'benthosMonthlyClimatology', &
!                'riverFluzSed', &
!                'benthos_monthly_flux', &
!                'benthosMonthlyForcing',  &
!                'riverFluxSed',  &
!                'linear',  &
!                forcingReferenceTimeMonthly,  &
!                forcingIntervalMonthly)

!      call MPAS_forcing_init_field( domain % streamManager, &
!                forcingGroupHead, &
!                'benthosMonthlyClimatology', &
!                'riverFluzPOC', &
!                'benthos_monthly_flux', &
!                'benthosMonthlyForcing',  &
!                'riverFluxPOC',  &
!                'linear',  &
!                forcingReferenceTimeMonthly,  &
!                forcingIntervalMonthly)


   ! TracerType: 1 for solutes, 0 for solids

      benthos_input%tracerType(:) = 1

   ! solids
      benthos_input%tracerType(benthos_indices%poca_ind) = 0
      benthos_input%tracerType(benthos_indices%pocb_ind) = 0
      benthos_input%tracerType(benthos_indices%pocc_ind) = 0
      benthos_input%tracerType(benthos_indices%pona_ind) = 0
      benthos_input%tracerType(benthos_indices%ponb_ind) = 0
      benthos_input%tracerType(benthos_indices%ponc_ind) = 0
      benthos_input%tracerType(benthos_indices%popa_ind) = 0
      benthos_input%tracerType(benthos_indices%popb_ind) = 0
      benthos_input%tracerType(benthos_indices%popc_ind) = 0
      benthos_input%tracerType(benthos_indices%mno2a_ind) = 0
      benthos_input%tracerType(benthos_indices%mno2b_ind) = 0
      benthos_input%tracerType(benthos_indices%feoh3a_ind) = 0
      benthos_input%tracerType(benthos_indices%feoh3b_ind) = 0
      benthos_input%tracerType(benthos_indices%fepa_ind) = 0
      benthos_input%tracerType(benthos_indices%fepb_ind) = 0
      benthos_input%tracerType(benthos_indices%caco3a_ind) = 0
      benthos_input%tracerType(benthos_indices%caco3b_ind) = 0
      benthos_input%tracerType(benthos_indices%camgco3_ind) = 0
      benthos_input%tracerType(benthos_indices%s_ind) = 0
      benthos_input%tracerType(benthos_indices%fes_ind) = 0
      benthos_input%tracerType(benthos_indices%fes2_ind) = 0

   ! prescribed initial concentration array -- modify later

      if (.not. config_do_restart) then

         benthosPhPrev3d(:,:) = 0.0_RKIND
         benthosTracersStorage(:,:) = 0.0_RKIND
         benthosTracers(benthos_indices%poca_ind,:,:) = 1060.0_RKIND ! *0.49 !  mmol/kg
         benthosTracers(benthos_indices%pocb_ind,:,:) = 1060.0_RKIND*0.15_RKIND !  mmol/kg
         benthosTracers(benthos_indices%pocc_ind,:,:) = 1060.0_RKIND*0.36_RKIND !  mmol/kg
         benthosTracers(benthos_indices%popa_ind,:,:) = benthosTracers(benthos_indices%poca_ind,:,:)*CtoP !
         benthosTracers(benthos_indices%popb_ind,:,:) = benthosTracers(benthos_indices%pocb_ind,:,:)*CtoP !
         benthosTracers(benthos_indices%popc_ind,:,:) = benthosTracers(benthos_indices%pocc_ind,:,:)*CtoP !
         benthosTracers(benthos_indices%pona_ind,:,:) = benthosTracers(benthos_indices%popa_ind,:,:)*CtoP/NtoP !
         benthosTracers(benthos_indices%ponb_ind,:,:) = benthosTracers(benthos_indices%popb_ind,:,:)*CtoP/NtoP !
         benthosTracers(benthos_indices%ponc_ind,:,:) = benthosTracers(benthos_indices%popc_ind,:,:)*CtoP/NtoP !
         benthosTracers(benthos_indices%mno2a_ind,:,:) = 0.1_RKIND !        mmol/kg  Reed et al
         benthosTracers(benthos_indices%feoh3a_ind,:,:) = 1.0_RKIND !    mmol/kg
         benthosTracers(benthos_indices%h2po4_ind,:,:) = 90.0_RKIND !   umol/L  Reed
         benthosTracers(benthos_indices%co2_ind,:,:) = 0.2_RKIND !   mM Krumins
         benthosTracers(benthos_indices%o2_ind,:,:) = 0.0_RKIND  !
         benthosTracers(benthos_indices%mno2a_ind,:,:) = 0.2_RKIND !  mmol/kg  Reed et al
         benthosTracers(benthos_indices%mn_ind,:,:) = 10.0_RKIND !  umol/L Reed
         benthosTracers(benthos_indices%fepa_ind,:,:) = 0.2_RKIND !  mmol/kg  Reed
         benthosTracers(benthos_indices%so4_ind,:,:) = 12.0_RKIND*1000.0_RKIND !  umol/L from mmol/L Reed
         benthosTracers(benthos_indices%h2s_ind,:,:) = 0.0_RKIND !   gets to 600 umol/L at bottom but theres a bottom flux
         benthosTracers(benthos_indices%hco3_ind,:,:) = 20.0_RKIND !  mM  Krumins et al
         benthosTracers(benthos_indices%s_ind,:,:) = 0.2_RKIND !  2000.0 ! 200.0 !  mmol/kg  Reed et al
         benthosTracers(benthos_indices%fepa_ind,:,:) =  0.5_RKIND !     mmol/kg
         benthosTracers(benthos_indices%co3_ind,:,:) = 66.75_RKIND !  66.75 ! 65 umol/kg * 1.027 g/m3   Sulpis et al 2018
         benthosTracers(benthos_indices%dic_ind,:,:) = benthosTracers(benthos_indices%co3_ind,:,:) + benthosTracers(benthos_indices%co2_ind,:,:) + benthosTracers(benthos_indices%hco3_ind,:,:) !
         benthosTracers(benthos_indices%alk_ind,:,:) = 2306.0_RKIND !  87 of the ocean conc. (:,:,benthos_indices%mmol/m3) same as ocean conc. Krumins et al
         benthosTracers(benthos_indices%caco3a_ind,:,:) = 200.0_RKIND*.63_RKIND !  ratio of PIC is 63:13:24...  Roughly  POC:PIC = (:,:,benthos_indices%depth)^(:,:,benthos_indices%-0.56)*64.3   so about 4.8:1
         benthosTracers(benthos_indices%caco3b_ind,:,:) = 200.0_RKIND*0.13_RKIND !
         benthosTracers(benthos_indices%camgco3_ind,:,:) = 200.0_RKIND*0.24_RKIND !
         if (config_use_step_initial_profiles) then
            do iBenthicTracers = 1,nBenthicTracers
               do iVertLevels = 1, int((nBenthicVertLevels+1)/3.0_RKIND)
                  benthosTracers(iBenthicTracers,iVertLevels,:) = benthosTracers(iBenthicTracers,iVertLevels,:)*0.1_RKIND
               end do
               do iVertLevels = int((nBenthicVertLevels+1)/3.0_RKIND*c2_benthos),nBenthicVertLevels+1
                  benthosTracers(iBenthicTracers,iVertLevels,:) = benthosTracers(iBenthicTracers,iVertLevels,:)*0.1_RKIND
               end do
            end do
         end if
      end if ! config_use_step_initial_profiles
      else
         call mpas_log_write('not associated benthosTracers', MPAS_LOG_CRIT)
      end if  !  associated(benthosTracers)

   !--------------------------------------------------------------------

   end subroutine ocn_benthos_tracer_init!}}}

!***********************************************************************

!***********************************************************************
!
!  routine get_benthosData
!
!> \brief   retrieve data needed to compute benthos deposition and runoff inputs
!> \author  Mathew Maltrud
!> \date    03/07/16
!> \details
!>  This routine calls mpas_forcing routines to acquire needed benthos forcing data and interpolates
!>    between time levels.  directly copied from ocn_get_shortWaveData.
!
!-----------------------------------------------------------------------

    subroutine ocn_get_benthos_data( streamManager, &
        domain, &
        simulationClock, &
        firstTimeStep) !{{{

        type (MPAS_streamManager_type), intent(inout) :: streamManager

        type (domain_type) :: domain
        type (MPAS_timeInterval_type) :: timeStepEco
        type (MPAS_clock_type) :: simulationClock

        logical, intent(in) :: firstTimeStep
        character(len=strKind), pointer :: config_dt
        real(kind=RKIND) :: dt

        type (mpas_pool_type), pointer :: forcingPool
        type (mpas_pool_type), pointer :: meshPool
        type (mpas_pool_type), pointer :: benthosAuxiliary
        type (mpas_pool_type), pointer :: benthosMonthlyForcing

        real (kind=RKIND), dimension(:), pointer :: &
         riverFluxNO3

! input flux components in benthosMonthlyForcing
        real (kind=RKIND), dimension(:), pointer :: &
         riverFluzNO3

        integer, pointer :: nCells
        integer :: iCell

        logical, pointer :: config_do_restart

        call MPAS_pool_get_config(domain%configs, 'config_dt', config_dt)
        call MPAS_pool_get_config(domain%configs, 'config_do_restart', config_do_restart)

        call mpas_set_timeInterval(timeStepEco,timeString=config_dt)
        call mpas_get_timeInterval(timeStepEco,dt=dt)

!maltrud debug
        if (firstTimeStep .and. config_do_restart) then
!          call MPAS_forcing_get_forcing(forcingGroupHead, &
!             'benthosMonthlyClimatology', streamManager, 0.0_RKIND)
        else
!          call MPAS_forcing_get_forcing(forcingGroupHead, &
!             'benthosMonthlyClimatology', streamManager, dt)
        end if

        call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'forcing', forcingPool)
        call mpas_pool_get_subpool(forcingPool, 'benthosAuxiliary', benthosAuxiliary)
        call mpas_pool_get_subpool(domain % blocklist % structs, 'benthosMonthlyForcing', benthosMonthlyForcing)

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(benthosAuxiliary, 'riverFluxNO3', riverFluxNO3)

        call mpas_pool_get_array(benthosMonthlyForcing, 'riverFluzNO3', riverFluzNO3)

        do iCell = 1, nCells
           riverFluxNO3(iCell)      = riverFluzNO3(iCell)
        end do

    end subroutine ocn_get_benthos_data!}}}

!***********************************************************************
!
!  routine ocn_benthos_forcing_write_restart
!
!> \brief   writes restart timestamp for benthos data to be read in on future restart
!> \author  Mathew Maltrud
!> \date    03/07/2016

!
!-----------------------------------------------------------------------

   subroutine ocn_benthos_forcing_write_restart(domain)!{{{

      type(domain_type) :: domain

!      call MPAS_forcing_write_restart_times(forcingGroupHead)

    end subroutine ocn_benthos_forcing_write_restart!}}}

end module ocn_tracer_benthos

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
