# Generated automatically by fu.pl
 
SRCS =	m_alloc_mods.f90 mo_aerosols.f90 mo_alai.f90 mo_avgrat.f90 \
	mo_base_types.f90 mo_buffer_fft.f90 mo_call_trans.f90 mo_cfc.f90 \
	mo_column.f90 mo_constants.f90 mo_control.f90 mo_convect_tables.f90 \
	mo_couple.f90 mo_cumulus_flux.f90 mo_decomposition.f90 \
	mo_diag_tendency.f90 mo_diagnostics.f90 mo_diagnostics_zonal.f90 \
	mo_diff.f90 mo_doctor.f90 mo_dsst.f90 mo_exception.f90 mo_f.f90 \
	mo_fft.f90 mo_field.f90 mo_filename.f90 mo_forecast_switches.f90 \
	mo_gaussgrid.f90 mo_global_op.f90 mo_grib.f90 mo_grid.f90 \
	mo_hdiff.f90 mo_hyb.f90 mo_io.f90 mo_io_tables.f90 mo_julian.f90 \
	mo_kind.f90 mo_landsea.f90 mo_legendre.f90 mo_linked_list.f90 \
	mo_longwave.f90 mo_machine.f90 mo_memory_base.f90 mo_memory_f.f90 \
	mo_memory_g1a.f90 mo_memory_g1b.f90 mo_memory_g2a.f90 \
	mo_memory_g2b.f90 mo_memory_g3a.f90 mo_memory_g3b.f90 \
	mo_memory_gl.f90 mo_memory_ls.f90 mo_memory_sp.f90 mo_midatm.f90 \
	mo_mpi.f90 mo_netCDF.f90 mo_nmi.f90 mo_nudging.f90 \
	mo_nudging_buffer.f90 mo_ozone.f90 mo_parallel_slt.f90 \
	mo_param_switches.f90 mo_parameters.f90 mo_physc1.f90 mo_physc2.f90 \
	mo_post.f90 mo_rad1.f90 mo_rad2.f90 mo_rad_switches.f90 \
	mo_radiation.f90 mo_radint.f90 mo_sc1.f90 mo_scan_buffer.f90 \
	mo_scanb.f90 mo_semi_impl.f90 mo_shortwave.f90 mo_skintem.f90 \
	mo_slt.f90 mo_soil_impl.f90 mo_sst.f90 mo_start_dataset.f90 \
	mo_stat_global.f90 mo_stat_zonal.f90 mo_stratiform.f90 \
	mo_test_trans.f90 mo_time_control.f90 mo_timeint.f90 \
	mo_tmp_buffer.f90 mo_tracer.f90 mo_transpose.f90 mo_truncation.f90 \
	mo_vegetation.f90 mo_year.f90

OBJS := $(SRCS:.f90=.o)

all: $(OBJS)

clean:
	rm -f $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

%.o: %.f90
	$(F90) $(F90FLAGS) -c $<

m_alloc_mods.o: mo_control.o mo_diff.o mo_hdiff.o mo_hyb.o mo_physc2.o \
	mo_post.o mo_rad2.o mo_semi_impl.o mo_stat_global.o mo_stat_zonal.o \
	mo_truncation.o
mo_alai.o: mo_control.o mo_decomposition.o mo_doctor.o mo_exception.o mo_io.o \
	mo_memory_g3a.o mo_mpi.o mo_rad_switches.o mo_timeint.o \
	mo_transpose.o
mo_avgrat.o: mo_control.o mo_decomposition.o mo_doctor.o mo_exception.o \
	mo_io.o mo_memory_g3a.o mo_mpi.o mo_rad_switches.o mo_timeint.o \
	mo_transpose.o
mo_buffer_fft.o: mo_decomposition.o
mo_call_trans.o: mo_buffer_fft.o mo_decomposition.o mo_doctor.o \
	mo_linked_list.o mo_memory_f.o mo_memory_g1a.o mo_memory_ls.o \
	mo_memory_sp.o mo_mpi.o mo_sc1.o mo_scan_buffer.o mo_test_trans.o \
	mo_transpose.o
mo_column.o: mo_control.o mo_decomposition.o mo_doctor.o mo_exception.o \
	mo_gaussgrid.o mo_memory_g1a.o mo_memory_g2a.o mo_memory_gl.o \
	mo_mpi.o mo_sc1.o mo_scan_buffer.o mo_start_dataset.o mo_transpose.o \
	$(INCLUDE)/columnctl.inc
mo_constants.o: mo_start_dataset.o
mo_convect_tables.o: mo_constants.o mo_control.o mo_exception.o
mo_couple.o: mo_constants.o mo_control.o mo_decomposition.o mo_doctor.o \
	mo_exception.o mo_filename.o mo_io.o mo_kind.o mo_landsea.o \
	mo_memory_g3a.o mo_mpi.o mo_nudging.o mo_start_dataset.o \
	mo_time_control.o mo_transpose.o mo_year.o
mo_cumulus_flux.o: mo_control.o
mo_decomposition.o: mo_doctor.o mo_exception.o mo_mpi.o
mo_diag_tendency.o: mo_control.o mo_doctor.o mo_exception.o mo_fft.o \
	mo_grib.o mo_legendre.o mo_memory_sp.o mo_post.o mo_start_dataset.o \
	mo_time_control.o mo_truncation.o mo_year.o
mo_diagnostics_zonal.o: mo_parameters.o
mo_diff.o: mo_parameters.o
mo_dsst.o: mo_control.o mo_couple.o mo_decomposition.o mo_doctor.o \
	mo_exception.o mo_kind.o mo_memory_g3a.o mo_mpi.o mo_transpose.o
mo_exception.o: mo_control.o mo_doctor.o mo_mpi.o mo_start_dataset.o
mo_fft.o: mo_doctor.o mo_mpi.o
mo_gaussgrid.o: mo_constants.o mo_control.o mo_parameters.o
mo_global_op.o: mo_decomposition.o mo_exception.o mo_mpi.o
mo_grib.o: mo_constants.o mo_control.o mo_decomposition.o mo_doctor.o \
	mo_exception.o mo_filename.o mo_gaussgrid.o mo_memory_g3b.o \
	mo_memory_gl.o mo_memory_sp.o mo_mpi.o mo_post.o mo_start_dataset.o \
	mo_time_control.o mo_tracer.o mo_transpose.o mo_year.o
mo_hdiff.o: mo_parameters.o
mo_hyb.o: mo_constants.o mo_control.o mo_parameters.o
mo_io.o: mo_control.o mo_decomposition.o mo_doctor.o mo_exception.o \
	mo_filename.o mo_io_tables.o mo_linked_list.o mo_memory_f.o \
	mo_memory_g1a.o mo_memory_g2a.o mo_memory_g3a.o mo_memory_gl.o \
	mo_mpi.o mo_netCDF.o mo_parameters.o mo_start_dataset.o mo_tracer.o \
	mo_transpose.o mo_year.o
mo_io_tables.o: mo_parameters.o
mo_landsea.o: mo_control.o mo_doctor.o
mo_legendre.o: mo_constants.o mo_control.o mo_decomposition.o mo_gaussgrid.o
mo_linked_list.o: mo_doctor.o mo_kind.o
mo_machine.o: mo_doctor.o mo_exception.o mo_kind.o mo_mpi.o
mo_memory_base.o: mo_doctor.o mo_kind.o mo_linked_list.o mo_netCDF.o
mo_memory_f.o: mo_kind.o mo_linked_list.o mo_memory_base.o mo_netCDF.o
mo_memory_g1a.o: mo_kind.o mo_linked_list.o mo_memory_base.o mo_netCDF.o
mo_memory_g1b.o: mo_kind.o mo_linked_list.o mo_memory_base.o
mo_memory_g2a.o: mo_kind.o mo_linked_list.o mo_memory_base.o mo_netCDF.o
mo_memory_g2b.o: mo_kind.o mo_linked_list.o mo_memory_base.o
mo_memory_g3a.o: mo_control.o mo_io_tables.o mo_kind.o mo_linked_list.o \
	mo_memory_base.o mo_netCDF.o
mo_memory_g3b.o: mo_control.o mo_io_tables.o mo_kind.o mo_linked_list.o \
	mo_memory_base.o mo_memory_g3a.o mo_netCDF.o
mo_memory_gl.o: mo_kind.o mo_linked_list.o mo_memory_base.o mo_netCDF.o
mo_memory_ls.o: mo_kind.o mo_linked_list.o mo_memory_base.o
mo_memory_sp.o: mo_kind.o mo_linked_list.o mo_memory_base.o
mo_midatm.o: mo_constants.o mo_control.o mo_decomposition.o mo_doctor.o \
	mo_exception.o mo_gaussgrid.o mo_io.o mo_memory_g3a.o mo_memory_g3b.o \
	mo_memory_sp.o mo_mpi.o mo_netCDF.o mo_param_switches.o \
	mo_start_dataset.o
mo_mpi.o: mo_doctor.o mo_kind.o
mo_netCDF.o: mo_control.o mo_doctor.o mo_exception.o mo_io_tables.o \
	mo_start_dataset.o mo_tracer.o
mo_nmi.o: mo_constants.o mo_control.o mo_doctor.o mo_exception.o mo_hyb.o \
	mo_memory_sp.o mo_nudging_buffer.o mo_rad_switches.o mo_semi_impl.o \
	mo_start_dataset.o mo_truncation.o $(INCLUDE)/nmictl.inc
mo_nudging.o: mo_constants.o mo_control.o mo_decomposition.o \
	mo_diag_tendency.o mo_doctor.o mo_exception.o mo_grib.o mo_kind.o \
	mo_memory_g3a.o mo_memory_g3b.o mo_memory_sp.o mo_mpi.o mo_nmi.o \
	mo_nudging_buffer.o mo_physc2.o mo_post.o mo_start_dataset.o \
	mo_time_control.o mo_transpose.o mo_truncation.o mo_year.o \
	$(INCLUDE)/ndgctl.inc
mo_parallel_slt.o: mo_decomposition.o mo_grid.o mo_mpi.o mo_transpose.o
mo_physc2.o: mo_parameters.o
mo_post.o: mo_kind.o mo_parameters.o
mo_rad2.o: mo_parameters.o
mo_scan_buffer.o: mo_decomposition.o mo_sc1.o mo_tracer.o
mo_semi_impl.o: mo_parameters.o
mo_skintem.o: mo_constants.o mo_control.o mo_param_switches.o mo_physc2.o \
	mo_radint.o
mo_soil_impl.o: mo_constants.o mo_control.o mo_diagnostics.o \
	mo_forecast_switches.o mo_gaussgrid.o mo_hdiff.o mo_hyb.o \
	mo_param_switches.o mo_parameters.o mo_physc1.o mo_physc2.o \
	mo_rad_switches.o mo_semi_impl.o mo_start_dataset.o mo_truncation.o
mo_sst.o: mo_constants.o mo_control.o mo_couple.o mo_decomposition.o \
	mo_doctor.o mo_dsst.o mo_exception.o mo_filename.o mo_io.o \
	mo_memory_g3a.o mo_memory_g3b.o mo_mpi.o mo_physc2.o \
	mo_rad_switches.o mo_start_dataset.o mo_timeint.o mo_transpose.o
mo_start_dataset.o: mo_doctor.o mo_mpi.o $(INCLUDE)/sdsctl.inc
mo_stat_global.o: mo_parameters.o
mo_stat_zonal.o: mo_parameters.o
mo_stratiform.o: mo_control.o mo_decomposition.o mo_doctor.o mo_exception.o \
	mo_kind.o mo_memory_g3a.o mo_mpi.o mo_transpose.o
mo_test_trans.o: mo_decomposition.o mo_doctor.o mo_exception.o mo_mpi.o \
	mo_transpose.o
mo_time_control.o: mo_constants.o mo_control.o mo_start_dataset.o mo_year.o
mo_timeint.o: mo_constants.o mo_control.o mo_parameters.o mo_start_dataset.o \
	mo_year.o
mo_tracer.o: mo_constants.o mo_control.o mo_doctor.o mo_exception.o \
	mo_gaussgrid.o mo_mpi.o mo_parameters.o mo_start_dataset.o \
	$(INCLUDE)/tractl.inc
mo_transpose.o: mo_decomposition.o mo_doctor.o mo_exception.o mo_mpi.o
mo_truncation.o: mo_doctor.o mo_mpi.o mo_parameters.o
mo_year.o: mo_julian.o mo_start_dataset.o

ifeq ($(ARCH), CRAY_PVP)
mo_buffer_fft.o:
	$(F90) $(F90FLAGS) -Ovector1 -c mo_buffer_fft.f90
mo_grib.o:
	$(F90) $(F90FLAGS) -Ovector1 -c mo_grib.f90
endif
