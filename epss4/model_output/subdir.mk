################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
F90_SRCS += \
../asseteuler.f90 \
../distribution.f90 \
../income.f90 \
../kinds.f90 \
../krusell_smith.f90 \
../laws_of_motion.f90 \
../main.f90 \
../meanshock_equilib.f90 \
../meanshock_wrapper.f90 \
../params_mod.f90 \
../policyfunctions.f90 \
../run_model_mod.f90 \
../save_results_mod.f90 \
../simulate_economy.f90 \
../statistics.f90 \
../types.f90 

OBJS += \
./asseteuler.o \
./distribution.o \
./income.o \
./kinds.o \
./krusell_smith.o \
./laws_of_motion.o \
./main.o \
./meanshock_equilib.o \
./meanshock_wrapper.o \
./params_mod.o \
./policyfunctions.o \
./run_model_mod.o \
./save_results_mod.o \
./simulate_economy.o \
./statistics.o \
./types.o 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.f90
	@echo 'Building file: $<'
	@echo 'Invoking: Intel(R) Intel(R) 64 Fortran Compiler'
	ifort -O3 -I"${MKLROOT}/include" -I"${MKLROOT}/include/intel64/lp64" -check none -c -assume realloc_lhs -no-prec-div -xHost -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

asseteuler.o: ../asseteuler.f90 kinds.o params_mod.o src_utilities/fun_lininterp.o

distribution.o: ../distribution.f90 kinds.o params_mod.o src_utilities/fun_locate.o

income.o: ../income.f90 kinds.o params_mod.o

kinds.o: ../kinds.f90

krusell_smith.o: ../krusell_smith.f90 kinds.o laws_of_motion.o params_mod.o policyfunctions.o save_results_mod.o simulate_economy.o src_utilities/aggregate_grids.o src_utilities/error_class.o src_utilities/interpolate_xgrid.o types.o

laws_of_motion.o: ../laws_of_motion.f90 kinds.o params_mod.o src_utilities/aggregate_grids.o types.o

main.o: ../main.f90 params_mod.o run_model_mod.o

meanshock_equilib.o: ../meanshock_equilib.f90 distribution.o income.o kinds.o laws_of_motion.o params_mod.o policyfunctions.o src_utilities/aggregate_grids.o src_utilities/error_class.o types.o

meanshock_wrapper.o: ../meanshock_wrapper.f90 distribution.o income.o kinds.o laws_of_motion.o meanshock_equilib.o params_mod.o src_utilities/aggregate_grids.o src_utilities/error_class.o src_utilities/fun_aggregate_diff.o src_utilities/fun_locate.o src_utilities/partial_sorting.o src_utilities/sub_broyden.o types.o

params_mod.o: ../params_mod.f90 kinds.o src_utilities/aggregate_grids.o src_utilities/fun_kronprod.o src_utilities/makegrid_mod.o src_utilities/markov_chain_approx.o src_utilities/markov_station_distr.o src_utilities/twostate_exact.o

policyfunctions.o: ../policyfunctions.f90 asseteuler.o income.o kinds.o laws_of_motion.o params_mod.o src_utilities/aggregate_grids.o src_utilities/error_class.o src_utilities/fun_lininterp.o src_utilities/fun_locate.o src_utilities/fun_zbrent.o src_utilities/makegrid_mod.o src_utilities/sub_zbrac.o src_utilities/sub_zbrak.o types.o

run_model_mod.o: ../run_model_mod.f90 distribution.o income.o krusell_smith.o laws_of_motion.o meanshock_wrapper.o params_mod.o save_results_mod.o src_utilities/aggregate_grids.o src_utilities/error_class.o src_utilities/markovchain_mod.o src_utilities/numrec_utils.o src_utilities/sub_alg_qn.o src_utilities/unformatted_io.o statistics.o types.o

save_results_mod.o: ../save_results_mod.f90 income.o kinds.o laws_of_motion.o params_mod.o src_utilities/aggregate_grids.o src_utilities/error_class.o statistics.o types.o

simulate_economy.o: ../simulate_economy.f90 distribution.o income.o kinds.o params_mod.o src_utilities/aggregate_grids.o src_utilities/fun_aggregate_diff.o src_utilities/fun_excessbonds.o src_utilities/fun_locate.o src_utilities/fun_zbrent.o src_utilities/partial_sorting.o types.o

statistics.o: ../statistics.f90 kinds.o params_mod.o

types.o: ../types.f90 kinds.o params_mod.o


