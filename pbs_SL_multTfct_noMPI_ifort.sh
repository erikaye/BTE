#!/bin/sh

# PBS script example for parallel jobs using Open MPI/GNU 4.1
# Jen Cobb, Matthew Inkman 09/25/11

# sets the number of nodes & processors per node (ppn)
# can also explicitly specify which nodes to use: nodes001:ppn=6+n010 [...]
#PBS -l nodes=1:ppn=4,walltime=15:00:00

# name of your job
#PBS -N pbs_test1

# Combine output and error streams into single file
#PBS -j oe

# # Name of queue to submit job to
#PBS -q reg_20

# sets the maximum amount of time the job can run for (hr:min:sec)
# PBS -l walltime=8:00:00

numproc=$(cat $PBS_NODEFILE | wc -l)

cd $PBS_O_WORKDIR
echo The number of processes is...$numproc
echo Working directory is...$(pwd)
ulimit -s

fortpath='/opt/intel/composer_xe_2011_sp1.8.273'

source $fortpath/bin/ifortvars.sh intel64
source $fortpath/bin/compilervars.sh intel64
echo $MKLROOT
# export KMP_AFFINITY=verbose,scatter

fileVar='1'
outName="output_np${numproc}_${fileVar}" 
echo $outName

# approx_old_v2 for same material;  AMM_v8 for diff materials
$fortpath/bin/intel64/ifort -O3 -r8 -c Constants_dp.f90 orderParams.f90 GetMatlParams.f90 GetCosExpCoeffs.f90 \
GetRTMats_approx_old_v3.f90 GetInterfaceBM_vSparse_v1.f90 GetInterfaceBM_PBC.f90 CheckBCs_BBbc.f90 CheckBCs_Pbc.f90 \
MultilayerT_PBc_fct_constL.f90 MultilayerT_PBc_fct.f90 MultilayerT_BBc_fct_noMPI.f90 MultilayerT_BBc_fct_noQ.f90 \
runCalcs_multTfct.f90 -i8 -openmp -I${MKLROOT}/include/intel64/ilp64 -I${MKLROOT}/include \
-traceback -check bounds -heap-arrays
# -ftz sets underflow to 0.  done when -O3 is used.
# default -mkl=parallel (tells compiler to link using threaded libraries in mkl)

$fortpath/bin/intel64/ifort -O3 -r8 Constants_dp.o orderParams.o GetMatlParams.o GetCosExpCoeffs.o GetRTMats_approx_old_v3.o \
GetInterfaceBM_vSparse_v1.o GetInterfaceBM_PBC.o CheckBCs_BBbc.o CheckBCs_Pbc.o MultilayerT_PBc_fct_constL.o MultilayerT_PBc_fct.o \
MultilayerT_BBc_fct_noMPI.o MultilayerT_BBc_fct_noQ.o runCalcs_multTfct.o -o $outName ${MKLROOT}/lib/intel64/libmkl_blas95_ilp64.a \
${MKLROOT}/lib/intel64/libmkl_lapack95_ilp64.a \
-L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_intel_thread \
-lmkl_core -liomp5 -lpthread -lm



# modify to include the correct mpirun and your executable
# time ./$outName SiMgO_inputs_hooke_1 > SiSiX_outputs_hooke_sp1_1.txt
# time ./$outName SiMgO_inputs_hooke_2 > SiSiX_outputs_hooke_sp1_2.txt
# time ./$outName SiMgO_inputs_hooke_3 > SiSiX_outputs_hooke_sp1_3.txt
# time ./$outName SiMgO_inputs_hooke_4 > SiSiX_outputs_hooke_sp1_4.txt
# time ./$outName SiMgO_inputs_hooke_5 > SiMgO_outputs_hooke_p5_N100.txt
# time ./$outName SiMgO_inputs_hooke_6 > SiMgO_outputs_hooke_p6_N100.txt
# time ./$outName SiMgO_inputs_hooke_7 > SiMgO_outputs_hooke_p7_N100.txt
# time ./$outName SiMgO_inputs_hooke_8 > SiMgO_outputs_hooke_p8_N100.txt

# time ./$outName SiSi_inputs_hooke_l > SiSi_outputs_hooke_spZ_l.txt
# time ./$outName SiSi_inputs_hooke_m > SiSi_outputs_hooke_spZ_m.txt

# time ./$outName SiSi_inputs_hooke_1 > SiSi_outputs_hooke_spZ_1.txt
# time ./$outName SiSi_inputs_hooke_2 > SiSi_outputs_hooke_spZ_2.txt
time ./$outName SiSi_inputs_hooke_3 > SiSi_outputs_hooke_spZ_3.txt
# time ./$outName SiSi_inputs_hooke_4 > SiSi_outputs_hooke_spZ_4.txt
# time ./$outName SiSi_inputs_hooke_5 > SiSi_outputs_hooke_spZ_5.txt
# time ./$outName SiSi_inputs_hooke_6 > SiSi_outputs_hooke_spZ_6.txt
# time ./$outName SiSi_inputs_hooke_7 > SiSi_outputs_hooke_spZ_7.txt
# time ./$outName SiSi_inputs_hooke_8 > SiSi_outputs_hooke_spZ_8.txt

# time ./$outName SiSi_inputs_hooke_a > SiSi_outputs_hooke_spZ_a.txt
# time ./$outName SiSi_inputs_hooke_b > SiSi_outputs_hooke_spZ_b.txt
# time ./$outName SiSi_inputs_hooke_c > SiSi_outputs_hooke_spZ_c.txt
# time ./$outName SiSi_inputs_hooke_d > SiSi_outputs_hooke_spZ_d.txt
# time ./$outName SiSi_inputs_hooke_e > SiSi_outputs_hooke_spZ_e.txt
# time ./$outName SiSi_inputs_hooke_f > SiSi_outputs_hooke_spZ_f.txt
# time ./$outName SiSi_inputs_hooke_g > SiSi_outputs_hooke_spZ_g.txt
# time ./$outName SiSi_inputs_hooke_h > SiSi_outputs_hooke_spZ_h.txt


# /opt/intel/impi/4.0.3.008/intel64/bin/mpirun -np 1 time ./output
