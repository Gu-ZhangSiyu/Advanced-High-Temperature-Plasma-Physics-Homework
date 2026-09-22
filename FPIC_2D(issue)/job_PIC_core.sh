#!/bin/bash
# build_and_submit.sh

set -euo pipefail

runfile="PIC_CORE_Zhang"

module load intel/2025.1

mpiifx -O3 -heap-arrays \
  pic_modules.f90 \
  boris_solver_yee.f90 \
  Domain_Decomposition_Module.f90 \
  Ghost_Cells_Module.f90 \
  diagnostics_solver.f90 \
  data_loader_module.f90 \
  restart_manager_module.f90 \
  Particle_Migration_Module.f90 \
  beam_injection_module.f90 \
  Utility_Functions.f90 \
  Particle_Operations.f90 \
  Load_Balance_Module.f90 \
  pic_core.f90 \
  -o "${runfile}"

mkdir -p dynamic_electrons_position
mkdir -p dynamic_ions_position
mkdir -p dynamic_electrons_velocity
mkdir -p dynamic_ions_velocity
mkdir -p dynamic_fields_Bx
mkdir -p dynamic_fields_By
mkdir -p dynamic_fields_Bz
mkdir -p dynamic_fields_Ex
mkdir -p dynamic_fields_Ey
mkdir -p dynamic_fields_Ez
mkdir -p diagnostics
JOB_SCRIPT="${runfile}_qsub_job.sh"
cat > "${JOB_SCRIPT}" << EOF
#!/bin/bash
#------- qsub option --------
#PBS -P NIFS26KISC037
#PBS -q A_S
#PBS -N ${runfile}
#PBS -o ${runfile}.log
#PBS -e ${runfile}.err
#PBS -l select=1:ncpus=128:mpiprocs=128:ompthreads=1
#PBS -l walltime=24:00:00
#PBS -m be
#PBS -M t262d003@gunma-u.ac.jp
#-------Program execution --------
module load intel/2025.1
cd \${PBS_O_WORKDIR}
date

export I_MPI_DEBUG=5
ulimit -s unlimited
ulimit -v unlimited
ulimit -l unlimited

mpirun ./${runfile}
date
EOF

qsub "${JOB_SCRIPT}"
