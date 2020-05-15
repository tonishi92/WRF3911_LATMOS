#!/bin/bash

# To start this job: qsub run_real.pbs

# Resources used
#PBS -l nodes=1:ppn=1
#PBS -q day 
#PBS -j eo
#PBS -k eo
#PBS -l mem=40G,vmem=40G

# Send a mail when the job ends
#PBS -m ae
#PBS -M tatsuo.onishi@latmos.ipsl.fr
cd $PBS_O_WORKDIR

# Link the necessary files from the SI directory

cat > prepare_for_run_ifort << EOF
module purge
module load intel/15.0.6.233 
module load openmpi/1.6.5-ifort
module load netcdf4/4.2.1.1-ifort
module load hdf5/1.8.14-ifort
export CFLAGS="-I/usr/lib64/openmpi/1.6.5-ifort/include -m64"
export LDFLAGS="-L/usr/lib64/openmpi/1.6.5-ifort/lib -lmpi"
export NETCDF=/opt/netcdf42/ifort
export PHDF5=/opt/hdf5/1.8.14/ifort
export HDF5=/opt/hdf5/1.8.14/ifort
EOF
source prepare_for_run_ifort
export MPI_LIB=-L$MPI_LIB
export TOPDIR=$PWD/..
export WRFDIR=$PWD
export TOOLDIR=$HOME/tools_ng
export SRCDIR=/data/$USER/sources_ng
export PATH=$TOOLDIR/bin:$PATH
export PATH=$TOOLDIR/lib:$PATH
export WRFIO_NCD_NO_LARGE_FILE_SUPPORT=0
export JASPERLIB=/usr/lib64
export JASPERINC=/usr/include/jasper
export WRF_EM_CORE=1
export WRF_NMM_CORE=0
export WRF_CHEM=1
export WRF_KPP=1
export YACC="$TOOLDIR/bin/yacc -d"
export FLEX_LIB_DIR=$TOOLDIR/lib
export HDF5_DISABLE_VERSION_CHECK=1


cd $WRFDIR
./clean -aa

cd $WRFDIR/chem/KPP/kpp/kpp-2.1/src/
flex scan.l

sed -i '
1 i \
#define INITIAL 0 \
#define CMD_STATE 1 \
#define INC_STATE 2 \
#define MOD_STATE 3 \
#define INT_STATE 4 \
#define PRM_STATE 5 \
#define DSP_STATE 6 \
#define SSP_STATE 7 \
#define INI_STATE 8 \
#define EQN_STATE 9 \
#define EQNTAG_STATE 10 \
#define RATE_STATE 11 \
#define LMP_STATE 12 \
#define CR_IGNORE 13 \
#define SC_IGNORE 14 \
#define ATM_STATE 15 \
#define LKT_STATE 16 \
#define INL_STATE 17 \
#define MNI_STATE 18 \
#define TPT_STATE 19 \
#define USE_STATE 20 \
#define COMMENT 21 \
#define COMMENT2 22 \
#define EQN_ID 23 \
#define INL_CODE 24
' $WRFDIR/chem/KPP/kpp/kpp-2.1/src/lex.yy.c

cd $WRFDIR
./configure <<EOF
15

EOF

./compile em_real 2>&1 |tee compile.log

