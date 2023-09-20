#!/bin/bash
#SBATCH --job-name=Thao
#SBATCH -p E5,Lake,Cascade
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2         # number of MPI processes per node
#SBATCH --cpus-per-task=8           # number of OpenMP threads per MPI process
#SBATCH --time=0-10:30:00
#
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=thi.nguyen1@ens-lyon.fr
module purge
module use /applis/PSMN/debian11/E5/modules/all
module load GCC/10.3.0
module load GSL
cd simuscale_APC
rm -rf build
mkdir build && cd build
cmake -DCMAKE_INSTALL_PREFIX=/home/tnguye22/bin ..
make simuscale
cd bin/
rm -rf ../../run/lymphocytes/Day_1_3decoratings_22.txt 
\cp ../../run/lymphocytes/param_3decoratings_9G_Day1.in ../../run/lymphocytes/param_3decoratings_9G.in
./simuscale
awk '$1 == 24.00 {print}' ../../run/lymphocytes/trajectory.txt > ../../run/lymphocytes/Day_1.txt
cut --complement -d ' ' -f2,4-14 ../../run/lymphocytes/Day_1.txt > ../../run/lymphocytes/Day_1_3decoratings_22.txt 
\cp ../../run/lymphocytes/param_3decoratings_9G_T.in ../../run/lymphocytes/param_3decoratings_9G.in 
./simuscale
