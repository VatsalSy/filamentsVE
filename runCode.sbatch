!/bin/bash

#SBATCH --partition=pm6-isw2,pm9-isw0,pm11-isw2
#SBATCH --job-name=die-swell_1
#SBATCH --account=ehpc-reg-2023r03-178
#SBATCH --qos=ehpc-reg-2023r03-178
#SBATCH --time=12:00:00

#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --ntasks-per-core=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=251G
#SBATCH -e job.%J.err
#SBATCH -o job.%J.out

source ~/.bash_shell


export OMP_NUM_THREADS=64
./filament_initialCondition