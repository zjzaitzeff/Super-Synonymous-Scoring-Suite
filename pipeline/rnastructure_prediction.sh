export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

inputPath=$1
outputPath=$2

./RNAstructure/exe/Fold-smp $inputPath $outputPath