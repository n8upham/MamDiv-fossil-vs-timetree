#BSUB -J BAMM_100Mgens_1mamPhy
#BSUB -n 10
#BSUB -R "span[]"  
#BSUB -q week
#BSUB -W 168:00
#BSUB -o outputs/BAMM_10Mgens_1mamPhy
#BSUB -N

echo "Slot count: $LSB_DJOB_NUMPROC"
/apps/hpc/Apps/BAMM/2.5.0/bin/bamm -c mamPhy_controlFile_2.txt
