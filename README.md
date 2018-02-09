# CHOP

FUNCTION: The pipeline creates bash scripts of all samples in a tab separated file. These bash scripts includes process of Alignment, Sorting, Quality Metrics, Addreplacegrps, Haplotype Caller.

Haplotype caller was used to call variants because the test sample didnot have tumor_matched_normal. In between the process it generated baits file for hsmetrics and interval list required for quality metrics and calling variants. The scripts folder generated should have all the scripts depending upon number of samples.

Once all scripts are generated bash scipt.sh cmd should be executed depecding upon avaialbility of resources, This pipeline is automated and can have a control on number of samples to be processed at a time.

* For SLURM change the header "#/bin/bash" of bash scripts to the slurm options and save it as ".slurm"

INPUT REQUIRED 
1. Tab separated file for fastq paired files
             
2. Tab separated file describing locations to respected variables (tool kit and other files) decribed in a template file.

Requirements :
Python/3.5.2 ; 
Java (version-1.8.x) ; 
BWA (version-0.7.x) ; 
GATK (version-3.6 or higher) ; 
Picard Tools (version-2.x.x)

Execution:
    
$ python/3.5.2 pipeline.py sample_file location_file

#For Bash:
$ bash script.sh 

#For SLURM:
$ sbatch script.slurm
