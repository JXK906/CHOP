#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  4 17:42:37 2018

@author: jxk906
"""
import sys
import os
import pandas as pd
import numpy as np
import subprocess as sub


''' Variant Calling Pipeline
INPUT REQUIRED : 1. Tab separated file for fastq paired files
                 2. Tab separated file describing locations to respected variables (tool kit and other files) decribed in a template file

FUNCTION: The pipeline creates bash scripts of all samples in a tab separated file. These bash scripts includes process of Alignment, Sorting, Quality Metrics, Addreplacegrps, HAplotype Caller.
haplotype caller was used as the test sample was a single caller. In between the process it generated baits file for hsmetrics and interval list required to call variants. The scripts folder 
generated should have all the scripts. 

Once all scripts are generated bash scipt.sh cmd should be executed depecding upon avaialbility of resources, This pipeline is automated and can have a conrol on number of samples to be processed at a time.
* CAN USE SLURM AS WELL
Requirements :
    Python/3.5.2
    Java (version-1.8.x)
    BWA (version-0.7.x)
    GATK.jar (version-3.6 or higher)
    Picard.jar file (version-2.x.x)

execution:
    
$ python/3.5.2 pipeline.py sample_file location_file
$ bash script.sh

OUTPUT asfter execution of bash script: 1. SAM file  2. sorted BAM file  3. Quality Metrics file  4.BAM with read groups replaced  5. VCF file  

'''
## CHOP PIPELINE

def pipeline(sample_file,location_file):
    sample_df = pd.read_table(sample_file,sep="\t", index_col=False,dtype=None) # reads sample file and generate list of all paired end reads
    sample_list = sample_df.values.tolist()
    locate = pd.read_table(location_file,sep="\t", index_col=False,dtype=None) # reads locations and provides dictionary of all variables
    locate_T = locate.transpose()
    locate_T.rename(columns=locate_T.iloc[0], inplace=True)
    locate_T.drop(['file_type'],inplace=True)
    pod = locate_T.to_dict("records")[0] ## Dictionary
    ## Generate bait interval file for HSmetrics (Required file)
    reference_dict = pod['Reference_Genome'].split(".")[0]+".dict" # reads a reference dictionary 
    df_refdict = pd.read_table(reference_dict,sep=" ",dtype=None, index_col=False, header=None)
    fil = pd.read_table(pod['RegionOfInterest_BED'], sep="\t", index_col=False,header=None, dtype=None) # Reads  ROI.bed file 
    bait = fil[[0,1,2,5,3]] # columns ['CHR', 'START', 'END', 'STRAND','ACESSIONID']
    bait.to_csv("/".join([pod['Project_Location'],'bait_tmp.txt']), sep="\t", index=False)
    _bait_ = pd.read_table("/".join([pod['Project_Location'],'bait_tmp.txt']), sep=" ", index_col=False,dtype=None) # creates a temp bait file having selected columns
    df_merge = pd.concat([df_refdict, _bait_]) 
    df_merge.to_csv("/".join([pod['Project_Location'],'bait.intervals']),header=False, sep =" ",index=False) ## FINAL BAIT Interval file
    BAIT_File= "/".join([pod['Project_Location'],'bait.intervals'])
    sub.call("rm %s"%"/".join([pod['Project_Location'],'bait_tmp.txt']),shell=True) # removal of temporary bait file
    ## MAIN FUNCTION starts here, generates all bash  scripts for all samples and creates individual directories for each process at the project location
    for pair in sample_list:
        script_dir = "/".join([pod['Project_Location'],'Scripts']) # Script directory
        if not os.path.exists(script_dir):
            os.mkdir(script_dir)
        print(pair)
        R1 =pair[0] # paired read 1
        R2 = pair[1] # paired read 2
        samp = R1.split(".")[0]
        print(samp)
        with open("/".join([script_dir,samp+".sh"]),"w") as bash:
            bash.write("#!/bin/bash\n")
            ## ALIGNMENT
            align_dir = "/".join([pod['Project_Location'],'Alignment'])
            if not os.path.exists(align_dir):
                os.mkdir(align_dir)
            Alignment = "%s mem -t 2 %s %s %s > %s \nwait"%(pod['BWA_Path'], pod['Reference_Genome'],"/".join([pod['Fastq_location'],R1]),"/".join([pod['Fastq_location'],R2]), "/".join([align_dir,samp+".sam"]))
            bash.write(Alignment)
            bash.write("\nwait")
            ## SORTING
            sort_dir = "/".join([pod['Project_Location'],'Sorting'])
            if not os.path.exists(sort_dir):
                os.mkdir(sort_dir)
            Sorting = "%s -Xmx8g -jar %s SortSam INPUT=%s OUTPUT=%s SORT_ORDER=coordinate CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT \nwait"%(pod['JAVA_Location'], pod['Picard.jar'],"/".join([align_dir,samp+".sam"]), "/".join([sort_dir,samp+".bam"]))
            bash.write(Sorting)
            bash.write("\nwait")
            ## HSMETRICS
            hsmet_dir = "/".join([pod['Project_Location'],'HSMETRICS'])
            if not os.path.exists(hsmet_dir):
                os.mkdir(hsmet_dir)
            hsmet = "%s -Xmx8g -jar %s CalculateHsMetrics I=%s O=%s R=%s BAIT_INTERVALS=%s TARGET_INTERVALS=%s \nwait\n"%(pod['JAVA_Location'], pod['Picard.jar'],"/".join([sort_dir,samp+".bam"]),"/".join([hsmet_dir,samp+".txt"]),pod['Reference_Genome'],BAIT_File,BAIT_File)
            bash.write(hsmet)
            bash.write("\n")
            bash.write("tail -n +7 %s > %s\nwait"%("/".join([hsmet_dir,samp+".txt"]),"/".join([hsmet_dir,samp+"Metrics.txt"])))
            bash.write("\nwait")
            ## ADDREPLACE_READGROUPs
            addreplacegrps_dir = "/".join([pod['Project_Location'],'Addreplacegrps'])
            if not os.path.exists(addreplacegrps_dir):
                os.mkdir(addreplacegrps_dir)
            Addrep = "%s -Xmx8g -jar %s AddOrReplaceReadGroups INPUT=%s OUTPUT=%s RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=%s CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT \nwait"%(pod['JAVA_Location'], pod['Picard.jar'],"/".join([sort_dir,samp+".bam"]),"/".join([addreplacegrps_dir,samp+".bam"]),samp)
            bash.write(Addrep)
            bash.write("\nwait")
            ## SINGLE SAMPLE HAPLOTYPE CALLER
            Var_call_dir = "/".join([pod['Project_Location'],"GATK_Caller"])
            if not os.path.exists(Var_call_dir ):
                os.mkdir(Var_call_dir)

            field2 = fil[[7]] 
            uniq = np.unique(field2.values).tolist()
            with open("/".join([pod['Project_Location'],"haplotype.interval_list"]),"w") as bedf:
                for i in uniq:
                    bedf.write(i)
                    bedf.write("\n")
            haplotypecal = "%s -Xmx8g -jar %s -R %s -T HaplotypeCaller -I %s -L %s -o %s\n"%(pod['JAVA_Location'], pod['GATK_Location'],pod['Reference_Genome'],"/".join([addreplacegrps_dir,samp+".bam"]),"/".join([pod['Project_Location'],"haplotype.interval_list"]),"/".join([Var_call_dir,samp+".vcf"]))
            bash.write(haplotypecal)
            bash.write("wait")
            
           
#pipeline("/mnt/pan/Data16/jxk906/PROJECTS/CHOP/samples.txt","/mnt/pan/Data16/jxk906/PROJECTS/CHOP/location_file")

if __name__=="__main__":
    sample_file=sys.argv[1]
    location_file = sys.argv[2]
pipeline(sample_file,location_file)
