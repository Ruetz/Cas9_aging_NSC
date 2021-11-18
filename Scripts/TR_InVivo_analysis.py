## Analyzing the in vivo screen data
#Path must include bowtie files, set to: PATH=$PATH/usr/local/bin:/usr/local/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/Users/tysonruetz/Desktop/Dropbox/TR_Cas9Paper_CodeCheck_BACKUP/bowtie
#Directory: cd Desktop/Dropbox/TR_Cas9Paper_CodeCheck



import subprocess

subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/TOC1dOB_S5_L001_R1_001.fastq.gz T10_24hr_OB Top10_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/TOC1Cr_S6_L001_R1_001.fastq.gz T10_24hr_Cer Top10_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/TOC1dCt_S7_L001_R1_001.fastq.gz T10_24hr_Cort Top10_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/T10_SF_5wk_SVZ.fastq.gz T10_5wk_SVZ Top10_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/TOC3wCer_S16_L001_R1_001.fastq.gz T10_5wk_Cer Top10_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/TOC3wCor_S15_L001_R1_001.fastq.gz T10_5wk_Cort Top10_safe100', shell=True)



subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/T10S-24_S13_L001_R1_001.fastq.gz IV9_T10S_24hr Top10_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/T10S-O1_S1_L001_R1_001.fastq.gz IV9_T10S_O1 Top10_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/T10S-O2_S2_L001_R1_001.fastq.gz IV9_T10S_O2 Top10_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/DSS-24_S15_L001_R1_001.fastq.gz IV9_DISS_24hr Dis_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/DSS-O1_S7_L001_R1_001.fastq.gz IV9_DISS_O1 Dis_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/3-CO-DIS_S12_L001_R1_001.fastq.gz Invivo_IV3_CO_DIS Dis_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/GRS-24_S17_L001_R1_001.fastq.gz IV9_GRS_24hr Gran_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/GR-S-Lib_S34_L001_R1_001.fastq.gz I6_GR_SV1 Gran_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/GRS-O2_S10_L001_R1_001.fastq.gz IV9_GRS_O2 Gran_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/GRS-O1_S9_L001_R1_001.fastq.gz IV9_GRS_O1 Gran_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/DPS-24_S16_L001_R1_001.fastq.gz IV9_DEP_24hr Dep_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/DPOCOB2_S21_L001_R1_001.fastq.gz IV6_DPS_O_OB2 Dep_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/DPS-O2_S12_L001_R1_001.fastq.gz IV9_DEP_O2 Dep_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/CGS-24_S14_L001_R1_001.fastq.gz IV9_CGP_24hr CGP1_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/CGS-O1_S5_L001_R1_001.fastq.gz IV9_CGP_O1 CGP1_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/CGS-O2_S6_L001_R1_001.fastq.gz IV9_CGP_O2 CGP1_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/CGOCSV1_S14_L001_R1_001.fastq.gz IV6_CGPS_O_SV CGP1_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/TOWTOB1_S1_L001_R1_001.fastq.gz WT_O_OB1 Top10_safe100', shell=True)
subprocess.call('python2 Scripts/makeCounts.py fastq_files/InVivo/TOWTOB2_S3_L001_R1_001.fastq.gz WT_O_OB2 Top10_safe100', shell=True)


subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_T10S_24hr_counts.csv Data/IV9_T10S_O1_counts.csv Top10_Sc1', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_T10S_24hr_counts.csv Data/IV9_T10S_O2_counts.csv Top10_Sc2', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_DISS_24hr_counts.csv Data/IV9_DISS_O1_counts.csv Glucose_HumanDisease_Sc1', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_DISS_24hr_counts.csv Data/Invivo_IV3_CO_DIS_counts.csv Glucose_Humandisease_Sc2', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_GRS_24hr_counts.csv Data/IV9_GRS_O1_counts.csv CytoplasmicRibonucleoproteinStructures_Sc1', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/I6_GR_SV1_counts.csv Data/IV9_GRS_O2_counts.csv CytoplasmicRibonucleoproteinStructures_Sc2', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_DEP_24hr_counts.csv Data/IV6_DPS_O_OB2_counts.csv Depleted_Sc1', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_DEP_24hr_counts.csv Data/IV9_DEP_O2_counts.csv Depleted_Sc2', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV6_CGPS_O_SV_counts.csv Data/IV9_CGP_O1_counts.csv PublishedRegulators_Sc1', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_CGP_24hr_counts.csv Data/IV9_CGP_O2_counts.csv PublishedRegulators_Sc2', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_T10S_24hr_counts.csv Data/WT_O_OB1_counts.csv WT_Top10_Sc1', shell=True)
subprocess.call('python2 Scripts/analyzeCounts.py Data/IV9_T10S_24hr_counts.csv Data/WT_O_OB2_counts.csv WT_Top10_Sc2', shell=True)


#
# subprocess.call('python2 Scripts/addPermutations.py Results/Top10_Sc1.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/Top10_Sc2.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/Glucose_HumanDisease_Sc1.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/Glucose_Humandisease_Sc2.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/CytoplasmicRibonucleoproteinStructures_Sc1.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/CytoplasmicRibonucleoproteinStructures_Sc2.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/Depleted_Sc1.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/Depleted_Sc2.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/PublishedRegulators_Sc1.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/PublishedRegulators_Sc2.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/WT_Top10_Sc1.csv 10000', shell=True)
# subprocess.call('python2 Scripts/addPermutations.py Results/WT_Top10_Sc2.csv 10000', shell=True)
