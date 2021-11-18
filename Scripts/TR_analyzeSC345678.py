#Purpose: Processing all Screen files with CasTLE analysis.
#Path must include bowtie files, set to: PATH=$PATH/usr/local/bin:/usr/local/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/Users/tysonruetz/Desktop/Dropbox/TR_Cas9Paper_CodeCheck/bowtie
#Directory: cd Desktop/Dropbox/TR_Cas9Paper_CodeCheck
#Terminal command: python2 Scripts/TR_analyzeSC345678.py

######################################
#Legend
# SC3 = Old screen1
# SC4 = Young screen 1
# SC5 = Young screen 2
# SC6 = Old screen 2
# SC7 = Young screen 3
# SC8 = Old screen 3
# Kp = Day 4 Ki67 positive (Kp) isolated cells
# RA = Day 14 reactivated (RA) cells
######################################

import subprocess

# #Make count files from fastq files (fq.gz)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC3_lib_FKDL190724827-1a-41_HWHG3CCXY_L2_1.fq.gz SC3_lib mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC3_O_Kp_FKDL190724827-1a-21_HWHG3CCXY_L2_1.fq.gz SC3_O_Kp mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC3_O_RA_FKDL190724827-1a-13_HWHG3CCXY_L2_1.fq.gz SC3_O_RA mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC4_Y_Kp_FKDL190724827-1a-38_HWHG3CCXY_L2_1.fq.gz SC4_Y_Kp mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC4_Y_RA_FKDL190724827-1a-37_HWHG3CCXY_L2_1.fq.gz SC4_Y_RA mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC5_Y_Kp_FKDL190724827-1a-6_HWHG3CCXY_L2_1.fq.gz SC5_Y_Kp mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC5_Y_RA_FKDL190724827-1a-5_HWHG3CCXY_L2_1.fq.gz SC5_Y_RA mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC6_O_Kp_FKDL190724827-1a-12_HWHG3CCXY_L2_1.fq.gz SC6_O_Kp mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC6_O_RA_FKDL190724827-1a-11_HWHG3CCXY_L2_1.fq.gz SC6_O_RA mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC7_Y_Kp_USPD16100718-16_H2V33CCX2_L5_1.fq.gz SC7_Y_Kp mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC7_Y_RA_USPD16100718-15_H2V33CCX2_L5_1.fq.gz SC7_Y_RA mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC8_O_Kp_USPD16100718-20_H2V33CCX2_L5_1.fq.gz SC8_O_Kp mm-Cas9-10,', shell=True)
# subprocess.call('python2 Scripts/makeCounts.py fastq_files/SC8_O_RA_USPD16100718-19_H2V33CCX2_L5_1.fq.gz SC8_O_RA mm-Cas9-10,', shell=True)


# #CasTLE Analyze Counts of each screen compared to starting library
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC3_O_Kp_counts.csv SC3_KpvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC3_O_RA_counts.csv SC3_RAvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC4_Y_Kp_counts.csv SC4_KpvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC4_Y_RA_counts.csv SC4_RAvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC5_Y_Kp_counts.csv SC5_KpvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC5_Y_RA_counts.csv SC5_RAvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC6_O_Kp_counts.csv SC6_KpvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC6_O_RA_counts.csv SC6_RAvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC7_Y_Kp_counts.csv SC7_KpvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC7_Y_RA_counts.csv SC7_RAvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC8_O_Kp_counts.csv SC8_KpvLib', shell=True)
# subprocess.call('python2 Scripts/analyzeCounts.py Data/SC3_lib_counts.csv Data/SC8_O_RA_counts.csv SC8_RAvLib', shell=True)


# # Add P-values with 100,000 permutations)
# subprocess.call('python Scripts/addPermutations.py Results/SC3_KpvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC3_RAvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC4_KpvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC4_RAvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC5_KpvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC5_RAvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC6_KpvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC6_RAvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC7_KpvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC7_RAvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC8_KpvLib.csv 100000', shell=True)
# subprocess.call('python Scripts/addPermutations.py Results/SC8_RAvLib.csv 100000', shell=True)
