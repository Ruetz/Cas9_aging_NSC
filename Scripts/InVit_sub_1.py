#InVitro Subscreen 1 analysis
##Path alread set to: export PATH="/usr/local/bin:/usr/local/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/Users/tysonruetz/Desktop/Dropbox/castle/bowtie"
##cd Desktop/Dropbox/castle

python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/48hr-lib_S1_L001_R1_001.fastq.gz InVitro_sub1_48hr Top_50
python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/O2-RA_S10_L001_R1_001.fastq.gz InVitro_sub1_O2_RA Top_50
python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/O3-RA_S12_L001_R1_001.fastq.gz InVitro_sub1_O3_RA Top_50
python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/O4-RA_S14_L001_R1_001.fastq.gz InVitro_sub1_O4_RA Top_50
python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/Y4-RA_S3_L001_R1_001.fastq.gz InVitro_sub1_Y4_RA Top_50
python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/Y5-RA_S5_L001_R1_001.fastq.gz InVitro_sub1_Y5_RA Top_50
python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/Y6-RA_S7_L001_R1_001.fastq.gz InVitro_sub1_Y6_RA Top_50
python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/Y7-RA_S8_L001_R1_001.fastq.gz InVitro_sub1_Y7_RA Top_50

python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/48hr-lib_S1_L001_R1_001.fastq.gz InVitro_sub1_48hr_Top10 Top10_safe100
python2 Scripts/makeCounts.py InVitro_subsc1_fastQ/O2-RA_S10_L001_R1_001.fastq.gz InVitro_sub1_O2_RA_Top10 Top10_safe100
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_48hr_Top10_counts.csv Data/InVitro_sub1_O2_RA_Top10_counts.csv InVitro_sub1_O2_RAvLib_Top10


python2 Scripts/analyzeCounts.py Data/InVitro_sub1_48hr_counts.csv Data/InVitro_sub1_O2_RA_counts.csv InVitro_sub1_O2_RAvLib
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_48hr_counts.csv Data/InVitro_sub1_O3_RA_counts.csv InVitro_sub1_O3_RAvLib
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_48hr_counts.csv Data/InVitro_sub1_O4_RA_counts.csv InVitro_sub1_O4_RAvLib
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_48hr_counts.csv Data/InVitro_sub1_Y4_RA_counts.csv InVitro_sub1_Y4_RAvLib
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_48hr_counts.csv Data/InVitro_sub1_Y5_RA_counts.csv InVitro_sub1_Y5_RAvLib
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_48hr_counts.csv Data/InVitro_sub1_Y6_RA_counts.csv InVitro_sub1_Y6_RAvLib
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_48hr_counts.csv Data/InVitro_sub1_Y7_RA_counts.csv InVitro_sub1_Y7_RAvLib

python2 Scripts/analyzeCounts.py Data/InVitro_sub1_O2_RA_counts.csv Data/InVitro_sub1_Y4_RA_counts.csv InVitro_sub1_O2vY4
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_O3_RA_counts.csv Data/InVitro_sub1_Y5_RA_counts.csv InVitro_sub1_O3vY5
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_O4_RA_counts.csv Data/InVitro_sub1_Y6_RA_counts.csv InVitro_sub1_O4vY6
python2 Scripts/analyzeCounts.py Data/InVitro_sub1_O4_RA_counts.csv Data/InVitro_sub1_Y7_RA_counts.csv InVitro_sub1_O4vY7


python Scripts/addPermutations.py Results/InVitro_sub1_O2_RAvLib.csv 100000
python Scripts/addPermutations.py Results/InVitro_sub1_O3_RAvLib.csv 100000
python Scripts/addPermutations.py Results/InVitro_sub1_O4_RAvLib.csv 100000
python Scripts/addPermutations.py Results/InVitro_sub1_Y4_RAvLib.csv 100000
python Scripts/addPermutations.py Results/InVitro_sub1_Y5_RAvLib.csv 100000
python Scripts/addPermutations.py Results/InVitro_sub1_Y6_RAvLib.csv 100000
python Scripts/addPermutations.py Results/InVitro_sub1_Y7_RAvLib.csv 100000

python2 Scripts/plotDist.py InVitro_sub1_Lib_diversity Data/InVitro_sub1_48hr_counts.csv Data/InVitro_sub1_O2_RA_counts.csv Data/InVitro_sub1_O3_RA_counts.csv Data/InVitro_sub1_O4_RA_counts.csv Data/InVitro_sub1_Y4_RA_counts.csv Data/InVitro_sub1_Y5_RA_counts.csv Data/InVitro_sub1_Y6_RA_counts.csv Data/InVitro_sub1_Y7_RA_counts.csv -l Lib_48hr Old1 Old2 Old3 Young1 Young2 Young3 Young4


python2 Scripts/plotVolcano5.py Results/InVitro_sub1_results/InVitro_sub1_O3_RAvLib.csv -xl -3 3 -yl 150 -n SLIT2 RSPH3A SLC2A4 PAPPA2 SPP1 SOCS1 MBNL1 RBPMS IER2 SLC45A4 RBM4 LRMDA ECSCR NPB B3GALNT2 EIF4E RASSF8 CCDC15 EDC3 PXDC1 AA986860 CNOT3 SLC2A12 GRB7 FAM53C VMN1R55 VMN1R107 CHCHD5 MIF SNRPB2 C1QTNF5 SORL1 CDKN1A ZFP541 DBPHT2 DIS3L2 RPLP1 RPL32 DEFA3 DEFA17 ID4 HES1 IGF1R BMPR2 GFAP ASCL1 OLIG2 BMPR1A SOX2

python Scripts/plotGenes.py Results/InVitro_sub1_results/InVitro_sub1_O2_RAvLib.csv RSPH3A
