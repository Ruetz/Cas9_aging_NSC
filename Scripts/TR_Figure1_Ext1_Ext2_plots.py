## Figure 1, Extended Data Figure 1 and 2, plots

import subprocess

##Figure 1f
subprocess.call('python2 Scripts/plotVolcano5.py Results/SC6_KpvLib.csv -f svg -yl 40 -xl -4.5 4.5 -n  CHCHD5 RPLP1 RPL32 MIF DEFA3 VMN1R107 SPTLC2 RSPH3A PXDC1 WDR35 INTU SLC2A4 SLC2A12 TRP53 SNAI1 FAM89B CDKN2A', shell=True)
subprocess.call('python2 Scripts/plotVolcano6.py Results/SC6_KpvLib.csv -yl 40 -xl -4.5 4.5 -n  CHCHD5 RPLP1 RPL32 MIF DEFA3 VMN1R107 SPTLC2 RSPH3A PXDC1 WDR35 INTU SLC2A4 SLC2A12 TRP53 SNAI1 FAM89B CDKN2A', shell=True)

#Figure 1g
subprocess.call('python2 Scripts/plotVolcano5.py Results/SC5_KpvLib.csv -f svg -yl 40 -xl -4 4 -n CREBBP CCDC169 TEX22 TOB2 ANAPC16 RAN SMC5 TCOF1 GINS2 CLP1 CDKN2A SNAI1 FAM89B TRP53', shell=True)
subprocess.call('python2 Scripts/plotVolcano6.py Results/SC5_KpvLib.csv -yl 40 -xl -4 4 -n CREBBP CCDC169 TEX22 TOB2 ANAPC16 RAN SMC5 TCOF1 GINS2 CLP1 CDKN2A SNAI1 FAM89B TRP53', shell=True)

#Extended Data Figure 1a/b
subprocess.call('python Scripts/plotDist.py SC345678_Day4_library_diversity Data/SC3_lib_counts.csv  Data/SC4_Y_Kp_counts.csv Data/SC3_O_Kp_counts.csv Data/SC5_Y_Kp_counts.csv Data/SC6_O_Kp_counts.csv Data/SC7_Y_Kp_counts.csv Data/SC8_O_Kp_counts.csv -l Library Young_1 Old_1 Young_2 Old_2 Young_3 Old_3', shell=True)
subprocess.call('python Scripts/plotDist.py SC345678_Day14_library_diversity Data/SC3_lib_counts.csv  Data/SC4_Y_RA_counts.csv Data/SC3_O_RA_counts.csv Data/SC5_Y_RA_counts.csv Data/SC6_O_RA_counts.csv Data/SC7_Y_RA_counts.csv Data/SC8_O_RA_counts.csv -l Library Young_1 Old_1 Young_2 Old_2 Young_3 Old_3', shell=True)

#Extended Data Figure 1e/f
subprocess.call('python2 Scripts/plotVolcano5.py Results/SC6_RAvLib.csv -f svg -yl 50 -xl -7.5 15 -n  CHCHD5 RPLP1 RPL32 MIF DEFA3 VMN1R107 SPTLC2 RSPH3A PXDC1 WDR35 INTU SLC2A4 SLC2A12 TRP53 SNAI1 FAM89B CDKN2A', shell=True)
subprocess.call('python2 Scripts/plotVolcano6.py Results/SC6_RAvLib.csv -yl 50 -xl -7.5 15 -n  CHCHD5 RPLP1 RPL32 MIF DEFA3 VMN1R107 SPTLC2 RSPH3A PXDC1 WDR35 INTU SLC2A4 SLC2A12 TRP53 SNAI1 FAM89B CDKN2A', shell=True)

subprocess.call('python2 Scripts/plotVolcano5.py Results/SC4_RAvLib.csv -f svg -yl 50 -xl -6 8 -n  CREBBP CCDC169 TEX22 TOB2 ANAPC16 RAN SMC5 TCOF1 GINS2 CLP1 CDKN2A SNAI1 FAM89B TRP53', shell=True)
subprocess.call('python2 Scripts/plotVolcano6.py Results/SC4_RAvLib.csv -yl 50 -xl -6 8 -n  CREBBP CCDC169 TEX22 TOB2 ANAPC16 RAN SMC5 TCOF1 GINS2 CLP1 CDKN2A SNAI1 FAM89B TRP53', shell=True)

#Extended Dta Figure 2a-c
subprocess.call('python2 Scripts/plotGenes_Teal.py Results/SC4_RAvLib.csv MDM2', shell=True)
subprocess.call('python2 Scripts/plotGenes_lavender.py Results/SC4_RAvLib.csv -x 15 -15 CREBBP', shell=True)
subprocess.call('python2 Scripts/plotGenes_Teal.py Results/SC4_RAvLib.csv PSMB1', shell=True)
subprocess.call('python2 Scripts/plotGenes_lavender.py Results/SC4_RAvLib.csv LRIG1', shell=True)
subprocess.call('python2 Scripts/plotGenes_lavender.py Results/SC6_KpvLib.csv Rsph3a', shell=True)
subprocess.call('python2 Scripts/plotGenes_Teal.py Results/SC6_KpvLib.csv Vmn1r107', shell=True)
