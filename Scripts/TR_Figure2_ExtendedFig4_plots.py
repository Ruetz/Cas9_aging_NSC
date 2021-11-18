##Figure 2 and Extended Figure 4 plots

import subprocess

# #Figure 2f
subprocess.call('python2 Scripts/plotGenes_Teal.py Results/Depleted_Sc1.csv -x 20 -20 VMN1R107', shell=True)
#Figure 2g
subprocess.call('python2 Scripts/plotGenes_Lavender.py Results/Glucose_HumanDisease_Sc1.csv SLC2A4 -x 15', shell=True)

#Extended Data Figure 4 a/b
subprocess.call('python2 Scripts/plotDist.py InVivo_T10_24hrs Data/IV9_T10S_24hr_counts.csv Data/T10_24hr_OB_counts.csv Data/T10_24hr_Cer_counts.csv Data/T10_24hr_Cort_counts.csv -l SVZ Olfactory_Bulb Cerebellum Cortex', shell=True)
subprocess.call('python2 Scripts/plotDist.py InVivo_T10_5wks Data/T10_5wk_SVZ_counts.csv Data/IV9_T10S_O1_counts.csv Data/T10_5wk_Cer_counts.csv Data/T10_5wk_Cort_counts.csv -l SVZ Olfactory_Bulb Cerebellum Cortex', shell=True)

#Extended Data Figure 4d-g
subprocess.call('python2 Scripts/plotGenes_Lavender.py Results/Top10_Sc2.csv -x 20 -20 RSPH3A', shell=True)
subprocess.call('python2 Scripts/plotGenes_Lavender.py Results/CytoplasmicRibonucleoproteinStructures_Sc2.csv -x 15 -15 EIF4E', shell=True)
subprocess.call('python2 Scripts/plotGenes_Teal.py Results/Depleted_Sc1.csv -x 30 -30 CHCHD5', shell=True)
subprocess.call('python2 Scripts/plotGenes_Lavender.py Results/PublishedRegulators_Sc1.csv BMPR2 -x 20', shell=True)
