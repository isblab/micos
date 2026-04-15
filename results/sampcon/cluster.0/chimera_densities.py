import chimera
from chimera import openModels
from chimera import runCommand

# Names of proteins/domains for which we have created densities
prots = ['MIC10.0_C-term_61-78',
        'MIC10.0_matrix_TM2_37-60',
        'MIC10.0_N-term_TM1_1-36',
        'MIC10.1_C-term_61-78',
        'MIC10.1_matrix_TM2_37-60',
        'MIC10.1_N-term_TM1_1-36',
        'MIC13.0_central_24-78',
        'MIC13.0_C-term_79-118',
        'MIC13.0_N-term_TM_1-23',
        'MIC19.0_central_175-185',
        'MIC19.0_CHCH_186-227',
        'MIC19.0_coiled-coil_59-174',
        'MIC19.0_N-term_1-58',
        'MIC19.1_coiled-coil_59-174',
        'MIC19.1_N-term_1-58',
        'MIC19.2_CHCH_186-227',
        'MIC60.0_coiled-coil_410-582',
        'MIC60.0_link_583-626',
        'MIC60.0_LBS1_627-648',
        'MIC60.0_LBS2_649-682',
        'MIC60.0_mitophilin_683-758',
        'MIC60.1_coiled-coil_410-582',
        'MIC60.1_link_583-626',
        'MIC60.1_LBS1_627-648',
        'MIC60.1_LBS2_649-682',
        'MIC60.1_mitophilin_683-758',
        'MIC60.2_coiled-coil_410-582',
        'MIC60.3_coiled-coil_410-582'
        ]

# Set visualization thresholds
threshold = {'MIC10.0_C-term_61-78':0.024,          #10%
        'MIC10.0_matrix_TM2_37-60':0.0211,          #10%
        'MIC10.0_N-term_TM1_1-36':0.0297,           #10%
        'MIC10.1_C-term_61-78':0.0714,              #39%
        'MIC10.1_matrix_TM2_37-60':0.0183,          #10%
        'MIC10.1_N-term_TM1_1-36':0.0158,           #10%
        'MIC13.0_central_24-78':0.0348,             #10%
        'MIC13.0_C-term_79-118':0.0285,             #10%
        'MIC13.0_N-term_TM_1-23':0.0132,            #10%
        'MIC19.0_central_175-185':0.0209,           #10%
        'MIC19.0_CHCH_186-227':0.0567,              #10%
        'MIC19.0_coiled-coil_59-174':0.0431,        #20%
        'MIC19.0_N-term_1-58':0.0133,               #17.3%
        'MIC19.1_coiled-coil_59-174':0.0291,        #23%
        'MIC19.1_N-term_1-58':0.016  ,              #26.8%
        'MIC19.2_CHCH_186-227':0.055,               #10%
        'MIC60.0_coiled-coil_410-582':0.0582,       #10%
        'MIC60.0_link_583-626':0.0319,              #10%
        'MIC60.0_LBS1_627-648':0.0385,              #10%
        'MIC60.0_LBS2_649-682':0.0262,              #10%
        'MIC60.0_mitophilin_683-758':0.0458,        #10%
        'MIC60.1_coiled-coil_410-582':0.0469,       #10%
        'MIC60.1_link_583-626':0.0293,              #10%
        'MIC60.1_LBS1_627-648':0.0369,              #10%
        'MIC60.1_LBS2_649-682':0.0152,              #7%
        'MIC60.1_mitophilin_683-758':0.0514,        #10%
        'MIC60.2_coiled-coil_410-582':0.0443,       #10%
        'MIC60.3_coiled-coil_410-582':0.0478        #10%
        }

# Color of each protein/domain
col = {'MIC10.0_C-term_61-78':'#e06633',
        'MIC10.0_matrix_TM2_37-60':'#d47a00',
        'MIC10.0_N-term_TM1_1-36':'#ffa100',
        'MIC10.1_C-term_61-78':'#e06633',
        'MIC10.1_matrix_TM2_37-60':'#d47a00',
        'MIC10.1_N-term_TM1_1-36':'#ffa100',
        'MIC13.0_central_24-78':'#2ed766',
        'MIC13.0_C-term_79-118':'#2E8B57',
        'MIC13.0_N-term_TM_1-23':'#90ee90',
        'MIC19.0_central_175-185':'#a67573',
        'MIC19.0_CHCH_186-227':'#a67573',
        'MIC19.0_coiled-coil_59-174':'#fa8072',
        'MIC19.0_N-term_1-58':'#ffb5b5',
        'MIC19.1_coiled-coil_59-174':'#ccd88f918490',
        'MIC19.1_N-term_1-58':'#ffffd6dae0fb',
        'MIC19.2_CHCH_186-227':'#a67573',
        'MIC60.0_coiled-coil_410-582':'sky blue',
        'MIC60.0_link_583-626':'#4ea2ff',
        'MIC60.0_LBS1_627-648':'#006bff',
        'MIC60.0_LBS2_649-682':'#0000ff',
        'MIC60.0_mitophilin_683-758':'#0c13a4',
        'MIC60.1_coiled-coil_410-582':'dodger blue',
        'MIC60.1_link_583-626':'#4ea2ff',
        'MIC60.1_LBS1_627-648':'#006bff',
        'MIC60.1_LBS2_649-682':'#0000ff',
        'MIC60.1_mitophilin_683-758':'#0c13a4',
        'MIC60.2_coiled-coil_410-582':'sky blue',
        'MIC60.3_coiled-coil_410-582':'dodger blue'
        }


runCommand('set bgcolor white')
i=0

#Read localization density by component, both samples together
for p in prots:
    runCommand('open LPD_'+p+'.mrc')
    runCommand('volume #'+str(i)+' step 1 ')
    runCommand('volume #'+str(i)+' level '+str(threshold[p]))
    runCommand('color '+col[p]+' #'+str(i))
    i += 1
