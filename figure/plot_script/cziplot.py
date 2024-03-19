level0_order=['Exc','Inh']
level1_order=['Exc_upper', 'Exc_inter', 'Exc_deep','Inh_CGE', 'Inh_MGE']
  
level2_order=[
              'L2-4IT_CUX2','L3-5IT_RORB_PLCH1','L4-5IT_RORB_TSHZ2','L4-5IT_RORB_LRRK1','L4-5IT_RORB_ARHGAP15', 
              'L56NP_TLE4_TSHZ2','L6IT_THEMIS_LINC00343','L6IT_THEMIS_CUX1', 'L6CT_TLE4_FAM95C','L6b_TLE4_NXPH4',
              'CGE_LAMP5', 'CGE_LAMP5_LHX6','CGE_VIP','CGE_ADARB2_ADAM33','CGE_PAX6',
              'MGE_SST','MGE_SST_CLMP','MGE_PVALB',
             ]
level2_order_full=[
              'L2-4IT_CUX2','L3-5IT_RORB_PLCH1','L4-5IT_RORB_TSHZ2','L4-5IT_RORB_LRRK1','L4-5IT_RORB_ARHGAP15', 
              'L56NP_TLE4_TSHZ2','L5ET_FEZF2_ADRA1A','L6IT_THEMIS_LINC00343','L6IT_THEMIS_CUX1', 'L6CT_TLE4_FAM95C','L6b_TLE4_NXPH4',
              'CGE_LAMP5', 'CGE_LAMP5_LHX6','CGE_VIP','CGE_ADARB2_ADAM33','CGE_PAX6',
              'MGE_SST','MGE_SST_CLMP','MGE_SST_NPY','MGE_PVALB', 'MGE_PVALB_COL15A1',
              'Glia_Astro','Glia_Oligo', 'Glia_Micro'
             ]

level3_order=[
            'L2-4IT_CUX2_LINC01331','L3-5IT_RORB_PLCH1','L4-5IT_RORB_TSHZ2','L4-5IT_RORB_GSN','L4-5IT_RORB_WHRN','L4-5IT_RORB_ARHGAP15',
            'L56NP_TLE4_TSHZ2','L6IT_THEMIS_LINC00343','L6IT_THEMIS_CUX1','L6CT_TLE4_FAM95C','L6b_TLE4_NXPH4',
            'CGE_LAMP5_NDNF','CGE_LAMP5_FREM1','CGE_LAMP5_LHX6',
            'CGE_VIP_ZBTB20','CGE_VIP_FGD5','CGE_VIP_DPF3','CGE_ADARB2_ADAM33','CGE_PAX6',
            'MGE_SST_RAB31','MGE_SST_CDH12','MGE_SST_CLMP',
            'MGE_PVALB_MYO5B',
             ]


sample_order=[
              'YM1B','YM2A','YM2B','YM3A', 'YM3B',
              'YF1A','YF1B','YF2A','YF2B',
              'AM1A', 'AM1B', 'AM2A','AM3A', 'AM3B',
              'AF1A', 'AF1B','AF2A','AF2B', 'AF3A','AF3B',
             ]
sample_order_rna=[
              'YM1B','YM2A','YM2B','YM3A', 'YM3B',
              'YF1A','YF1B','YF2A','YF2B',
              'AM1A', 'AM1B', 
              'AF1A', 'AF1B','AF2B', 'AF3A','AF3B',
             ]

donor_order=['YM1', 'YM2','YM3','YF1','YF2','AM1', 'AM2','AM3','AF1', 'AF2','AF3']
donor_order_rna=['YM1', 'YM2','YM3','YF1','YF2','AM1','AF1', 'AF2','AF3']
demo_order=['YM','YF','AM','AF']
age_order = ['Y','A']
sex_order = ['M','F']
 
sample_palette={'AF1A':'#ED004B','AF1B':'#ED004B','AF2A':"#F500AA",'AF2B':'#F500AA','AF3A':'#9700AE','AF3B':'#9700AE',
                'AM1A':'#00d5ff','AM1B':'#00d5ff','AM2A':'#0090FE','AM2B':'#0090FE',"AM3A":'#001BD0','AM3B':'#001BD0',
                'YF1A':'#FFD000','YF1B':'#FFD000','YF2A':"#FF9200",'YF2B':'#FF9200',
                'YM1A':'#92F000','YM1B':'#92F000','YM2A':'#0FA100','YM2B':'#0FA100',"YM3A":'#005F18','YM3B':'#005F18'
               }
donor_palette={'AF1':'#ED004B','AF2':"#F500AA",'AF3':'#9700AE',
                    'AM1':'#00d5ff','AM2':'#0090FE',"AM3":'#001BD0',
                    'YF1':'#FFD000','YF2':"#FF9200",
                    'YM1':'#92F000','YM2':'#0FA100',"YM3":'#005F18'
                   }
demo_palette={'YM':'#8ac0e5','YF':'#f5a4c4','AM':'#47719b','AF':'#ba4869'}
age_palette={'A':'#384fa2','Y':'#f8bb1f'}
sex_palette={'F':'#F16090','M':'#61A3D7'}
                
level2_palette={'L2-4IT_CUX2':'#64bc46','L3-5IT_RORB_PLCH1':'#31b6b1',
                'L4-5IT_RORB_ARHGAP15':'#24a371','L4-5IT_RORB_TSHZ2':'#6ff2d2','L4-5IT_RORB_LRRK1':'#088523',
                'L56NP_TLE4_TSHZ2':'#298b98','L5ET_FEZF2_ADRA1A':'#0c4657',
                'L6IT_THEMIS_LINC00343':'#bdc123','L6IT_THEMIS_CUX1':'#5e4382',
                'L6CT_TLE4_FAM95C':'#1f6466','L6b_TLE4_NXPH4':'#234466',
                'CGE_LAMP5':'#8e5a63','CGE_LAMP5_LHX6':'#632727',
                'CGE_VIP':'#e891d0','CGE_ADARB2_ADAM33':'#b36c76','CGE_PAX6':'#f4d3f5',
                'MGE_PVALB':'#ea2c44', 'MGE_PVALB_COL15A1':'#d82764',
                'MGE_SST':'#d97807','MGE_SST_CLMP':'#e89c04','MGE_SST_NPY':'#ebd716',
                'Glia_Oligo':'lightgrey','Glia_Astro':'darkgrey','Glia_Micro':'grey'}


celltype_simdic={
 'L2-4IT_CUX2':'L2-4IT',
 'L3-5IT_RORB_PLCH1':'L3-5IT',
 'L4-5IT_RORB_TSHZ2':'L4-5IT TSHZ2',
 'L4-5IT_RORB_LRRK1':'L4-5IT LRRK1',
 'L4-5IT_RORB_ARHGAP15':'L4-5IT ARHGAP15',
 'L56NP_TLE4_TSHZ2':'L56NP',
 'L6IT_THEMIS_LINC00343':'L6IT LINC00343',
 'L6IT_THEMIS_CUX1':'L6IT CUX1',
 'L6CT_TLE4_FAM95C':'L6CT',
 'L6b_TLE4_NXPH4':'L6b',
 'CGE_LAMP5':'CGE LAMP5',
 'CGE_LAMP5_LHX6':'CGE LAMP5 LHX6',
 'CGE_VIP':'CGE VIP',
 'CGE_ADARB2_ADAM33':'CGE ADAM33',
 'CGE_PAX6':'CGE PAX6',
 'MGE_SST':'MGE SST',
 'MGE_SST_CLMP':'MGE SST CLMP',
 'MGE_PVALB':'MGE PVALB'
}


