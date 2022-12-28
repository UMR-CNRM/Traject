# -*- coding: utf-8 -*-
"""
Test case 2
Tracking extratropical cyclone Alex in two PEARP members using VDG algorithm
Required inputs :
    - inputdef file ./inputs/indef_testcase2.json
    - algodef file ./inputs/algo_testcase2.json
    - data files /cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/PEARP/grid*.grib (as specified in inputdef file)
    - Reference track of Alex ./inputs/ALEX_ANA.json

Expected outputs (in ./test2/):
    - track file track_test_case2*.json
    - figure file track_test_case2*.png

Only for members 008 et 010
- for basetime 2020100106, track for mb008 starts at 2020100112,  track for mb010 starts at 2020100106
- for basetime 2020100112, track for mb008 starts at 2020100112,  track for mb010 starts at 2020100112

"""

from traject import *

#Epygram environment
import epygram
epygram.init_env()
os.environ["ECCODES_SAMPLES_PATH"]=("/home/common/epygram/ext/eccodes/share/eccodes/samples")

#Directory
repin='./inputs/'
repout='./test2/'

#Prepare input files
#Extraction from hendrix using vortex
#Extract data fields in single grib files in /cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/PEARP/
#PrepareFiles(repin+"algo_testcase2.json",repin+"indef_testcase2.json",{"vortex":"","extract":"/cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/PEARP/","filter":"./tmp/case2/","dirout":repout},{'start':"2020100106",'final':"2020100112",'step':"06"},termtraj={'final':48,'step':3},members=[8,10])

#only the filter for testing
#PrepareFiles(repin+"algo_testcase2.json",repout+"indef_xtract.json",{"filter":"./tmp/case2/","dirout":repout},{'start':"2020100106",'final':"2020100112",'step':"06"},termtraj={'final':48,'step':3},members=[8,10])

#a/ Compute track using the extracted grib files
ltraj=track(repin+"algo_testcase2.json",repout+"indef_xtract.json",{'start':"2020100106",'final':"2020100112",'step':"06"},termtraj={'final':48,'step':3},members=[8,10],reftraj=repin+"ALEX_ANA.json",outfile=repout+'track_test_case2a_v'+str(traject_version)+'.json',plotfile=repout+'track_test_case2a_v'+str(traject_version)+'.png')

#b/ Compute track using the domain-filtered data
#ltraj=track(repin+"algo_testcase2.json",repout+"indef_filter.json",{'start':"2020100106",'final':"2020100112",'step':"06"},termtraj={'final':48,'step':3},reftraj=repin+"ALEX_ANA.json",outfile=repout+'track_test_case2b_v'+str(traject_version)+'.json',plotfile=repout+'track_test_case2b_v'+str(traject_version)+'.png',members=[8,10])

if False:
    #Read output track
    ltraj1=Read(repout+'track_test_case2a_v'+str(traject_version)+'.json')
    for ivi in range(len(ltraj1)):
        print("Reading the track afterwards :")
        print(ltraj1[ivi].__dict__)
    ltraj2=Read(repout+'track_test_case2b_v'+str(traject_version)+'.json')
    for ivi in range(len(ltraj2)):
        print("Reading the track afterwards :")
        print(ltraj2[ivi].__dict__)
