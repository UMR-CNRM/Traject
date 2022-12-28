# -*- coding: utf-8 -*-
"""
Test case 3

Apply Ayrault algorithm on ARPEGE analyses during the period of Alex (Sept. & Oct. 2020)
Required inputs :
    - inputdef file ./inputs/indef_testcase3.json
    - algodef file : ./inputs/algo_testcase3.json
    - data files /cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPana/grid*.grib (as specified in inputdef file)

Expected outputs (in ./inputs/):
    - track file track_test_case3.json
    - figure file track_test_case3.png

"""

from traject import *

#Epygram environment
import epygram
epygram.init_env()
os.environ["ECCODES_SAMPLES_PATH"]=("/home/common/epygram/ext/eccodes/share/eccodes/samples")

#Directory
repin='./inputs/'
repout='./test3/'

#Prepare input files
#Extraction from hendrix using vortex
#Extract data fields in single grib files in /cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPana/
#PrepareFiles(repin+"algo_testcase3.json",repin+"indef_testcase3.json",{"vortex":"","extract":"/cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPana/","filter":"./tmp/case3/","dirout":repout},{'start':"2020092500",'final':"2020100500",'step':"06"})

#only the filter for testing
#PrepareFiles(repin+"algo_testcase3.json",repout+"indef_xtract.json",{"filter":"./tmp/case3/","dirout":repout},{'start':"2020092500",'final':"2020100500",'step':"06"})

#a/ Compute track using the extracted grib files
ltraj=track(repin+"algo_testcase3.json",repout+"indef_xtract.json",{'start':"2020092500",'final':"2020100500",'step':"06"}, outfile=repout+'track_test_case3a_v'+str(traject_version)+'.json',plotfile=repout+'track_test_case3a_v'+str(traject_version)+'.png')

#b/ Compute track using the domain-filtered data
#ltraj=track(repin+"algo_testcase3.json",repout+"indef_filter.json",{'start':"2020092500",'final':"2020100500",'step':"06"},outfile=repout+'track_test_case3b_v'+str(traject_version)+'.json',plotfile=repout+'track_test_case3b_v'+str(traject_version)+'.png')

if False:
    #Read output track
    #ltraj1=Read(repout+'track_test_case3a.json')
    #for ivi in range(len(ltraj1)):
    #    print("Reading the track afterwards :")
    #    print(ltraj1[ivi].__dict__)
    ltraj2=Read(repout+'track_test_case3b_v'+str(traject_version)+'.json')
    for ivi in range(len(ltraj2)):
        print("Reading the track afterwards :")
        print(ltraj2[ivi].__dict__)

