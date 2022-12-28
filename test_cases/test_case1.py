# -*- coding: utf-8 -*-
"""
Test case 1
Tracking extratropical cyclone Alex in an ARPEGE forecast using MSLP algorithm
Required inputs :
    - inputdef file ./inputs/idef_testcase1.json
    - algodef file ./inputs/algo_testcase1.json
    - data files /cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPfc/grid*.grib (as specified in inputdef file) that contain at least mslp field
    - Reference track of Alex ./test_cases/inputs/ALEX_ANA.json

Expected outputs (in ./test_cases/test1/):
    - track file track_test_case1.json
    - figure file track_test_case1.png

"""

from traject import *

#Epygram environment
import epygram
epygram.init_env()
os.environ["ECCODES_SAMPLES_PATH"]=("/home/common/epygram/ext/eccodes/share/eccodes/samples")

#Directory
repin='./inputs/'
repout='./test1/'

#Prepare input files
#Extraction from hendrix using vortex
#Extract mslp data in single grib files in /cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPfc/mslp/
#Domain filtering in ./tmp/
#PrepareFiles(repin+"algo_testcase1.json",repin+"indef_testcase1.json",{"vortex":"","extract":"/cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPfc/mslp/","filter":"./tmp/case1/","dirout":repout},{'start':"2020100112",'final':"2020100112",'step':"06"},termtraj={'final':48,'step':3})

#only filtering
#PrepareFiles(repin+"algo_testcase1.json",repout+"indef_xtract.json",{"filter":"./tmp/case1/","dirout":repout},{'start':"2020100112",'final':"2020100112",'step':"06"},termtraj={'final':48,'step':3})

#a/ Compute track using the extracted mslp grib files
ltraj=track(repin+"algo_testcase1.json",repout+"indef_xtract.json",{'start':"2020100112",'final':"2020100112",'step':"00"},termtraj={'final':48,'step':3},reftraj=repin+"ALEX_ANA.json",outfile=repout+'track_test_case1a_v'+str(traject_version)+'.json',plotfile=repout+'track_test_case1a_v'+str(traject_version)+'.png')
#exit()

#b/ Compute track using the domain-filtered data
#ltraj=track(repin+"algo_testcase1.json",repout+"indef_filter.json",{'start':"2020100112",'final':"2020100112",'step':"00"},termtraj={'final':48,'step':3},reftraj=repin+"ALEX_ANA.json",outfile=repout+'track_test_case1b_v'+str(traject_version)+'.json',plotfile=repout+'track_test_case1b_v'+str(traject_version)+'.png')

if False:
    #Read output tracks
    ltraj2=Read(repout+'track_test_case1a_v'+str(traject_version)+'.json')
    for ivi in range(len(ltraj2)):
        print("Reading the track afterwards :")
        print(ltraj2[ivi].__dict__)
        print("list of points:")
        for ivj in range(ltraj2[ivi].nobj):
            print(ltraj2[ivi].traj[ivj].__dict__)

    ltraj3=Read(repout+'track_test_case1b_v'+str(traject_version)+'.json')
    for ivi in range(len(ltraj3)):
        print("Reading the track afterwards :")
        print(ltraj3[ivi].__dict__)
        print("list of points:")
        for ivj in range(ltraj3[ivi].nobj):
            print(ltraj3[ivi].traj[ivj].__dict__)

