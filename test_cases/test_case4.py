# -*- coding: utf-8 -*-
"""
Test case 4

Apply Ayrault algorithm on ERA5 analyses during January 2020 over the Mediterranean - about 15 cyclones should be found
Required inputs :
    - inputdef file ./inputs/indef_testcase4.json
    - algodef file ./inputs/algo_testcase4.json
    - data files /cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPana/grid*.grib (as specified in inputdef file)

The varalgo and filtering parameters have been tuned by Benjamin Doiteau during his PhD.

Expected outputs (in ./test4/):
    - track file track_test_case4.json
    - figure file track_test_case4.png

"""

from traject import *

#Epygram environment
import epygram
epygram.init_env()
os.environ["ECCODES_SAMPLES_PATH"]=("/home/common/epygram/ext/eccodes/share/eccodes/samples")

#Directory
repin='./inputs/'
repout='./test4/'

timetraj={'start':"2020010100",'final':"2020020100",'step':"03"}
#Prepare input files - filter
#PrepareFiles(repin+"algo_testcase4.json",repin+"indef_testcase4.json",{"filter":"./tmp/case4/","dirout":repout},timetraj)

#a/ Compute track using the native nc files
ltraj=track(repin+"algo_testcase4.json",repin+"indef_testcase4.json",timetraj,outfile=repout+'track_testcase4a_v'+str(traject_version)+'.json',plotfile=repout+'track_testcase4a_v'+str(traject_version)+'.png')

#b/ Compute track using the domain-filtered data
#ltraj=track(repin+"algo_testcase4.json",repout+"indef_filter.json",timetraj,outfile=repout+'track_testcase4b_v'+str(traject_version)+'.json',plotfile=repout+'track_testcase4b_v'+str(traject_version)+'.png')

if False:
    ltraj=Read(repout+'track_testcase4a_v'+str(traject_version)+'.json')
    for ivi in range(len(ltraj)):
        print("Reading the track afterwards :")
        print(ltraj[ivi].__dict__)
    ltraj2=Read(repout+'track_testcase4b_v'+str(traject_version)+'.json')
    for ivi in range(len(ltraj)):
        print("Reading the track afterwards :")
        print(ltraj2[ivi].__dict__)

