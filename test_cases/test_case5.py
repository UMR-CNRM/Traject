# -*- coding: utf-8 -*-
"""
Test case 5

Apply VDG algorithm on ARPEGE forecast on a tropical cyclone (Cheneso, January 2023).
Required inputs :
    - inputdef file ./inputs/indef_testcase5.json
    - algodef file ./inputs/algo_testcase5.json
    - reference input file (that has been processed from https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.last3years.list.v04r00.csv)
    - data files /cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPana/grid*.grib (as specified in inputdef file)

Expected outputs (in ./test5/):
    - track file track_test_case5.json
    - figure file track_test_case5.png

A single track should be found.
"""

from traject import *

#Epygram environment
import epygram
epygram.init_env()
os.environ["ECCODES_SAMPLES_PATH"]=("/home/common/epygram/ext/eccodes/share/eccodes/samples")

#Directory
repin='./inputs/'
repout='./test5/'

timeref={'start':"2023011500",'final':"2023013000",'step':"06"}
diags=["mslp_min","ff10m_max"]

#Reference input file (data from https://www.ncei.noaa.gov/products/international-best-track-archive)
ibfile=repin+"ibtracs_CHENESO.json"

#If you want to apply this case on another cyclone, then get input reference file from IbTracs (https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.last3years.list.v04r00.csv) and uncomment the following lines:
#domtraj={"lonmin":-30.0,"lonmax":90.0,"latmin":-40.0,"latmax":-5.0} #same as indef_testcase5.json
#lref=ConvertIBTRACS("./ibtracs.last3years.list.v04r00.csv",timeref,domtraj,diags=diags)
#cycname="CHENESO"
#mslp_thr=1010.0 #Threshold (hPa) to keep the track
#ibfile="ibtracs_"+cycname+".json"
#for traj in Select(lref,{"name":cycname}):
#    print(traj.name,traj.basetime)
#    mini = traj.tmin("mslp_min")/100.0
#    if mini<mslp_thr:
#        print(traj.name, "Minimum mslp: ",mini, " hPa -- Kept ")
#        Write(traj,ibfile)


########################################
#Tracking ARPEGE forecasts at different dates
timefc={'start':"2023011800",'final':"2023011900",'step':"06"}

#Prepare input files
#Extraction from hendrix using vortex
#Extract mslp data in single grib files in /cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPfc/vdg/
#Domain filtering in ./tmp/
#PrepareFiles(repin+"algo_testcase5.json",repin+"indef_testcase5.json",{"vortex":"","extract":"/cnrm/recyf/NO_SAVE/Data/users/plu/TRAJECT/ARPfc/vdg/","filter":"./tmp/case5/","dirout":repout},timefc,termtraj={'final':48,'step':3})
#Clean_should()

#only filtering
#PrepareFiles(repin+"algo_testcase5.json",repout+"indef_xtract.json",{"filter":"./tmp/case5/","dirout":repout},timefc,termtraj={'final':48,'step':3})

#a/ Compute track using the extracted grib files
ltraj=track(repin+"algo_testcase5.json",repout+"indef_xtract.json",timefc,termtraj={'final':48,'step':3},reftraj=ibfile,outfile=repout+'track_test_case5a_v'+str(traject_version)+'.json',plotfile=repout+'track_test_case5a_v'+str(traject_version)+'.png')

#b/ Compute track using the extracted grib files
#ltraj=track(repin+"algo_testcase5.json",repout+"indef_filter.json",timefc,termtraj={'final':48,'step':3},reftraj=ibfile,outfile=repout+'track_test_case5b_v'+str(traject_version)+'.json',plotfile=repout+'track_test_case5b_v'+str(traject_version)+'.png')

if False:
    ltraj=Read(repout+'track_testcase5a_v'+str(traject_version)+'.json')
    for ivi in range(len(ltraj)):
        print("Reading the track afterwards :")
        print(ltraj[ivi].__dict__)
    ltraj2=Read(repout+'track_testcase5b_v'+str(traject_version)+'.json')
    for ivi in range(len(ltraj)):
        print("Reading the track afterwards :")
        print(ltraj2[ivi].__dict__)

#Plotting tracks
from visu_traject import *
lref=Read(ibfile)
lname=[traj.name for traj in lref]
fcfile=repout+'track_test_case5a_v'+str(traject_version)+'.json'
lfc=Read(fcfile)
print(lname)

for name in lname:
    ltraj=Select(lref,{"name":name})
    ltraj.extend(Select(lfc,{"name":name}))
    legend=["Ref"]
    for ivi in range(1,len(ltraj)):
        legend.append(ltraj[ivi].basetime)
    fig=map_plot(ltraj,diag="",typeplot="line",colormap="gist_rainbow",opt="zoom",leg=legend)
    for obj in ltraj[0].traj:
        plt.plot(obj.lonc,obj.latc,'o',color='red')
    fig.savefig(repout+"tracks_"+name+".png")
    #Plot mslp
    fig2=time_plot(ltraj,diag="mslp_min_p",typeplot="line",colormap="gist_rainbow",leg=legend)
    fig2.savefig(repout+"mslp_"+name+".png")
    plt.close("all")

