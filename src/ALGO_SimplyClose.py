#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Traject/SimplyClose algorithm

https://github.com/UMR-CNRM/Traject
Traject is disributed under Cecill-C license
(For rights and obligations, see LICENSE.txt)

UMR-CNRM, Meteo-France & CNRS
"""
# =============================================================================
#                           ALGORITHME DE SUIVI
# =============================================================================

''' A simple algorithm that tracks an object from an initial minimum value of the input
field (for instance mslp), and from time to time by taking the closest
minimim value at a distance below max_dist'''

#------------------------------------------------------------------------------
                                # Importations
#------------------------------------------------------------------------------

import Inputs,Tools
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from datetime import timedelta, datetime
from traject import *

#------------------------------------------------------------------------------
                                # Variables globales 
#------------------------------------------------------------------------------
                                
def track(algo,indf,linst,lfile,**kwargs):

    #Initialisation of domtraj, diag_parameter, deltat depending on the input parameters
    domtraj, res, deltat, Hn, diag_parameter, subnproc = Tools.init_algo(algo,indf,linst,lfile,**kwargs)

    #Tracking variables used by the algorithm
    max_dist = max(algo.varalgo["maxdist"],algo.varalgo["maxspeed"]*deltat)
                     #maximum distance (km) between 2 points that belong to a track
                     #maxdist is in km, maxspeed in km/h
    trackpar, signtrack = Tools.get_parsign(algo.specfields["track"], Hn)
                     #parameter whose extremum is followed and associated sign (-1=min ; 1=max)
    track_parameter = [trackpar] #List of parameters that are declared and used in every objects

    ltraj=[]

    if "basetime" in kwargs:
        basetime=kwargs["basetime"]
    else:
        basetime=datetime.strftime(linst[0],format=time_fmt)

    if "last_start" in algo.varalgo: #List of fc instants where to search for starting point
        linst2 = [inst for inst in linst if Tools.comp_difftime(linst[0],inst)<=algo.varalgo["last_start"]]
    else:
        linst2 = linst

    ss=0.0
    rd=0.0
    if "ss" in algo.varalgo:
        ss=algo.varalgo["ss"] #Size of the square for search of extremas to compute diagnostics (in degrees)
    if "rd" in algo.varalgo:
        rd=algo.varalgo["rd"] #Radius of the circle for search of extremas to compute diagnostics (in km)

    #Read reference track and find the list of tracks valid at the first instant and in the domain
    if "reftraj" in kwargs:
        lref = Tools.get_reftraj([linst2[0]],domtraj,kwargs["reftraj"])
    else:
        lref=[]
    if len(lref)==0:
        print("WARNING - No reference track found could be found (required for SimplyClose tracking algorithm)")

    #Loops on points found in the reference tracks
    for iref in range(len(lref)):
        if "traj" in locals():
            del traj
        reftraj=lref[iref]

        it0=-1
        gook=False
        while not gook and it0<len(linst2)-1: #Search for starting point
            it0=it0+1
            print(linst2[it0])
            refobj, gook =reftraj.find_inst(linst2[it0]) #Attention, delta t différent !!
            fld = Inputs.extract_data(lfile[it0],linst[it0],indf,trackpar,domtraj,res[trackpar],basetime,subnproc)
            lon,lat,val,dist,gook = Tools.find_locextr(signtrack,fld,refobj[0].lonc,refobj[0].latc)
            gook = gook and(dist<=max_dist)
            #print("Reference: ",linst[0],refobj[0].lonc,refobj[0].latc)
            #print(gook,lon,lat,val,dist)

        #------------------------------------------------------------------------------
                        # First step: Start the track from the closest point to reftraj
        #------------------------------------------------------------------------------
        if gook:
            print("A starting point has been found at time : ",linst2[it0])
            #Creation of the track
            traj = DefTrack(algo.classobj,basetime=basetime)
            traj.name = reftraj.name
            objectm0 = DefObject(algo.classobj, [], track_parameter,lonc=lon,latc=lat,time=linst2[it0])
            Tools.make_diags(diag_parameter,objectm0,ss,rd,lfile[0],linst2[0],indf,domtraj,Hn,res,basetime,lon,lat,subnproc)
            traj.add_obj(objectm0)

        #------------------------------------------------------------------------------
                                    # TIME LOOP
        #------------------------------------------------------------------------------
                                    
        it=it0

        while gook and it<len(linst)-1:
            it = it + 1
            print("Looking for a point at instant ", linst[it])
            instm1=linst[it-1]
            inst=linst[it]
            objm1,fnd=traj.find_inst(instm1)

            fld = Inputs.extract_data(lfile[it],linst[it],indf,trackpar,domtraj,res[trackpar],basetime,subnproc)
            lon,lat,val,dist,gook = Tools.find_locextr(signtrack,fld,objm1[0].lonc,objm1[0].latc)

            gook=(dist<=max_dist)
            if gook:
                #CREATE POINT ON THE TRACK
                objectm = DefObject(algo.classobj, [], track_parameter,lonc=lon,latc=lat,time=inst)

                #DIAGNOSTIC PARAMETERS
                Tools.make_diags(diag_parameter,objectm,ss,rd,lfile[it],linst[it],indf,domtraj,Hn,res,basetime,lon,lat,subnproc)
                traj.add_obj(objectm)

        if "traj" in locals():
            ltraj.append(traj)

    return ltraj
