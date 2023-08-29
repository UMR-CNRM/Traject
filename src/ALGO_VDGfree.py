#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Traject/VDG algorithm

https://github.com/UMR-CNRM/Traject
Traject is disributed under Cecill-C license
(For rights and obligations, see LICENSE.txt)

UMR-CNRM, Meteo-France & CNRS

"""
# =============================================================================
#                           ALGORITHME DE SUIVI
# =============================================================================

''' VDG algorithm
inspired by Van Der Grijn, 2002, Tech Memo 386 ECMWF '''

#------------------------------------------------------------------------------
                                # Importations
#------------------------------------------------------------------------------

import Inputs,Tools
#import cartopy.crs as ccrs
#import cartopy.feature as cfeature
#import matplotlib.pyplot as plt
from datetime import timedelta, datetime
from traject import *
import concurrent.futures

DEBUGGING=False

def track(algo,indf,linst,lfile,**kwargs):

    #Initialisation of domtraj, res, deltat depending on the input parameters
    domtraj, res, deltat, Hn, diag_parameter, lparam, subnproc = Tools.init_algo(algo,indf,linst,lfile,**kwargs)

    #Tracking variables used by the algorithm
    parfilt=algo.parfilt #Effective resolution of the fields
    ss=algo.varalgo["ss"] #Size of the square for search of minimas (in degrees)
    w = algo.varalgo["w"] #Parameter to balance steering flow and past movement to find the guess

    #uvmean_box ?
    if "uvmean_box" in algo.varalgo:
        uvmean_box=algo.varalgo["uvmean_box"]
    else:
        uvmean_box=0

    filtapply=Tools.check_filtapply(algo,indf,linst,lfile,res)
             # 1 if filtering must be applied in the routine ; 0 if the input data has already been filtered
    steering_levels = algo.varalgo["steering_levels"]
             #hPa levels taken into account for the computation of steering flow

    trackpar, signtrack = Tools.get_parsign(algo.specfields["track"], Hn)
                     #parameter whose extremum is tracked and associated sign (-1=min ; 1=max)
                     #for cyclone : rv850

    pairpar, signpair = Tools.get_parsign(algo.specfields["pairing"], Hn)
                     #parameter whose extremum is paired to and associated sign (-1=min ; 1=max)
                     #for cyclone : mslp

    thr_track = algo.varalgo["thr_track"] #Minimum value of track parameter (rv850 for instance) to keep the track at the end

    thr_pair = algo.varalgo["thr_pairing"] #Minimum value of paired core

    track_parameter = [trackpar,pairpar,"olon","olat","u_speed","v_speed","u_steer","v_steer"] #olon, olat : coordinates of the trackpar extrema, used for tracking
    if "thr_core" in algo.varalgo:
        thr_core = algo.varalgo["thr_core"] #Minimum value of track parameter (btir for instance) to detect the initial cores
    else:
        thr_core = thr_track

    rd=0.0
    if "rd" in algo.varalgo:
        rd=algo.varalgo["rd"] #Radius of the circle for search of extremas to compute diagnostics (in km)

    radmax = algo.varalgo["radmax"] #Radius (km) inside which only one minimum is kept (the highest vorticity value) 

    ltraj=[]

    if "basetime" in kwargs:
        basetime=kwargs["basetime"]
    else:
        basetime=datetime.strftime(linst[0],format=time_fmt)

    if "last_start" in algo.varalgo: #List of fc instants where to search for starting point
        linst2 = [inst for inst in linst if Tools.comp_difftime(linst[0],inst)<=algo.varalgo["last_start"]]
    else:
        linst2 = linst

##################*

    it=0
    exclude=False
    ltraj0=[]
    while it<len(linst)-1:
        print("Instant - ",linst[it])
        
        if it<len(linst2)-1:
            dict_fld = Tools.load_fld(lparam,lfile[it],linst[it],indf,algo,domtraj,res,basetime,subnproc,parfilt=parfilt,filtapply=filtapply) #input fields at the given instant
            #Finding all extremas (only for instants in linst2)
            lobj = Search_allcores(algo.classobj,dict_fld=dict_fld,inst=linst[it],trackpar=trackpar,signtrack=signtrack,track_parameter=track_parameter,dist=radmax, thr=thr_core)

            print("Testing "+str(len(lobj))+" initial objects - Parallelisation on subnproc="+str(subnproc)+" threads")
            execore = concurrent.futures.ThreadPoolExecutor(max_workers=subnproc)
            outf2 = []

            #Read input fields once
            trackfld=[]
            pairfld=[]
            condfld=[]

            for obj in lobj:
                #Look for potential initial cores
                outf2.append(execore.submit(init_object,obj,dict_fld,linst[it],pairpar,Hn,algo,domtraj,track_parameter,diag_parameter,ss,thr_pair,True,True,ltraj0,radmax))

            for out in outf2:
                obj, obj2, isclose = out.result()
                if not isclose:
                    #Create track
                    print("A starting point has been found at time : ",linst[it])
                    print("Point found:",obj2.lonc,obj2.latc)
                    traj = DefTrack(algo.classobj,basetime=basetime)

                    #Initialisations of obj2 variables and addition to traj
                    obj2.traps["olon"]=obj.lonc
                    obj2.traps["olat"]=obj.latc
                    u_steer, v_steer = Tools.comp_steering(obj2,steering_levels,uvmean_box,dict_fld,pos="o")

                    obj2.traps["u_steer"] = u_steer
                    obj2.traps["v_steer"] = v_steer
                    obj2.traps["u_speed"] = u_steer #we cannot determine a movement at step 0, so we take wind_steer
                    obj2.traps["v_speed"] = v_steer
                    obj2.traps[trackpar] = obj.traps[trackpar]

                    traj.add_obj(obj2)
                    ltraj0.append(traj)

        #Add point at next step for all tracks (for all instants)
        it = it + 1
        dict_fld = Tools.load_fld(lparam,lfile[it],linst[it],indf,algo,domtraj,res,basetime,subnproc,parfilt=parfilt,filtapply=filtapply) #input fields at the given instant
        for traj in ltraj0:
            print("Looking for a point at instant ", linst[it])
            instm1=linst[it-1]
            inst=linst[it]
            objm,fnd=traj.find_inst(instm1) #Object at the previous instant
            if fnd:
                objm1=objm[0]

                #Computation of a guess (advection by the steering flow and displacement speed at objm1)
                #We advect olon, olat
                obj_guess = objm1.advect(inst, (1-w)*objm1.traps["u_steer"]+ w*objm1.traps["u_speed"],
                    (1-w)*objm1.traps["v_steer"]+ w*objm1.traps["v_speed"],pos="o")

                obj = obj_guess.search_core(dict_fld,linst[it],trackpar,Hn,track_parameter,[],ss=ss,thr_param=thr_core)
                obj2 = None
                if obj is not None: # Pairing
                    obj2 = obj.search_core(dict_fld,linst[it],pairpar,Hn,track_parameter,diag_parameter,ss=ss,thr_param=thr_pair,pairing=True,smooth=True)
                    if obj2 is not None:
                        isok, exclude = obj2.conditiontype(dict_fld,algo,domtraj,init=False)

                if obj is not None and obj2 is not None and not exclude:
                    #Initialisations of obj2 variables and addition to traj
                    obj2.traps["olon"]=obj.lonc
                    obj2.traps["olat"]=obj.latc
                    u_steer, v_steer = Tools.comp_steering(obj2,steering_levels,uvmean_box,dict_fld,pos="o")
                    u_speed, v_speed = obj2.comp_mvt(objm1,pos="o")
                    obj2.traps["u_steer"] = u_steer
                    obj2.traps["v_steer"] = v_steer
                    obj2.traps["u_speed"] = u_speed
                    obj2.traps["v_speed"] = v_speed
                    obj2.traps[trackpar] = obj.traps[trackpar]
                    traj.add_obj(obj2)

    #Finalisation
    stt=(signtrack==1)
    for traj in ltraj0:
        bmax=False
        maxv=0.0

        if traj.nobj>0:
            bmax=False
            maxv=0.0
            for obj in traj.traj:
                print("Final "+trackpar+":",obj.traps[trackpar])
                #condition on min/max value
                if (stt and obj.traps[trackpar]>=thr_track) or ((not stt) and obj.traps[trackpar]<=thr_track):
                    maxv=obj.traps[trackpar]
                    bmax=True
        if bmax:
            ltraj.append(traj)
        else:
            print("This trajectory does not reach the tracking threshold - We skip it")

    return ltraj

def init_object(obj,dict_fld,inst,param,Hn,algo,domtraj,track_parameter,diag_parameter,ss,thr_pair,pairing,smooth,ltraj0,radmax):

    isclose=True
    obj2 = obj.search_core(dict_fld,inst,param,Hn,track_parameter,diag_parameter,ss=ss,thr_param=thr_pair,pairing=pairing,smooth=smooth)

    if obj2 is not None:
        isok, exclude = obj2.conditiontype(dict_fld,algo,domtraj,init=True)
        if obj2 is not None and not exclude:
            #Test if obj2 is close to an pre-existing traj
            isclose=False
            for traj in ltraj0:
                objt,fnd=traj.find_inst(inst)
                if fnd:
                    objt1=objt[0]
                    if Tools.comp_length(obj2.lonc,obj2.latc,objt1.lonc,objt1.latc)<radmax:
                        isclose=True

    return obj, obj2, isclose
