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

    if "ss_pair" in algo.varalgo: #Optional ss value for pairing
        ss_pair = algo.varalgo["ss_pair"]
    else:
        ss_pair = ss

    if "mintlen" in algo.varalgo:
        mintlen=algo.varalgo["mintlen"] #Minimum length required (in hours) to keep the track at the end
    else:
        mintlen = 0

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

    dico_ker={} #Dictionary that contains the object candidates at the different instants
    dict_fld={} #Dictionary of fields

#=========================================================================#
    ''' Step 1: extraction of candidate points at every time step (fill in dico_ker) '''
#=========================================================================#    

    it=0
    while it<len(linst) and Inputs.check_file(lfile[it],indf,trackpar,trackpar):

        #We could introduce here PARALLELISATION (PoolProcess with subnproc) ... then report it in ALGO_AYRAULT
        inst=linst[it]

        dict_fld[str(it)] = Tools.load_fld(lparam,lfile[it],linst[it],indf,algo,domtraj,res,basetime,subnproc,parfilt=parfilt,filtapply=filtapply) #input fields at the given instant

        #Finding all extremas
        lobj0 = Search_allcores(algo.classobj,dict_fld=dict_fld[str(it)],inst=inst,trackpar=trackpar,signtrack=signtrack,track_parameter=track_parameter,dist=radmax, thr=thr_core)

        #END PARALLELISATION - Output: lobj

        dico_ker[str(it)] = lobj0

        it = it + 1

    it=0
    exclude=False
    ltraj0=[]

#=========================================================================#
    ''' Step 2: Loop on time steps '''
#=========================================================================#    

    while it<len(linst)-1:
        print("Instant - ",linst[it])

        #Read input fields once
        trackfld=[]
        pairfld=[]
        condfld=[]

        #LOOK FOR STARTING TRACKS
        for obj in dico_ker[str(it)]:
            #Look for potential initial cores
            obj, obj2, isclose = init_object(obj,dict_fld[str(it)],linst[it],pairpar,Hn,algo,domtraj,track_parameter,diag_parameter,ss_pair,thr_pair,True,True,ltraj0,radmax)

            if not isclose:
                #Create track
                print("A starting point has been found at time : ",linst[it])
                print("Point found:",obj2.lonc,obj2.latc)
                traj = DefTrack(algo.classobj,basetime=basetime)

                #Initialisations of obj2 variables and addition to traj
                obj2.traps["olon"]=obj.lonc
                obj2.traps["olat"]=obj.latc
                u_steer, v_steer = Tools.comp_steering(obj2,steering_levels,uvmean_box,dict_fld[str(it)],pos="o")

                obj2.traps["u_steer"] = u_steer
                obj2.traps["v_steer"] = v_steer
                obj2.traps["u_speed"] = u_steer #we cannot determine a movement at step 0, so we take wind_steer
                obj2.traps["v_speed"] = v_steer
                obj2.traps[trackpar] = obj.traps[trackpar]

                traj.add_obj(obj2)
                ltraj0.append(traj)

        #Add point at next step for all tracks 
        it = it + 1
        lobj = [] #Potential following object for each track
        lobj2 = [] 

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

                obj = obj_guess.search_kernel(dico_ker[str(it)], ss) #Il ne trouve pas le point ici : pb ?
                obj2 = None

                if obj is not None: # Pairing
                    obj2 = obj.search_core(dict_fld[str(it)],linst[it],pairpar,Hn,track_parameter,diag_parameter,ss=ss_pair,thr_param=thr_pair,pairing=True,smooth=True)
                    if obj2 is not None:
                        isok, exclude = obj2.conditiontype(dict_fld[str(it)],algo,domtraj,init=False)

                if obj is not None and obj2 is not None and not exclude:
                    lobj.append(obj)
                    lobj2.append(obj2)
                else:
                    lobj.append(None)
                    lobj2.append(None)
            else:
                lobj.append(None)
                lobj2.append(None)

        #ICI, ATTRIBUTION DES OBJETS AUX TRAJECTOIRES
        #POUR CHAQUE OBJET :
        # - ATTRIBUER A LA TRAJECTOIRE LA PLUS LONGUE
        for iob in range(len(lobj)):
            if lobj[iob] is not None:
                samobj = [iob] #List of tracks that match the same object
                for job in range(len(lobj)):
                    if lobj[iob] == lobj[job] and not iob==job:
                        samobj.append(job)
                maxltraj = -1
                optobj = -1
                for job in samobj:
                    if Tools.comp_difftime(ltraj0[job].traj[-1].time, instm1)==0 and ltraj0[job].nobj > maxltraj:
                        maxltraj = ltraj0[job].nobj
                        optobj = job
                #Add obj to track
                if optobj > -1:
                    obj = lobj[optobj]
                    obj2 = lobj2[optobj]
                    obj2.traps["olon"]=obj.lonc
                    obj2.traps["olat"]=obj.latc
                    u_steer, v_steer = Tools.comp_steering(obj2,steering_levels,uvmean_box,dict_fld[str(it)],pos="o")
                    u_speed, v_speed = obj2.comp_mvt(ltraj0[optobj].traj[-1],pos="o")
                    obj2.traps["u_steer"] = u_steer
                    obj2.traps["v_steer"] = v_steer
                    obj2.traps["u_speed"] = u_speed
                    obj2.traps["v_speed"] = v_speed
                    obj2.traps[trackpar] = obj.traps[trackpar]
                    ltraj0[optobj].add_obj(obj2)
                    #Delete other objects
                    for job in samobj:
                        lobj[job] = None

    #Finalisation
    stt=(signtrack==1)
    for traj in ltraj0:
        bmax=False
        maxv=0.0

        if traj.nobj>0:
            blen=False
            bmax=False
            maxv=0.0
            for obj in traj.traj:
                #condition on min/max value
                if (stt and obj.traps[trackpar]>=thr_track) or ((not stt) and obj.traps[trackpar]<=thr_track):
                    maxv=obj.traps[trackpar]
                    bmax=True

        if traj.tlen(unit="h")>=mintlen: #Condition on trajectory length
            blen=True

        if bmax and blen:
            ltraj.append(traj)
        else:
            print("This trajectory does not reach the tracking threshold or the length - We skip it")

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

