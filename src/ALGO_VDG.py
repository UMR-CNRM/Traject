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

DEBUGGING=False

def track(algo,indf,linst,lfile,**kwargs):

    #Initialisation of domtraj, res, deltat depending on the input parameters
    domtraj, res, deltat, Hn, diag_parameter, subnproc = Tools.init_algo(algo,indf,linst,lfile,**kwargs)

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

    thr_track = algo.varalgo["thr_track"] #Minimum absolute value of tracked core
    thr_pair = algo.varalgo["thr_pairing"] #Minimum value of paired core

    track_parameter = [trackpar,pairpar,"olon","olat","u_speed","v_speed","u_steer","v_steer"] #olon, olat : coordinates of the trackpar extrema, used for tracking

    rd=0.0
    if "rd" in algo.varalgo:
        rd=algo.varalgo["rd"] #Radius of the circle for search of extremas to compute diagnostics (in km)

    ltraj=[]

    if "basetime" in kwargs:
        basetime=kwargs["basetime"]
    else:
        basetime=datetime.strftime(linst[0],format=time_fmt)

    if "last_start" in algo.varalgo: #List of fc instants where to search for starting point
        linst2 = [inst for inst in linst if Tools.comp_difftime(linst[0],inst)<=algo.varalgo["last_start"]]
    else:
        linst2 = linst

    #Read reference track and find the list of tracks valid at any instant in linst2 and in the domain
    if "reftraj" in kwargs:
        lref = Tools.get_reftraj(linst2,domtraj,kwargs["reftraj"])
    else:
        lref=[]
    if len(lref)==0:
        print("WARNING - No reference track found could be found (required for VDG tracking algorithm)")

    #Loops on points found in the reference tracks
    for iref in range(len(lref)):
        if "traj" in locals():
            del traj
        reftraj=lref[iref]

        #------------------------------------------------------------------------------
                        # First step: Start the track from the closest point to reftraj
        #------------------------------------------------------------------------------

        it0=-1
        gook=False
        obj2=None
        exclude=False
        while not gook and it0<len(linst2)-1:
            it0=it0+1
            print(it0)
            refobj, gook1 =reftraj.find_inst(linst2[it0]) #Attention, delta t différent !!
            if gook1: # Tracking main core
                obj = refobj[0].search_core(lfile[it0],linst[it0],trackpar,Hn,indf,algo,domtraj,res,parfilt,filtapply,track_parameter,basetime,subnproc,[],ss=ss,thr_param=thr_track)
                obj2 = None
                if obj is not None: # Pairing
                    obj2 = obj.search_core(lfile[it0],linst[it0],pairpar,Hn,indf,algo,domtraj,res,parfilt,filtapply,track_parameter,basetime,subnproc,diag_parameter,ss=ss,thr_param=thr_pair,pairing=True,smooth=True)
                    isok, exclude = obj2.conditiontype(lfile[it0],linst[it0],indf,algo,domtraj,res,parfilt,filtapply,basetime,subnproc,init=True)
            gook = obj is not None and obj2 is not None and not exclude

        if gook: #A 
            #Creation of the track
            print("A starting point has been found at time : ",linst2[it0])
            print("Reference point:",refobj[0].lonc,refobj[0].latc)
            print("Point found:",obj2.lonc,obj2.latc)
            traj = DefTrack(algo.classobj,basetime=basetime)
            traj.name = reftraj.name

            #Initialisations of obj2 variables and addition to traj
            obj2.traps["olon"]=obj.lonc
            obj2.traps["olat"]=obj.latc
            u_steer, v_steer = Tools.comp_steering(obj2,steering_levels,uvmean_box,lfile[it0],linst2[it0],indf,res,domtraj,basetime,parfilt,filtapply,subnproc,pos="o")

            obj2.traps["u_steer"] = u_steer
            obj2.traps["v_steer"] = v_steer
            obj2.traps["u_speed"] = u_steer #we cannot determine a movement at step 0, so we take wind_steer
            obj2.traps["v_speed"] = v_steer
            obj2.traps[trackpar] = obj.traps[trackpar]

            print(obj2.nameobj)
            traj.add_obj(obj2)

        #------------------------------------------------------------------------------
                                    # TIME LOOP
        #------------------------------------------------------------------------------
        it=it0
        gook= obj2 is not None and not exclude

        while gook and it<len(linst)-1 and Inputs.check_file(lfile[it+1],indf,trackpar,trackpar):
            it = it + 1
            print("Looking for a point at instant ", linst[it])
            instm1=linst[it-1]
            inst=linst[it]
            objm,fnd=traj.find_inst(instm1) #Object at the previous instant
            objm1=objm[0]

            #Computation of a guess (advection by the steering flow and displacement speed at objm1)
            #We advect olon, olat
            obj_guess = objm1.advect(inst, (1-w)*objm1.traps["u_steer"]+ w*objm1.traps["u_speed"],
                (1-w)*objm1.traps["v_steer"]+ w*objm1.traps["v_speed"],pos="o")

            obj = obj_guess.search_core(lfile[it],linst[it],trackpar,Hn,indf,algo,domtraj,res,parfilt,filtapply,track_parameter,basetime,subnproc,[],ss=ss,thr_param=thr_track)
            obj2 = None
            if obj is not None: # Pairing
                obj2 = obj.search_core(lfile[it],linst[it],pairpar,Hn,indf,algo,domtraj,res,parfilt,filtapply,track_parameter,basetime,subnproc,diag_parameter,ss=ss,thr_param=thr_pair,pairing=True,smooth=True)
                isok, exclude = obj2.conditiontype(lfile[it],linst[it],indf,algo,domtraj,res,parfilt,filtapply,basetime,subnproc,init=False)
            gook = obj is not None and obj2 is not None and not exclude

            if gook:
                #Initialisations of obj2 variables and addition to traj
                obj2.traps["olon"]=obj.lonc
                obj2.traps["olat"]=obj.latc
                u_steer, v_steer = Tools.comp_steering(obj2,steering_levels,uvmean_box,lfile[it],linst[it],indf,res,domtraj,basetime,parfilt,filtapply,subnproc,pos="o")
                u_speed, v_speed = obj2.comp_mvt(objm1,pos="o")
                obj2.traps["u_steer"] = u_steer
                obj2.traps["v_steer"] = v_steer
                obj2.traps["u_speed"] = u_speed
                obj2.traps["v_speed"] = v_speed
                obj2.traps[trackpar] = obj.traps[trackpar]

                traj.add_obj(obj2)

            #DEBUG - Plot champs et des points candidats (pour voir)
            if DEBUGGING:

                flda = Inputs.extract_data(lfile[it],linst[it],indf,trackpar,domtraj,res[trackpar],basetime,subnproc,filtrad=parfilt[trackpar]*filtapply)
                Tools.plot_F_LL(flda,[objm1.traps["olon"]],[objm1.traps["olat"]],trackpar+datetime.strftime(linst[it],format=time_fmt),nl=20,\
                    vect=[(1-w)*objm1.traps["u_steer"]+ w*objm1.traps["u_speed"],(1-w)*objm1.traps["v_steer"]+ w*objm1.traps["v_speed"]])
                fldb = Inputs.extract_data(lfile[it],linst[it],indf,pairpar,domtraj,res[trackpar],basetime,subnproc,filtrad=parfilt[trackpar]*filtapply)
                Tools.plot_F_LL(fldb,[objm1.lonc],[objm1.latc],pairpar+datetime.strftime(linst[it],format=time_fmt),nl=20)

        if "traj" in locals():
            stt=(signtrack==1)
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

