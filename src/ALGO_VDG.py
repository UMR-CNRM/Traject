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

        #print("reftraj")
        #print(reftraj.__dict__)

        #------------------------------------------------------------------------------
                        # First step: Start the track from the closest point to reftraj
        #------------------------------------------------------------------------------

        it0=-1
        gook=False
        while not gook and it0<len(linst2)-1:
            it0=it0+1
            print(linst2[it0])
            refobj, gook1 =reftraj.find_inst(linst2[it0]) #Attention, delta t différent !!
            if gook1:
                lon,lat,val,valpair,olon,olat,gook = search_pairing(trackpar,signtrack,thr_track,pairpar,signpair,\
                thr_pair,basetime,parfilt,filtapply,refobj[0].lonc,refobj[0].latc,lfile[it0],linst2[it0],indf,domtraj,res,ss,subnproc)

        if gook: #A 
            #Creation of the track
            print("A starting point has been found at time : ",linst2[it0])
            print("Reference point:",refobj[0].lonc,refobj[0].latc)
            print("Point found:",lon,lat)
            traj = DefTrack(algo.classobj,basetime=basetime)
            traj.name = reftraj.name
            objectm0 = DefObject(algo.classobj, [], track_parameter,lonc=lon,latc=lat,time=linst2[it0])
            objectm0.traps["olon"]=olon
            objectm0.traps["olat"]=olat
            u_steer, v_steer = Tools.comp_steering(objectm0.traps["olon"],objectm0.traps["olat"],steering_levels,uvmean_box,\
                    lfile[it0],linst2[it0],indf,res,domtraj,basetime,parfilt,filtapply,subnproc)

            objectm0.traps["u_steer"] = u_steer
            objectm0.traps["v_steer"] = v_steer
            objectm0.traps["u_speed"] = u_steer #we cannot determine a movement, so we take wind_steer
            objectm0.traps["v_speed"] = v_steer
            objectm0.traps[trackpar] = val
            objectm0.traps[pairpar] = valpair
            #if pairpar in diag_parameter: #Add pairing diagnostic in traj #Depreciated
            #    obj.diags.append(pairpar)
            #    setattr(obj,pairpar,valpair)

            Tools.make_diags(diag_parameter,objectm0,lfile[it0],linst2[it0],indf,domtraj,Hn,res,basetime,olon,olat,subnproc,parfilt=parfilt,filtapply=filtapply) #Add diagnostics in traj
            traj.add_obj(objectm0)

        #------------------------------------------------------------------------------
                                    # TIME LOOP
        #------------------------------------------------------------------------------
                                    
        it=it0

        while gook and it<len(linst)-1 and Inputs.check_file(lfile[it+1],indf,trackpar,trackpar):
            it = it + 1
            print("Looking for a point at instant ", linst[it])
            instm1=linst[it-1]
            inst=linst[it]
            print(inst,lfile[it])
            objm,fnd=traj.find_inst(instm1) #Object at the previous instant
            objm1 = objm[0]

            #Computation of a guess (advection by the steering flow and displacement speed at objm1)
            #We advect olon, olat
            obj_guess = objm1.advect(inst, (1-w)*objm1.traps["u_steer"]+ w*objm1.traps["u_speed"],
                (1-w)*objm1.traps["v_steer"]+ w*objm1.traps["v_speed"],pos="o")
            #print("STEERING WIND:",objm1.traps["u_steer"],objm1.traps["v_steer"])
            #print("SPEED WIND:",objm1.traps["u_speed"],objm1.traps["v_speed"])
            #print(objm1.__dict__)
            #print(obj_guess.__dict__)

            lon,lat,val,valpair,olon,olat,gook = search_pairing(trackpar,signtrack,thr_track,pairpar,signpair,\
                thr_pair,basetime,parfilt,filtapply,obj_guess.lonc,obj_guess.latc,lfile[it],linst[it],indf,domtraj,res,ss,subnproc)

            #Some plot to debug
                  #DEBUG - Plot champs et des points candidats (pour voir)
            if DEBUGGING:

                flda = Inputs.extract_data(lfile[it],linst[it],indf,trackpar,domtraj,res[trackpar],basetime,subnproc,filtrad=parfilt[trackpar]*filtapply)
                Tools.plot_F_LL(flda,[objm1.traps["olon"]],[objm1.traps["olat"]],trackpar+datetime.strftime(linst[it],format=time_fmt),nl=20,\
                    vect=[(1-w)*objm1.traps["u_steer"]+ w*objm1.traps["u_speed"],(1-w)*objm1.traps["v_steer"]+ w*objm1.traps["v_speed"]])
                fldb = Inputs.extract_data(lfile[it],linst[it],indf,pairpar,domtraj,res[trackpar],basetime,subnproc,filtrad=parfilt[trackpar]*filtapply)
                Tools.plot_F_LL(fldb,[objm1.lonc],[objm1.latc],pairpar+datetime.strftime(linst[it],format=time_fmt),nl=20)

            if gook:
                #CREATE POINT ON THE TRACK
                objectm = DefObject(algo.classobj, [], track_parameter,lonc=lon,latc=lat,time=inst)
                objectm.traps["olon"]=olon
                objectm.traps["olat"]=olat
                u_steer, v_steer = Tools.comp_steering(objectm.traps['olon'],objectm.traps["olat"],steering_levels,uvmean_box,\
                        lfile[it],linst[it],indf,res,domtraj,basetime,parfilt,filtapply,subnproc)
                objectm.traps["u_steer"] = u_steer
                objectm.traps["v_steer"] = v_steer
                u_speed, v_speed = objectm.comp_mvt(objm1,pos="o")
                objectm.traps["u_speed"] = u_speed
                objectm.traps["v_speed"] = v_speed
                objectm.traps[trackpar] = val
                objectm.traps[pairpar] = valpair

                # WRITE DIAGNOSTIC PARAMETERS
                #if pairpar in diag_parameter: #Add pairing diagnostic in traj #Depreciated
                #    obj.diags.append(pairpar)
                #    setattr(obj,pairpar,valpair)
                Tools.make_diags(diag_parameter,objectm,lfile[it],linst[it],indf,domtraj,Hn,res,basetime,olon,olat,subnproc,parfilt=parfilt,filtapply=filtapply) #Add diagnostics in traj
                traj.add_obj(objectm)

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

def search_core(trackpar,signtrack,thr_track,pairpar,signpair,thr_pair,basetime,parfilt,filtapply,g_lon,g_lat,filin,inst,indf,domtraj,res,ss,subnproc):
    #Search the core from the guess position g_lon, g_lat obtained after advection from the previous instant
    #trackpar is the parameter which is tracked from one instant to the other (for instance : rv850)
    #pairpar is the core of the structure (for instance: mslp)
    #ss is the square (in degrees) inside which the core is searched from the advected structure

    lon=missval
    lat=missval

    #Extract trackpar field
    flda = Inputs.extract_data(filin,inst,indf,trackpar,domtraj,res[trackpar],basetime,subnproc,filtrad=parfilt[trackpar]*filtapply)
    lon0,lat0,val,dist,gook1 = Tools.find_locextr(signtrack, flda,g_lon,g_lat,ss=ss,thr=thr_track)
    #print(trackpar+" kernel: ", val, "at location :", lon0, lat0)

    if gook1:
        fldb = Inputs.extract_data(filin,inst,indf,pairpar,domtraj,res[pairpar],basetime,subnproc,filtrad=parfilt[pairpar]*filtapply)
        lon,lat,val,dist,gook2 = Tools.find_locextr(signpair, fldb,lon0,lat0,ss=ss,thr=thr_pair)
        #print(pairpar+" kernel: ", val, "at location :", lon, lat)
        #tests
        #tlon,tlat,tval=Tools.find_allmin(fldb,thr=100000.00)
        #print("guess:",lon0,lat0,'+/-',ss/2)
        #print("All mslp min:")
        #for ivi in range(len(tval)):
        #    print(tlon[ivi],tlat[ivi],tval[ivi])
        #print(gook2)
    else:
        gook2 = False

    if gook2:
        #LISSAGE DU MINIMA (AVEC  PROCHES VOISIN DE lon,lat)
        lon, lat, val = Tools.smooth_extr(signpair,fldb,lon,lat)
    else:
        val=missval

    return lon,lat,val,lon0,lat0,gook2

def search_pairing(trackpar,signtrack,thr_track,pairpar,signpair,thr_pair,basetime,parfilt,filtapply,g_lon,g_lat,filin,inst,indf,domtraj,res,ss,subnproc):
    #Search the core from the guess position g_lon, g_lat obtained after advection from the previous instant
    #trackpar is the parameter which is tracked from one instant to the other (for instance : rv850)
    #pairpar is the core of the structure (for instance: mslp)
    #ss is the square (in degrees) inside which the core is searched from the advected structure
    #Compared to search_core, this is a soft pairing: if no pairing core is found, we use the position of the tracking core

    lon=missval
    lat=missval
    val=missval
    if signpair==1:
        valpair=missval
    else:
        valpair=-missval

    #Extract trackpar field
    flda = Inputs.extract_data(filin,inst,indf,trackpar,domtraj,res[trackpar],basetime,subnproc,filtrad=parfilt[trackpar]*filtapply)
    #print("kernel guess :", g_lon, g_lat)
    lon0,lat0,val,dist,gook1 = Tools.find_locextr(signtrack, flda,g_lon,g_lat,ss=ss,thr=thr_track)
    #print(trackpar+" kernel: ", val, "at location :", lon0, lat0)

    if gook1:
        fldb = Inputs.extract_data(filin,inst,indf,pairpar,domtraj,res[pairpar],basetime,subnproc,filtrad=parfilt[pairpar]*filtapply)
        lon,lat,valpair,dist,gook2 = Tools.find_locextr(signpair, fldb,lon0,lat0,ss=ss,thr=thr_pair)
        #print(pairpar+" kernel: ", valpair, "at location :", lon, lat)
        #tests
        #tlon,tlat,tval=Tools.find_allmin(fldb,thr=100000.00)
        #print("guess:",lon0,lat0,'+/-',ss/2)
        #print("All mslp min:")
        #for ivi in range(len(tval)):
        #    print(tlon[ivi],tlat[ivi],tval[ivi])
        #print(gook2)
    else:
        gook2 = False

    if gook2:
        #LISSAGE DU MINIMA (AVEC  PROCHES VOISIN DE lon,lat)
        lon, lat, valpair = Tools.smooth_extr(signpair,fldb,lon,lat)
    else:
        lon = lon0
        lat = lat0
        #Pas utile d'ajouter un lissage sur trackpar

    return lon,lat,val,valpair,lon0,lat0,gook2
