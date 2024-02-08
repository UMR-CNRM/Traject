#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Traject/Ayrault algorithm

https://github.com/UMR-CNRM/Traject
Traject is disributed under Cecill-C license
(For rights and obligations, see LICENSE.txt)

UMR-CNRM, Meteo-France & CNRS
"""
# =============================================================================
#                           ALGORITHME DE SUIVI
# =============================================================================

''' AYRAULT algorithm
'''

#------------------------------------------------------------------------------
                                # Importations
#------------------------------------------------------------------------------

import Inputs,Tools
import cartopy.crs as ccrs #DEBUG
import cartopy.feature as cfeature #DEBUG
import matplotlib.pyplot as plt #DEBUG
from datetime import timedelta, datetime
import copy
from traject import *

DEBUGGING=False

def track(algo,indf,linst,lfile,**kwargs):

    #Initialisation of domtraj, diag_parameter, Hn, deltat depending on the input parameters
    domtraj, res, deltat, Hn, diag_parameter, lparam, subnproc = Tools.init_algo(algo,indf,linst,lfile,**kwargs)

    #Tracking variables used by the algorithm
    parfilt=algo.parfilt #Effective resolution of the fields
    radmax = algo.varalgo["radmax"] #Radius (km) inside which only one minimum is kept (the highest vorticity value) 
    distmax = algo.varalgo["speedmax"]*deltat #max distance (km) between 2 kernels at successive instants
    niter=algo.varalgo["niter"] #Number of iterations
    pdmin=algo.varalgo["pdmin"] #Quality of advection
    AVOr=algo.varalgo["AVOr"] #Absolute reference value used for quality of advection (unit : same as trackpar)
    pvmin=algo.varalgo["pvmin"] #Quality of vorticity
    steering_levels = algo.varalgo["steering_levels"]
             #hPa levels taken into account for the computation of steering flow
    #uvmean_box
    if "uvmean_box" in algo.varalgo:
        uvmean_box = algo.varalgo["uvmean_box"]
    else:
        uvmean_box = 0

    if "mintlen" in algo.varalgo:
        mintlen=algo.varalgo["mintlen"] #Minimum length required (in hours) to keep the track at the end
    else:
        mintlen = 0

    accbool = algo.varalgo["accbool"] #if the acceleration criterion is applied or not
    if accbool:
        accfct = algo.varalgo["accfct"] #factor applied to the acceleration criterion (the higher accfct, the higher tolerance; accfct=infty is equivalent to accbool=False)
    else:
        accfct=0.0

    if "basetime" in kwargs:
        basetime=kwargs["basetime"]
    else:
        basetime=datetime.strftime(linst[0],format=time_fmt)

    filtapply=Tools.check_filtapply(algo,indf,linst,lfile,res)
             # 1 if filtering must be applied in the routine ; 0 if the input data has already been filtered

    trackpar, signtrack = Tools.get_parsign(algo.specfields["track"], Hn)
                     #parameter whose extremum is tracked and associated sign (-1=min ; 1=max)
                     #for cyclone : rv850
    thr_track = algo.varalgo["thr_track"] #Minimum value of track parameter (rv850 for instance) to keep the track at the end
    reseff = parfilt[trackpar] #Effective resolution (kilometer) for the detection of vortices (used in the iterations)
    track_parameter = [trackpar,"u_steer1","v_steer1","u_steer2","v_steer2","Zqual","Isuiv"]

    if "thr_core" in algo.varalgo:
        thr_core = algo.varalgo["thr_core"] #Minimum value of track parameter (btir for instance) to detect the cores to be tracked
    else:
        thr_core = thr_track

    if "pairing" in algo.specfields:
        pairpar, signpair = Tools.get_parsign(algo.specfields["pairing"], Hn)
                     #parameter whose extremum is paired to and associated sign (-1=min ; 1=max)
                     #for cyclone : mslp
        track_parameter.append(pairpar)
        thr_pairing = algo.varalgo["thr_pairing"] #threshold value of the pairing parameter (mslp for instance) to keep the minimum core
    else:
        pairpar=""
        signpair=0
        thr_pairing = 0.0

    lev1=int(steering_levels[0])
    lev2=int(steering_levels[1])
        #u_steer1, v_steer1: associated to the 1st level in steering_levels
        #u_steer2, v_steer2: associated to the 2nd level in steering_levels
    
    ss=0.0
    rd=0.0
    if "ss" in algo.varalgo:
        ss=algo.varalgo["ss"] #Size of the square for search of extremas to compute diagnostics (in degrees)
    if "rd" in algo.varalgo:
        rd=algo.varalgo["rd"] #Radius of the circle for search of extremas to compute diagnostics (in km)

    #Initialisation of variables to start the loops
    ltraj=[] #Output list of trajectories
    ltraj0=[] #Output list of trajectories
    dico_ker={} #Dictionary that contains the object candidates at the different instants
    zquald1=0.001
    zquald2=0.001

#=========================================================================#
    ''' Step 1: extraction of candidate points at every time step (fill in dico_ker) '''
#=========================================================================#    

    it = 0
    dict_fld = {}
    while it<len(linst) and Inputs.check_file(lfile[it],indf,trackpar,trackpar):

        lobj0, dict_fld0 = init_cores(it,lparam,lfile[it],linst[it],indf,algo,domtraj,res,basetime,subnproc,parfilt,filtapply,trackpar,signtrack,track_parameter,radmax,thr_core, uvmean_box, lev1, lev2)
        dict_fld[str(it)] = dict_fld0
        dico_ker[str(it)] = lobj0
        it = it + 1

#=========================================================================#
    ''' Step 2: building of trajectories by pairing successive objects'''
    #Isuiv is the track number for the different objects - it is initialized at zero values in Step 1
    #Zqual is the total quality factor - it is initialized in Step 1 (at zero value)
#=========================================================================#

    nsuiv = 0 #counter on track numbers

    #Initialization of Isuiv at time it=0
    for obj in dico_ker['0']:
        obj.traps["Zqual"] = 0.0
        if obj.traps["Isuiv"]==0:
            nsuiv = nsuiv + 1
            obj.traps["Isuiv"]=nsuiv

    #MAIN LOOP
    for jter in range(niter):

        #Number of successive instants used for matching test
        if jter==0:
            inbinst=2
        else:
            inbinst=4
        itdeb=int(inbinst/2)-1
        itfin=len(linst)-itdeb-1

        print("------------------------------")
        print("------- ITERATION "+str(jter+1)+" ----------")
        print("from " + datetime.strftime(linst[itdeb],"%Y%m%d%H") + " to " +  datetime.strftime(linst[itfin],"%Y%m%d%H"))
        print("------------------------------")
                
        it=itdeb
        while it<itfin and Inputs.check_file(lfile[it],indf,trackpar,trackpar): #Loop on instants from instant 2 
            inst=linst[it]
            print(inst,lfile[it])

            #loop on the kernels at it+1
            for obj2 in dico_ker[str(it+1)]:
                obj2.traps["Zqual"]=0.0 

                #loop on the kernels at it
                for obj1 in dico_ker[str(it)]:
                    if (obj1.traps["Isuiv"]>0) and (Tools.comp_length(obj1.lonc,obj1.latc,obj2.lonc,obj2.latc)<distmax):
                    #To save computation time, the obj1 too far away from obj2 or the non-tracked obj1 are not processed
                        if DEBUGGING:
                            print('-')
                            print("From Instant t=", inst, "(it="+str(it)+")", " to Instant t+1=",linst[it+1],"(it+1="+str(it+1)+")")
                            print("--Iteration "+str(jter+1)+"-- Processing object1 no "+str(obj1.traps["Isuiv"])+"/("+str(obj1.lonc)+","+str(obj1.latc),")"+\
                                " to object2 no "+str(obj2.traps["Isuiv"])+"/("+str(obj2.lonc)+","+str(obj2.latc),")")
                        zquald = Tools.adv_quality(trackpar,AVOr,obj1, obj2,res[trackpar],pdmin) #advection quality
                        zqualv = Tools.avo_quality(trackpar,obj1,obj2,pvmin) #vorticity quality

                        #Determining adjacent objects (after 1st iteration)
                        if jter>0 and accbool: 
                            b0=False
                            b3=False
                            for obj in dico_ker[str(it-1)]:
                                if obj.traps["Isuiv"]==obj1.traps["Isuiv"]:
                                    b0=True
                                    obj0=obj
                            for obj in dico_ker[str(it+2)]:
                                if obj.traps["Isuiv"]==obj2.traps["Isuiv"]:
                                    b3=True
                                    obj3=obj
                            #Compute pseudo-position quality
                            zquald1=0.001
                            zquald2=0.001
                            if b3:
                                zquald1 =  Tools.pseudo_pos_quality(obj2,obj3,obj1,res[trackpar]*accfct,pdmin)
                            if b0:
                                zquald2 =  Tools.pseudo_pos_quality(obj1,obj0,obj2,res[trackpar]*accfct,pdmin)
                            if not b3:
                                zquald1=zquald2
                            if not b0:
                                zquald2=zquald1

                        if jter==0 or not accbool:
                            zqualtot = zquald * zqualv
                        else:
                            zqualtot = zquald * zqualv * zquald1 * zquald2

                        if (zqualtot > obj2.traps["Zqual"]): #On apparie obj2 a obj1
                            if DEBUGGING:
                                print("Total quality:",zqualtot,zquald1*zquald2,)
                                print("Pairing (",obj1.lonc,obj1.latc,") at t with ((",obj2.lonc,obj2.latc,") at t+dt - Isuiv=",obj1.traps["Isuiv"],obj2.traps["Isuiv"])
                                print("With wind "+str(steering_levels[0])+"hPa:", obj1.traps["u_steer1"], obj1.traps["v_steer1"], " - "+\
                                    str(steering_levels[1])+"hPa:" , obj1.traps["u_steer2"], obj1.traps["v_steer2"])
                            obj2.traps["Isuiv"] = obj1.traps["Isuiv"]
                            obj2.traps["Zqual"] = zqualtot

                if obj2.traps["Zqual"]==0.0:
                    obj2.traps["Isuiv"] = 0

            #Choix de la meilleure qualite en amont
            for obj1 in dico_ker[str(it)]:
                zq=0.0
                objq=obj1
                for obj2 in dico_ker[str(it+1)]:
                    if (obj2.traps["Isuiv"] == obj1.traps["Isuiv"]) and (obj2.traps["Zqual"]>zq):
                        zq = obj2.traps["Zqual"]
                        objq = obj2
                for obj2 in dico_ker[str(it+1)]:
                    if (obj2.traps["Isuiv"] == obj1.traps["Isuiv"]) and not (obj2==objq):
                        obj2.traps["Isuiv"] = 0
                        obj2.traps["Zqual"] = 0.0

            #Numbering of all unpaired kernels at time it+1
            for obj2 in dico_ker[str(it+1)]:
                if obj2.traps["Isuiv"]==0:
                    nsuiv = nsuiv + 1
                    obj2.traps["Isuiv"]=nsuiv

            it = it + 1

#=========================================================================#
    ''' Step 3: Selection and output trajectories '''
#=========================================================================#

    #Loop on nsuiv
    st=(signtrack==1)
    inam=0
    for isuiv in range(1,nsuiv+1):
        traj=DefTrack(algo.classobj,basetime=basetime)
        for it in range(len(dico_ker)):
            for obj in dico_ker[str(it)]:
                if obj.traps["Isuiv"]==isuiv:
                    traj.add_obj(obj)
        if traj.nobj>0:
            blen=False
            bmax=False
            maxv=0.0
            for obj in traj.traj:
                #condition on min/max value
                if (st and obj.traps[trackpar]>=thr_track) or ((not st) and obj.traps[trackpar]<=thr_track):
                    maxv=obj.traps[trackpar]
                    bmax=True

            if traj.tlen(unit="h")>=mintlen: #Condition on trajectory length
                blen=True

            print('-------------------')
            print("Output: Examining traj no "+str(isuiv))
            print("Max value along the track: "+str(maxv))
            print("Duration (h) of the track: "+str(traj.tlen(unit="h")))
            if blen and bmax:
                ltraj0.append(traj)

    for traj in ltraj0:
        traj2=DefTrack(algo.classobj,basetime=basetime)
        inam=inam+1
        #Give name to traj
        traj2.name="AY-"+basetime+ '-' +str(inam)
        ltraj.append(traj2)

    for it in range(len(linst)):
        for ir in range(len(ltraj0)):
            #Correction of position if a mslp minimum is found - Computation of diagnostics
            for io in range(ltraj0[ir].nobj):
                obj=ltraj0[ir].traj[io]
                traj2=ltraj[ir]
                if Tools.comp_difftime(ltraj0[ir].traj[io].time,linst[it])==0:
                    if not pairpar == "":
                        obj_pair = obj.search_core(dict_fld[str(it)],linst[it],pairpar,Hn,track_parameter,diag_parameter,ss=ss,thr_param=thr_pairing,pairing=True,smooth=True)
                        if obj_pair is not None:
                            obj_pair.traps[trackpar] = obj.traps[trackpar]
                            traj2.add_obj(copy.deepcopy(obj_pair))
                        else:
                            Tools.make_diags(diag_parameter,obj,dict_fld[str(it)],obj.lonc,obj.latc) #Add diagnostics in traj
                            traj2.add_obj(copy.deepcopy(obj))
                    else:
                        Tools.make_diags(diag_parameter,obj,dict_fld[str(it)],obj.lonc,obj.latc) #Add diagnostics in traj
                        traj2.add_obj(copy.deepcopy(obj))

    return ltraj

def init_cores(it,lparam,fic,inst,indf,algo,domtraj,res,basetime,subnproc,parfilt,filtapply,trackpar,signtrack,track_parameter,radmax,thr_core,uvmean_box,lev1,lev2):

    dict_fld0 = Tools.load_fld(lparam,fic,inst,indf,algo,domtraj,res,basetime,subnproc,parfilt=parfilt,filtapply=filtapply) #input fields at the given instant
    lobj0 = Search_allcores(algo.classobj,dict_fld=dict_fld0,inst=inst,trackpar=trackpar,signtrack=signtrack,track_parameter=track_parameter,dist=radmax, thr=thr_core)

    #Computation of steering flow parameters
    u_steer1, v_steer1 = Tools.comp_steering(lobj0,[lev1],uvmean_box,dict_fld0)
    u_steer2, v_steer2 = Tools.comp_steering(lobj0,[lev2],uvmean_box,dict_fld0)

    for ivi in range(len(lobj0)):
        lobj0[ivi].traps["u_steer1"] = u_steer1[ivi]
        lobj0[ivi].traps["v_steer1"] = v_steer1[ivi]
        lobj0[ivi].traps["u_steer2"] = u_steer2[ivi]
        lobj0[ivi].traps["v_steer2"] = v_steer2[ivi]
        lobj0[ivi].traps["Zqual"] = 0.0
        lobj0[ivi].traps["Isuiv"] = 0

    return lobj0, dict_fld0
