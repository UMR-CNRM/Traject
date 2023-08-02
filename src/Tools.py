#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Traject/Utilities and Tools

https://github.com/UMR-CNRM/Traject
Traject is disributed under Cecill-C license
(For rights and obligations, see LICENSE.txt)

UMR-CNRM, Meteo-France & CNRS

"""
#=========================================================================#
#            GENERIC ROUTINES FOR THE TRACKING ALGORITHMS                 #
#=========================================================================#

#--------------------------------------------------------------------------#
#                            Importations                                  #
#--------------------------------------------------------------------------#
import numpy as np
from datetime import datetime
import cartopy.crs as ccrs
import matplotlib as mpl
import matplotlib.pyplot as plt
from datetime import timedelta, datetime
import Inputs, traject
#mpl.use('Agg')
import pyproj
from shapely.geometry import Point, Polygon
import copy

#--------------------------------------------------------------------------#
#                              Variables                                   #
#--------------------------------------------------------------------------#
r = 6371e3  # Radius of Earth (m)
omega = 7.2921 * 10**(-5)  # Angular rotation speed of Earth
maxval = 1e12 #(maximum value, infinity)

#--------------------------------------------------------------------------#
#                         General functions                                #
#--------------------------------------------------------------------------#
def init_algo(algo,indf,linst,lfile,**kwargs):
    #Initialisation of input variables needed for all algorithms

    #Time step (hours)
    deltat = abs(linst[1]-linst[0]).seconds/3600.0

    #Check and complete domtraj
    if algo.domtraj is None:
        varname=list(algo.parfilt.keys())
        domtraj, res, lons, lats = Inputs.extract_domain(lfile[0],linst[0],indf,varname[0]) #domtraj is the total grid
    else:
        domtraj=algo.domtraj
    
    #subnproc
    subnproc=1
    if algo.parallel is not None:
        if "subnproc" in algo.parallel:
            subnproc=algo.parallel["subnproc"]

    #Compute dictionary of resolutions
    dres = get_parres(algo,indf,lfile[0],linst[0])

    #Northern Hemisphere or not ?
    if domtraj["latmin"]>0.0 and domtraj["latmax"]>0.0:
        Hn = True
    elif domtraj["latmin"]<0.0 and domtraj["latmax"]<0.0:
        Hn = False
    else:
        if abs(domtraj["latmin"])>abs(domtraj["latmax"]):
            Hn=False
        else:
            Hn=True
        print("The Equator crosses the domain...")
        print("There may be some problems if detection is based on vorticity...")
        print("We continue anyway ...")

    if "pairing" in algo.specfields:
        pairpar, signpair = get_parsign(algo.specfields["pairing"], Hn)
    else:
        pairpar=""

    ldiag=[]
    if "ss" in algo.varalgo:
        ss=algo.varalgo["ss"]
    else:
        ss=traject.missval
    if "rd" in algo.varalgo:
        rd=algo.varalgo["rd"]
    else:
        rd=traject.missval
    for diagname in algo.diag_parameter:
        if diagname==pairpar and not pairpar=="": #Special case to diagnose pairpar
            diag=diagdef(pairpar,Hn,ss,rd)
            ldiag.append(diag)
        else:
            diag=diagdef(diagname,Hn,ss,rd)
            ldiag.append(diag)

    return domtraj, dres, deltat, Hn, ldiag, subnproc

#--------------------------------------------------------------------------#

def get_parres(algo,indf,file1,inst1):
    #Defines the grid resolution of parameter fields, depending on algo and indf
    #Output : a dictionary giving the resolution of every required parameters

    dres={}
    lpar=[k for k in algo.parfilt] #list of required variables
    reskeys=[k for k in algo.parres]

    if reskeys[0]=="all": #then all fields apply the same rule
        if algo.parres["all"]<1e-6:
            #then all fields keep their original resolution
            for par in lpar:
                domtraj, res, lons, lats = Inputs.extract_domain(file1,inst1,indf,par)
                dres[par] = res
                algo.parres[par] = res #new
        else:
            #then all fields have same original resolution given by algo.parres
            for par in lpar:
                dres[par]=algo.parres["all"]

    else: #the values should be given for all fields
        for par in lpar:
            if par in reskeys:
                if algo.parres[par]<1e-6:
                    #then the field keeps its original resolution
                    domtraj, res, lons, lats = Inputs.extract_domain(file1,inst1,indf,par)
                    dres[par]=res
                    algo.parres[par] = res #new
                else:
                    #then resolution is given by algo.parres
                    dres[par]=algo.parres[par]
            else:
                print("Error in Tools.get_parres when reading algodef.parres")
                print("--> ")
                print("algo.parres:",algo.parres)
                print("algo.parfilt:",algo.parfilt)
                print("Please declare a resolution in parres for all parameters in parfilt")
                exit()

    return dres

#--------------------------------------------------------------------------#

def get_parsign(par, Hn):
    #Input : parameter and hemisphere
    #Output : parameter name and sign to define if the minimum (-1)
        #or the maximum (1)should be tracked
    #If the input parameter name has a "_max" or "_min" suffix,
    #then the sign value is forced
    #If not, then the sign value is guessed from the parameter name and the 
    #hemisphere

    par2 = par.split("_")

    if par2[-1]=="min" or par2[-1]=="max":
    #First case: the sign (min or max) is forced in the parameter name
        parname='_'.join(par2[0:-1])
        if par2[-1]=="min":
            parsign=-1
        elif par2[-1]=="max":
            parsign=1
    else:
    #Second case: the sign in not given in the parameter name ... so we guess it here
        parname=par
        if parname=="mslp" or parname=="btir":
            parsign=-1
        elif parname[0:2]=="rr" or parname[0:2]=="ff":
            parsign=1
        elif parname[0:2]=="rv" or parname[0:2]=="av":
            if Hn:
                parsign=1
            else:
                parsign=-1
        else:
            print("ABORT in Tools.get_parsign: please add parameter name : "+par+" in this routine")
            exit()

    return parname, parsign
#--------------------------------------------------------------------------#

def check_filtapply(algo,indf,linst,lfile,dres):
    #Checks if the filter (or change of resolution) must be applied in the algorithm --> filtapply=1
    #or if it has been applied and prepared before --> filtapply =0
    #Output : filtapply (0 or 1)

    filtapply=1

    #check resolutions (between the grid and the expected value in dres)
    checkres=True
    lpar=[k for k in algo.parfilt] #list of required variables

    infilt={} #Get input filter
    if "filtered" in indf.__dict__.keys():
        for par in lpar:
            if par in indf.filtered.keys():
                infilt[par]=indf.filtered[par]
            else:
                infilt[par]=0
    else:
        for par in lpar:
            infilt[par]=0

    for par in lpar:
        domtraj, res, lons, lats = Inputs.extract_domain(lfile[0],linst[0],indf,par)
        checkres = checkres and abs(dres[par]-res)<1e-6

    #check filters (between indf and algo)
    checkfilt=True
    for par in lpar:
        if algo.parfilt[par]>1e-6:
            checkfilt = checkfilt and abs(infilt[par]-algo.parfilt[par])<1e-6

    if checkres and checkfilt:
        filtapply=0

    return filtapply
#--------------------------------------------------------------------------#
def get_reftraj(linst,domtraj,refnam):
    #Reads the reference track information (if provided)
    #Reference track is required by some algorithms
    #Output : list of tracks which have at least a common point with linst)

    lreftraj=[]
    if isinstance(refnam,str): #Filename
        lt = traject.Read(refnam)
    elif isinstance(refnam,list): #List of tracks
        lt = refnam
    else: #Single track
        lt = [refnam]
    #lt : list of tracks

    for i in range(len(lt)): #Loop on tracks in the first file
        traki=lt[i]
        chkall=False
        indom=False
        for inst in linst:
            objs,chk=traki.find_inst(inst) #find valid instants in the tracks
            if chk:
                chkall=True
                obj=objs[0]
                if obj.check_dom(domtraj):
                    indom=True
        if chkall and indom:
            lreftraj.append(traki)

    return lreftraj


#--------------------------------------------------------------------------#
###Some lat/lon matrix management
#--------------------------------------------------------------------------#
def ll_to_ij(lons,lats,lon1,lat1):
    #Converts point (lon1,lat1) coordinates to the closest (xi,yi) index
    #chk is True if the point is in the domain
    #lons, lats must be in increasing order

    lons2=[lons[ivi] for ivi in range(len(lons))]
    lon12=lon1
    if lons[0]>lons[-1]: #ordering lons if needed
        for ivi in range(len(lons)):
            if lons[ivi]>180.0:
                lons2[ivi]=lons[ivi]-360.0
            if lon1>180.0:
                lon12=lon1-360.0

    resx = lons2[1] - lons2[0]
    resy = lats[1] - lats[0]

    xi = int(round((lon12 - lons2[0])/resx))
    yi = int(round((lat1 - lats[0])/resy))

    chk = (xi>-1 and xi<len(lons)) and (yi>-1 and yi<len(lats))

    return xi, yi, chk

#--------------------------------------------------------------------------#
###Time management
#--------------------------------------------------------------------------#

def get_time_index(linst,inst0):
    #Get the index it in linst that corresponds to time inst0

    it=int(traject.missval)

    for ivi in range(len(linst)):
        if inst0==datetime.strftime(linst[ivi],Inputs.time_fmt):
            it = ivi

    return it

def get_validtime(bt,lt=0):
    #returns valid time from basetime bt and lead time lt (hours)
    #output is in string fmt

    if isinstance(bt,str):
        tt= datetime.strptime(bt,Inputs.time_fmt)+timedelta(hours=lt)
    else:
        tt= bt + timedelta(hours=lt)

    vtime = datetime.strftime(tt,Inputs.time_fmt)

    return vtime

def comp_difftime(bt1,bt2):
    #returns time difference between instant bt1 and bt2 (bt2-bt1)
    #output is integer format (hours)

    if isinstance(bt1,str):
        tt1= datetime.strptime(bt1,Inputs.time_fmt)
    else:
        tt1= bt1

    if isinstance(bt2,str):
        tt2= datetime.strptime(bt2,Inputs.time_fmt)
    else:
        tt2= bt2

    tdiff = 24*(tt2-tt1).days + ((tt2-tt1).seconds/3600.0)

    return int(tdiff)

def set_time_limits(ltraj):

    mint = datetime.strptime(ltraj[0].traj[0].time, Inputs.time_fmt)
    maxt = datetime.strptime(ltraj[0].traj[ltraj[0].nobj-1].time, Inputs.time_fmt)
    for traj in ltraj:
        min1 = datetime.strptime(traj.traj[0].time, Inputs.time_fmt)
        max1 = datetime.strptime(traj.traj[traj.nobj-1].time, Inputs.time_fmt)
        if min1 < mint:
            mint = min1
        if max1 > maxt:
            maxt = max1

    return mint, maxt

def get_dom_limits(ltraj, diag, res):
    #From a list of trajectories, outputs the rectangular limits around the list of trajectories
    #The limits can be determined from the diagnostic values or by the object itself (diag="")

    llon=[]
    llat=[]
    if diag=="":
        for traj in ltraj:
            for obj in traj.traj:
                llon.append(obj.lonc)
                llat.append(obj.latc)
    else:
        dd=guess_diag(diag,True)
        if dd.kind[0]=='c':
            for traj in ltraj:
                for obj in traj.traj:
                    ivi=0
                    ld=obj.__dict__[diag]
                    while ivi<len(ld):
                        llon.append(ld[ivi])
                        llat.append(ld[ivi+1])
                        ivi = ivi + 3
        else:
            for traj in ltraj:
                for obj in traj.traj:
                    llon.append(obj.__dict__[diag][0])
                    llat.append(obj.__dict__[diag][1])

    #Compute length of grid and grid boundaries
    dyy = comp_length(np.mean(llon),np.mean(llat),np.mean(llon),np.mean(llat)+res)
    lonmin = res * int( (min(llon)/res ) )
    lonmax = res * int( (max(llon)/res ) ) + res
    latmin = res * int( (min(llat)/res ) )
    latmax = res * int( (max(llat)/res ) ) + res

    return {"lonmin":lonmin,"lonmax":lonmax,"latmin":latmin,"latmax":latmax}


#--------------------------------------------------------------------------#
#                    Calculs sur les paramètres                            #               
#--------------------------------------------------------------------------#


def comp_steering(lobj,steering_levels,uvmean_box,filin,inst,indf,res,domtraj,basetime,parfilt,filtapply,subnproc,pos=""):
    #Compute steering flow at levels steering_levels at the location (lon,lat) - lon, lat can be single or list of values
    #If filtrad > 0.0, the steering flow is computed at equivalent resolution filtrad (in km),
    #filin is the input file, and indf the inputdef data
    #uvmean_box : if >0, the size of the box (degrees) to compute the mean value of u,v
    #pos : position where to compute (default=lonc,latc)
    #Output : u,v values of the steering flow (or list, if lon, lat are list)

    nlev = len(steering_levels)
    u_steer0 = []
    v_steer0 = []

    if isinstance(lobj,list):
        npts=len(lobj) #list of values lon,lat
        uniq=False
        tobj=lobj
    else:
        npts=1 #unique value
        uniq=True
        tobj = [lobj]

    #Extraction of total steering field
    lev=steering_levels[0]
    if nlev==1 and lev<1e-6: #valeur nulle
        u_steer0=[0.0 for ivi in range(npts)]
        v_steer0=[0.0 for ivi in range(npts)]
    else:
        ch='u'+str(int(lev))
        fu = Inputs.extract_data(filin,inst,indf,ch,domtraj,res[ch],basetime,subnproc,filtrad=parfilt[ch]*filtapply)
        ch='v'+str(int(lev))
        fv = Inputs.extract_data(filin,inst,indf,ch,domtraj,res[ch],basetime,subnproc,filtrad=parfilt[ch]*filtapply)
        for ivi in range(1,nlev):
            lev=steering_levels[ivi]
            ch='u'+str(int(lev))
            fu.operation("+",Inputs.extract_data(filin,inst,indf,ch,domtraj,res[ch],basetime,subnproc,filtrad=parfilt[ch]*filtapply))
            ch='v'+str(int(lev))
            fv.operation("+",Inputs.extract_data(filin,inst,indf,ch,domtraj,res[ch],basetime,subnproc,filtrad=parfilt[ch]*filtapply))
        fu.operation("/",nlev)
        fv.operation("/",nlev)

        for ivi in range(npts):
            #Apply mean value of u,v in uvmean_box
            if pos=="o":
                lon1=tobj[ivi].traps["olon"]
                lat1=tobj[ivi].traps["olat"]
            else:
                lon1=tobj[ivi].lonc
                lat1=tobj[ivi].latc
            u_steer0.append(np.mean(box_values(fu,lon1,lat1,uvmean_box)))
            v_steer0.append(np.mean(box_values(fv,lon1,lat1,uvmean_box)))

    if uniq==1:
        u_steer=u_steer0[0]
        v_steer=v_steer0[0]
        #print("STEERING FLOW WITH BOX: ", u_steer,v_steer)
        #print("STEERING FLOW AT POINT: ", fu.getvalue_ll(tlon[ivi],tlat[ivi],interpolation="linear"), fv.getvalue_ll(tlon[ivi],tlat[ivi],interpolation="linear"))
    else:
        u_steer=u_steer0
        v_steer=v_steer0
    
    return u_steer,v_steer

def box_values(fld,lon,lat,ssbox):
    #Computes the box of values of fld
    #of size ssbox (degrees) around lon,lat
    #Output : a numpy array or a list

    #Dimensions of fld
    X = fld.geometry.dimensions['X']
    Y = fld.geometry.dimensions['Y']

    res=get_res(fld)

    if ssbox==0:
        box=[fld.getvalue_ll(lon,lat,interpolation="linear")]
    else:
        ssbox=ssbox+2*res
        maxsizeX=int(round((ssbox/res)/2.0))
        maxsizeY=int(round((ssbox/res)/2.0))
    
        if (fld.geometry.point_is_inside_domain_ll(lon+2*res, lat+2*res) and fld.geometry.point_is_inside_domain_ll(lon-2*res, lat-2*res)): 

            pts=fld.geometry.nearest_points(lon,lat,request=dict(n="2*2"))

            imin=max(1,pts[0][0]-maxsizeX)
            imax=min(pts[2][0]+maxsizeX,X-1)
            jmin=max(1,pts[0][1]-maxsizeY)
            jmax=min(pts[3][1]+maxsizeY,Y-1)
            box=fld.data[jmin:jmax,imin:imax]
        else:
            box=[fld.getvalue_ll(lon,lat,interpolation="linear")]

    return box

#--------------------------------------------------------------------------#
#                    Calculs des critères de qualité                       #                                     
#--------------------------------------------------------------------------#
def adv_quality(trackpar, AVOr, obj1, obj2, res, pdmin):
    
    '''Calcule la qualité d'advection.

    Entrées:
        obj1, obj2 : les deux objets à t et t+dt respectivement.
        res : la resolution de la grille.
        pdmin : le rapport d'erreur spatiale.

    Sortie : la qualité d'advection.'''
    
    #dyy = res*comp_length(obj1.lonc,obj1.latc,obj1.lonc,obj1.latc+1)
    dyy = comp_length(obj1.lonc,obj1.latc,obj1.lonc,obj1.latc+res)

    # Récupération des données d'AVO
    AVO1 = obj1.traps[trackpar]
    AVO2 = obj2.traps[trackpar]

    if trackpar[0:2]=="rv": #conversion to absolute vorticity
        AVO1 = AVO1+2*omega*np.sin(obj1.latc*np.pi/180.0)
        AVO2 = AVO2+2*omega*np.sin(obj2.latc*np.pi/180.0)
   
    # Récupération des données de vent
    u1 = (obj1.traps["u_steer1"]+obj2.traps["u_steer1"])/2.0
    v1 = (obj1.traps["v_steer1"]+obj2.traps["v_steer1"])/2.0
    u2 = (obj1.traps["u_steer2"]+obj2.traps["u_steer2"])/2.0
    v2 = (obj1.traps["v_steer2"]+obj2.traps["u_steer2"])/2.0

    # Calcul des advections
    objt1 = obj1.advect(obj2.time,u1,v1)
    lon_advt1 = objt1.lonc
    lat_advt1 = objt1.latc
    objt2 = obj1.advect(obj2.time,u2,v2)
    lon_advt2 = objt2.lonc
    lat_advt2 = objt2.latc

    # Calcul des erreurs
    d1 = dyy*((AVO1 + AVO2)/AVOr)
    d2 = comp_length(obj2.lonc, obj2.latc, lon_advt1, lat_advt1)
    d3 = comp_length(obj2.lonc, obj2.latc, lon_advt2, lat_advt2)
    d4 = comp_length(lon_advt1, lat_advt1, lon_advt2, lat_advt2)
    
    err= max(0.0,d2 + d3 -d4 -d1)
    
    # Calcul de la qualité
    quality = max(0.0,(dyy/(dyy + err) -pdmin)/ (1-pdmin))
    #print("SEUIL DISTANCE :",(dyy*(1-pdmin)/pdmin),"(pdmin=",pdmin,", dy=",dyy,")")
    #print("adv_quality: ", quality, "(",d2, d3, -d4, -d1, err ,")")
    #if quality > 0.0:
        #print("adv quality:",quality)
        #print([d2,d3,-d1,-d4])

    return quality

#--------------------------------------------------------------------------#
def avo_quality(trackpar,obj1,obj2, pvmin):

    '''Calcule l'erreur du tourbillon absolu et renvoie la qualité de cette
    erreur.'''
    '''Entrées : 
        obj1, obj2 : les deux objets à t et t+dt respectivement.
        pvmin : le rapport d'erreur de tourbillon


    Sortie : la qualité d'AVO.'''
    
    # Récupération des données d'AVO
    AVO1 = obj1.traps[trackpar]
    AVO2 = obj2.traps[trackpar]

    if trackpar[0:2]=="rv": #conversion to absolute vorticity
        AVO1 = AVO1+2*omega*np.sin(obj1.latc*np.pi/180.0)
        AVO2 = AVO2+2*omega*np.sin(obj2.latc*np.pi/180.0)

    # Calcul de la variation d'AVO
    var = AVO1/AVO2
    if abs(var) > 1:
        var = 1/var
    # Calcul de la qualité
    quality = max(0,(var -pvmin)/(1 -pvmin))

    #if quality > 0.0:
    #    print("avo quality:",quality)
        
    return quality

#--------------------------------------------------------------------------#
def pseudo_pos_quality(obja, objb, objc, res, pdmin):
    
    '''Computes 2nd order position quality

    Entrées:
        obja, objb, objc : trois objets consécutifs
        res : la resolution de la grille.
        pdmin : le rapport d'erreur spatiale.

    Sortie : la qualité d'advection.'''
    
    #dyy = res*comp_length(objb.lonc,objb.latc,objb.lonc,objb.latc+1)
    dyy = comp_length(objb.lonc,objb.latc,objb.lonc,objb.latc+res)

    # Calcul des erreurs
    #print("POS :",dyy,comp_length(2*obja.lonc-objb.lonc, 2*obja.latc-objb.latc, objc.lonc, objc.latc))
    err = max(0.0,comp_length(2*obja.lonc-objb.lonc, 2*obja.latc-objb.latc, objc.lonc, objc.latc) - dyy)
    
    # Calcul de la qualité
    quality = max(0.0,(dyy/(dyy + err) -pdmin)/ (1-pdmin))

    return quality

#--------------------------------------------------------------------------#
#                            Min/max searches                              #
#--------------------------------------------------------------------------#               
def local_minima(array2d):

    '''La fonction renvoie les coordonnées et valeurs de l'ensemble des minima
    locaux dans un tableau

    Entrée : un tableau de données
    Sortie : un tuple contenant les indices des minimums locaux '''

    a = ((array2d <= np.roll(array2d,  1, 0)) &
         (array2d <= np.roll(array2d, -1, 0)) &
         (array2d <= np.roll(array2d,  1, 1)) &
         (array2d <= np.roll(array2d, -1, 1)))

    minx = 0  # Indice de la 1ère ligne
    maxx = np.shape(array2d)[0]-1  # Indice de la dernière ligne
    miny = 0  # Indice de la première colonne
    maxy = np.shape(array2d)[1]-1  # Indice de la dernière colonne

    # Pour ne pas traiter les bords de la grille
    # Attribue aux éléments de la première ligne la valeur False pour filtrer
    a[minx, :] = False
    a[maxx, :] = False  # Idem pour les éléments de la dernière ligne
    a[:, miny] = False  # Et les éléments de la première colonne
    a[:, maxy] = False  # Et la dernière colonne

    # les indices (x et y) correspondants

    indices = a.nonzero()  # Renvoie les indices des minimums locaux

    return indices
#--------------------------------------------------------------------------#
def comp_length(lon1, lat1, lon2, lat2):

    '''Computes the distance between two points on the globe.

    Input: point coordinates (longitudes, latitudes).
    Output : distance (km) '''

    #Conversion to the same coordinate system
    if lon1 > 180:
        lon1 = lon1 - 360
    if lon2 > 180:
        lon2 = lon2 - 360

    #Conversion to radians
    lon1 = lon1 * np.pi / 180
    lon2 = lon2 * np.pi / 180
    lat1 = lat1 * np.pi / 180
    lat2 = lat2 * np.pi / 180

    a = np.sin(lat1) * np.sin(lat2)

    b = np.cos(lat1) * np.cos(lat2) * np.cos(lon2 - lon1)
    
    ab=min(1.0,a+b)
    ab=max(-1.0,ab)
    #if (ab)>1.0 or (ab)<-1.0:
    #    print(a,b,a+b,np.arccos(a+b))
    #    print("problem with arccos??")

    return (r/1000.0)*np.arccos(ab)

def comp_deriv(fld,lon,lat,ss,domtraj):
    #Computes derivative of fld in the 4 directions at a distance of ss(degrees)
    #If the distance gets out of the domain, we compute the derivative at the frontier
    #Output : list of 4 derivative values

    deriv=[]
    #East
    if abs(domtraj["lonmax"]-lon)>ss:
        deriv.append(fld.getvalue_ll(lon+ss,lat)-fld.getvalue_ll(lon,lat))
    else:
        deriv.append((fld.getvalue_ll(domtraj["lonmax"],lat)-fld.getvalue_ll(lon,lat))*(ss/abs(domtraj["lonmax"]-lon)))
    #West
    if abs(domtraj["lonmin"]-lon)>ss:
        deriv.append(fld.getvalue_ll(lon+ss,lat)-fld.getvalue_ll(lon,lat))
    else:
        deriv.append((fld.getvalue_ll(domtraj["lonmin"],lat)-fld.getvalue_ll(lon,lat))*(ss/abs(domtraj["lonmin"]-lon)))
    #North
    if abs(domtraj["latmax"]-lat)>ss:
        deriv.append(fld.getvalue_ll(lon,lat+ss)-fld.getvalue_ll(lon,lat))
    else:
        deriv.append((fld.getvalue_ll(lon,domtraj["latmax"])-fld.getvalue_ll(lon,lat))*(ss/abs(domtraj["latmax"]-lat)))
    #South
    if abs(domtraj["latmin"]-lat)>ss:
        deriv.append(fld.getvalue_ll(lon,lat-ss)-fld.getvalue_ll(lon,lat))
    else:
        deriv.append((fld.getvalue_ll(lon,domtraj["latmin"])-fld.getvalue_ll(lon,lat))*(ss/abs(domtraj["latmin"]-lat)))

    return deriv
#--------------------------------------------------------------------------#
def get_res(fld):
    #Input : an epygram field on a regular latlon grid
    #Ouput: resolution of the grid

    (lons, lats) = fld.geometry.get_lonlat_grid()
    ll=lons[0,:].tolist()

    return abs(ll[1]-ll[0])

#--------------------------------------------------------------------------#
def find_allmin(fld,dist=maxval,thr=0.0):
    #Finds the local minima in field fld, that are below thr
    #if several minima are inside a dist radius (km),
    #we keep the one with highest absolute value

    indexf = local_minima(fld.data)
    tlat=[]
    tlon=[]
    tval=[]

    for ivj in range(np.shape(indexf)[1]):
        lon,lat = fld.geometry.ij2ll(indexf[1][ivj],indexf[0][ivj])
        val = fld.getvalue_ll(lon,lat)

        if val<thr:

            if (dist<maxval-1):
                indc=[] #list of indices that are in the radius
                for ivi in range(len(tlat)):
                    if comp_length(lon,lat,tlon[ivi],tlat[ivi])<dist:
                        indc.append(ivi)
                    #indc: liste des indices qui sont dans le rayon dist

                if len(indc)>0:
                    vmax = val
                    imax = -1
                    for ivi in range(len(indc)):
                        #if abs(tval[indc[ivi]])>vmax:
                        if tval[indc[ivi]]<vmax:
                            imax = indc[ivi]
                            vmax = abs(tval[indc[ivi]])
                    #vmax, imax: valeur maximal

                    if imax==-1: #le max est le nouveau - on supprime les valeurs dans indc
                        nsup=len(indc)
                        for ivi in range(nsup-1,-1,-1):
                            del tlat[indc[ivi]]
                            del tlon[indc[ivi]]
                            del tval[indc[ivi]]
                        tlat.append(lat)
                        tlon.append(lon)
                        tval.append(val)
                    #else: (imax>=0) # le max est deja existant - on ajoute rien

                else: #len(indc>0) : pas de max detecte dans le rayon
                    tlat.append(lat)
                    tlon.append(lon)
                    tval.append(val)

            else: #pas de critere de distance
                tlat.append(lat)
                tlon.append(lon)
                tval.append(val)

    return tlon,tlat,tval

#--------------------------------------------------------------------------#
def find_allmax(fld,dist=maxval,thr=0.0):
        #Finds the local maxima in field fld that are above thr
        #if several minima are inside a dist radius (km),
        #we keep the one with highest absolute value

    fld.data=-fld.data
    tlon,tlat,tval = find_allmin(fld,dist=dist,thr=thr)
    fld.data=-fld.data

    tval1 = [-tval[i] for i in range(len(tval))]

    return tlon,tlat,tval1

#--------------------------------------------------------------------------#
def find_allextr(sign, fld, dist=maxval, thr=0.0):
    #Finds the local maxima if sign=1, minima if sign=-1, respectively
    
    if sign==-1:
        tlon, tlat, tval = find_allmin(fld,dist=dist,thr=thr)
    elif sign==1:
        tlon, tlat, tval = find_allmax(fld,dist=dist,thr=thr)
    else:
        print('Error in Tools.find_allextr: wrong sign (should be -1 or 1)')

    return tlon,tlat,tval

def list_abovebelow(sign,fld,lon0,lat0,thr,ss):
    #List all the points in fld that are above (sign=1) or below (sign=-1) thr in a square of size ss centred on lon0,lat0

    lon=[]
    lat=[]
    val=[]

    #Dimensions of fld
    X = fld.geometry.dimensions['X']
    Y = fld.geometry.dimensions['Y']
    res=get_res(fld)

    maxsize=int(round((ss/res)/2.0))

    if (fld.geometry.point_is_inside_domain_ll(lon0+res, lat0+res) and
       fld.geometry.point_is_inside_domain_ll(lon0-res, lat0-res)): 

        pts=fld.geometry.nearest_points(lon0,lat0,request=dict(n="2*2"))

        imin=max(1,pts[0][0]-maxsize)
        imax=min(pts[2][0]+maxsize,X-1)
        jmin=max(1,pts[0][1]-maxsize)
        jmax=min(pts[3][1]+maxsize,Y-1)
        tab=fld.data[jmin:jmax,imin:imax]

        for ivj in range(jmin,jmax):
            for ivi in range(imin,imax):
                if sign*tab[ivj-jmin,ivi-imin]>sign*thr:
                    val.append(tab[ivj-jmin,ivi-imin])
                    x,y = fld.geometry.ij2ll(ivi,ivj)
                    lon.append(x)
                    lat.append(y)

    return lon, lat, val
#--------------------------------------------------------------------------#
def find_amin(fld, lon0 ,lat0 ,ss=0 , thr=maxval ):
    #Finds the minimum absolute value and position in field fld to the point (lon0,lat0)
    #In a square of size ss degrees (if ss not specified or 0, the whole
    #domain is searched for) and which value is below thr
    #Output : gook : boolean to verify a minimum has been found
    #Output : lat and lon are gridpoints of fld 
    #(we do not recommend to apply smooth_min afterwards because the result is not a local minimum)

    gook=False
    val=0.0
    dmin = traject.missval
    valmin = traject.missval
    lon1=traject.missval
    lat1=traject.missval

    #Dimensions of fld
    X = fld.geometry.dimensions['X']
    Y = fld.geometry.dimensions['Y']

    res=get_res(fld)

    if ss==0:
        maxsizeX=int(X/2)-1
        maxsizeY=int(Y/2)-1
    else:
        #ss=ss+2*res
        maxsizeX=int(round((ss/res)/2.0))
        maxsizeY=int(round((ss/res)/2.0))
    
    if (fld.geometry.point_is_inside_domain_ll(lon0+2*res, lat0+2*res) and
       fld.geometry.point_is_inside_domain_ll(lon0-2*res, lat0-2*res)): 

        pts=fld.geometry.nearest_points(lon0,lat0,request=dict(n="2*2"))

        valmin = maxval

        imin=max(1,pts[0][0]-maxsizeX)
        imax=min(pts[2][0]+maxsizeX,X-1)
        jmin=max(1,pts[0][1]-maxsizeY)
        jmax=min(pts[3][1]+maxsizeY,Y-1)
        tab=fld.data[jmin:jmax,imin:imax]

        valmin = np.amin(tab,axis=None)
        indexf = np.unravel_index(np.argmin(tab,axis=None), tab.shape)
        
        if np.size(indexf) > 0: #Finds the local minimum with lower value
            lon,lat = fld.geometry.ij2ll(indexf[1]+imin,indexf[0]+jmin)
            dist=comp_length(lon0,lat0,lon,lat)
            val = fld.getvalue_ll(lon,lat)
            if (val<thr):
                gook=True
                dmin=dist
                valmin=val
                lon1=lon
                lat1=lat

    return lon1, lat1, valmin, dmin, gook


def find_locmin(fld, lon0 ,lat0 ,ss=0 , thr=maxval ):
    #Finds the closest local minimum in field fld to the point (lon0,lat0)
    #In a square of size ss degrees (if ss not specified or 0, the whole
    #domain is searched for)
    #and which value is below thr
    #Output : gook : boolean to verify a minimum has been found
    #Output : lat and lon are gridpoints of fld (to interpolate, use smooth_min() afterwards

    val=0.0
    lon1=traject.missval
    lat1=traject.missval

    #Dimensions of fld
    X = fld.geometry.dimensions['X']
    Y = fld.geometry.dimensions['Y']

    res=get_res(fld)

    if ss==0:
        maxsize=min(int(X/2),int(Y/2)) 
    else:
        #The function local_minima exclude the boundary zones
        ss=ss+2*res
        maxsize=(ss/res)/2.0
    
    if (fld.geometry.point_is_inside_domain_ll(lon0+2*res, lat0+2*res) and
       fld.geometry.point_is_inside_domain_ll(lon0-2*res, lat0-2*res)): 

        pts=fld.geometry.nearest_points(lon0,lat0,request=dict(n="2*2"))

        dmin=maxval
        valmin=maxval
        outok = False

        ivi = int(maxsize)
        #print("ss=",ss,"maxsize=",maxsize,"ivi=",ivi)
        imin=max(1,pts[0][0]-ivi)
        imax=min(pts[2][0]+ivi,X-1)
        jmin=max(1,pts[0][1]-ivi)
        jmax=min(pts[3][1]+ivi,Y-1)
        #print("imin,imax,jmin,jmax",imin,imax,jmin,jmax)
        #print("llmin,llmax",fld.geometry.ij2ll(imin,jmin),fld.geometry.ij2ll(imax,jmax))
        tab=fld.data[jmin:jmax,imin:imax]
        indexf = local_minima(tab)
        
        if np.size(indexf) > 0: #Finds the closest local minimum
            outok=True
            nfound = np.shape(indexf)[1]
            for ivj in range(nfound):
                lon,lat = fld.geometry.ij2ll(indexf[1][ivj]+imin,indexf[0][ivj]+jmin)
                dist=comp_length(lon0,lat0,lon,lat)
                val = fld.getvalue_ll(lon,lat)
                #print("LONLAT2",lon,lat,val,tab[indexf[0],indexf[1]])
                if (val<thr) and (dist<dmin):
                    dmin=dist
                    valmin=val
                    lon1=lon
                    lat1=lat
                    #print("Potential kernel : ",lon1,lat1,valmin, dist)

        gook = (dmin<maxval) and (valmin<thr)

    else:
        gook = False #too close to the domain boundaries
        dmin = traject.missval
        valmin = traject.missval

    ### MODIF FP
    #DEBUG - Plot champs et des points candidats (pour voir)
    if False:
    #if True:
        proj = ccrs.PlateCarree()
        plon,plat=fld.geometry.get_lonlat_grid()
        fig1=plt.figure(figsize=(30,15))
        ax1 = plt.axes(projection = proj)
        ax1.set_extent([0.0, 10.0, 39.0, 46.0])
        ax1.coastlines()
        #clev=[x*100 for x in range(950,1050,2)]
        clev=[x for x in range(-20,0,2)]
        ax1.contour(plon,plat,fld.data,clev)
        ax1.scatter(lon0,lat0,marker="o",c="b")
        ax1.scatter(lon1,lat1,marker="o",c="r")
        fig1.savefig('mslp-fig-.png')
        plt.close()
    ### MODIF FP

    return lon1, lat1, valmin, dmin, gook

#--------------------------------------------------------------------------#
def find_locmax(fld, lon0, lat0, ss=0, thr=maxval):
    #Finds the closest local max in field fld to the point (lon0,lat0)
    #which value is above thr
    
    fld.data=-fld.data
    if thr==maxval:
        thr=-maxval
    lon1, lat1, val, dmin, gook = find_locmin(fld,lon0,lat0,ss=ss,thr=-thr)
    fld.data=-fld.data

    return lon1, lat1, -val, dmin, gook

#--------------------------------------------------------------------------#
def find_amax(fld, lon0, lat0, ss=0, thr=maxval):
    #Finds a local max in field fld close to the point (lon0,lat0)
    #which value is above thr and with the maximum value
    
    fld.data=-fld.data
    if thr==maxval:
        thr=-maxval
    lon1, lat1, val, dmin, gook = find_amin(fld,lon0,lat0,ss=ss,thr=-thr)
    fld.data=-fld.data

    return lon1, lat1, -val, dmin, gook

#--------------------------------------------------------------------------#
def find_locextr(sign, fld, lon0, lat0, ss=0, thr=maxval):
    #Finds the closest local max if sign=1, min if sign=-1, respectively
    
    if sign==-1:
        lon1, lat1, valmin, dmin, gook = find_locmin(fld,lon0,lat0,ss=ss,thr=thr)
    elif sign==1:
        lon1, lat1, valmin, dmin, gook = find_locmax(fld,lon0,lat0,ss=ss,thr=thr)
    else:
        print('Error in Tools.find_locextr: wrong sign (should be -1 or 1)')

    return lon1, lat1, valmin, dmin, gook

#--------------------------------------------------------------------------#
def find_aextr(sign, fld, lon0, lat0, ss=0, thr=maxval):
    #Finds the local max if sign=1, min if sign=-1, respectively
    #which value is above thr and with the highest value
    
    if sign==-1:
        lon1, lat1, valmin, dmin, gook = find_amin(fld,lon0,lat0,ss=ss,thr=thr)
    elif sign==1:
        lon1, lat1, valmin, dmin, gook = find_amax(fld,lon0,lat0,ss=ss,thr=thr)
    else:
        print('Error in Tools.find_extr: wrong sign (should be -1 or 1)')

    return lon1, lat1, valmin, dmin, gook

#--------------------------------------------------------------------------#
def smooth_min(fld,glon,glat):
    #Finds the subgrid minima in fld,
    #given glon, glat is a grid minima 
    #Using Tech Memo ECMWF 386 formula

    res=get_res(fld)

    pm = fld.getvalue_ll(glon,glat)
    pmA = fld.getvalue_ll(glon-res,glat)
    pmB = fld.getvalue_ll(glon,glat+res)
    pmC = fld.getvalue_ll(glon,glat-res)
    pmD = fld.getvalue_ll(glon+res,glat)

    d1 = pmD-pmA
    d2 = pmB-pmC

    s1 = pmA+pmD-2*pm
    s2 = pmB+pmC-2*pm

    if not s1*s2==0:
        lon = glon - 0.5*res*d1/s1
        lat = glat - 0.5*res*d2/s2
        val = pm -(d1*d1/s1 + d2*d2/s2)/8.0
    else:
        lon = glon
        lat = glat
        val = pm

    return lon, lat, val

#--------------------------------------------------------------------------#
def smooth_max(fld, glon, glat):
    #Finds the subgrid maxima in fld,
    
    fld.data=-fld.data
    lon, lat, val = smooth_min(fld, glon, glat)
    fld.data=-fld.data

    return lon, lat, -val

#--------------------------------------------------------------------------#
def smooth_extr(sign, fld, glon, glat):
    #Finds the subgrib maxima if sign=1, min if sign=-1, respectively
    
    if sign==-1:
        lon, lat, val = smooth_min(fld,glon,glat)
    elif sign==1:
        lon, lat, val = smooth_max(fld,glon,glat)
    else:
        print('Error in Tools.find_locextr: wrong sign (should be -1 or 1)')

    return lon, lat, val

#--------------------------------------------------------------------------#
def maskrad(lons,lats,lonc,latc,rad):
    #Returns the matrix tab (x=lons, y=lats) which are inside (tab=1)
    #or outside (tab=0) the circle of centre (lonc,latc) of radius rad (km)
    ##Input:
    #lons, lats: longitudes and latitudes of the domain on which to compute the mask
    #(lonc,latc) : central point
    #rad: distance to apply around the central point
    #Output: tab: array of shape (len(lons), len(lats)), 1 if in the mask, 0 otherwise

    polygonSides = 180 #number of angles to define the circle

    aeqd_proj = "+proj=aeqd +lat_0="+str(latc)+" +lon_0="+str(lonc)+" +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"
    latlon_proj = "+proj=latlon"
    proj = pyproj.Transformer.from_crs(pyproj.CRS.from_string(aeqd_proj),pyproj.CRS.from_proj4(latlon_proj))
    angles = np.linspace(0, 2 * np.pi, polygonSides, endpoint=False)
    points_list = [proj.transform( np.sin(a) * rad*1000.0, np.cos(a) *rad*1000.0) for a in angles]

    pol = Polygon(points_list)

    nlon = len(lons)
    nlat = len(lats)
    #print("Tools.maskrad dimensions: nlon=",nlon," nlat=",nlat," - nlon*nlat=",nlon*nlat)
    tab = np.zeros((nlat,nlon),dtype=int)
    for ivx in range(nlon):
        for ivy in range(nlat):
            if pol.contains(Point(lons[ivx],lats[ivy])):
                tab[ivy,ivx] = 1

    return tab

#--------------------------------------------------------------------------#
def bufferrad(lons,lats,tab1,rad):
    #Returns the matrix tab (x=lons, y=lats) which are inside (tab=1)
    #or outside (tab=0) an extension of tab1 by the radius rad (km)
    ##Input:
    #lons, lats: longitudes and latitudes of the domain on which to compute the mask
    #tab : first guess (0 or 1) around which the radius is extended
    #rad: distance to apply around the existing zone in tab1 (km)
    #Output: tab: array of shape (len(lons), len(lats)), 1 if in the mask, 0 otherwise

    nlat=len(lats)
    nlon=len(lons)
    listseuils=[0.5]
    latlon_proj = "+proj=latlon"
    tab2=tab1
    reslon=lons[1]-lons[0]
    reslat=lats[1]-lats[0]

    #print("Tools.bufferrad dimensions: nlon=",nlon," nlat=",nlat," - nlon*nlat=",nlon*nlat)
    #Find existing contour in tab
    if np.max(tab2) > np.min(tab2):
        cs = plt.contour(tab2,levels=listseuils,cmap='jet')
        paths=cs.collections[0].get_paths()
        #print("Tools.bufferrad ; number of objects: npaths=",len(paths))
        for p in paths:
            polygon=Polygon(p.vertices)

            #Projection du polygone en x,y
            xt, yt = polygon.exterior.coords.xy
            x0,y0=polygon.centroid.xy
            x=x0[0]
            y=y0[0] 
            lonc=lons[int(max(0,min(x,nlon-1)))]
            latc=lats[int(max(0,min(y,nlat-1)))]

            tdx=comp_length(lonc,latc,lonc+reslon,latc)
            tdy=comp_length(lonc,latc,lonc,latc+reslat)
            marg = 3 # margin
            minxt=max(0,int(min(xt)-rad/tdx-marg))
            maxxt=min(nlon-1,int(max(xt+rad/tdx+1+marg)))
            minyt=max(0,int(min(yt)-rad/tdy-marg))
            maxyt=min(nlat-1,int(min(yt+rad/tdy)+1+marg))
            #print("max/min: ", minxt, maxxt,minyt,maxyt)

            aeqd_proj = "+proj=aeqd +lat_0="+str(latc)+" +lon_0="+str(lonc)+" +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m"
            proj = pyproj.Transformer.from_crs(pyproj.CRS.from_proj4(latlon_proj),pyproj.CRS.from_string(aeqd_proj))
            xy_list = [proj.transform( interp(x,int(x),int(x)+1,lons[int(x)],lons[int(x)]+reslon),
                interp(y,int(y),int(y)+1,lats[int(y)],lats[int(y)]+reslat)  ) for x,y in p.vertices]

            #Buffer polygon x,y (km)
            xy_pol = Polygon(xy_list).buffer(rad*1000.0).exterior.xy

            #Retour en latlon
            proj_inv = pyproj.Transformer.from_crs(pyproj.CRS.from_string(aeqd_proj),pyproj.CRS.from_proj4(latlon_proj))
            pol = Polygon([proj_inv.transform(xy_pol[0][ivi], xy_pol[1][ivi]) for ivi in range(len(xy_pol[0]))])

            for ivx in range(minxt,maxxt):
                for ivy in range(minyt,maxyt):
                    if pol.contains(Point(lons[ivx],lats[ivy])):
                        tab2[ivy,ivx] = 1

    return tab2


def interp(x,x1,x2,lon1,lon2):

    return lon1 + (lon2-lon1)*(x-x1)/(x2-x1)
#--------------------------------------------------------------------------#

#--------------------------------------------------------------------------#
#                      Computation of diagnostics                          #
#--------------------------------------------------------------------------#
#--------------------------------------------------------------------------#

class diagdef:
#Class that defines a diagnostic
#and its representation as a string

    def __init__(self,ds,Hn,ss,rd):
        #Initialize diag
        #If ds has the format PARAMETER_KIND_AREA, then all information is provided
        #If ds has the format PARAMETER_AREA, then Hn is required to guess KIND

        ds2=ds.split("_")
        
        if len(ds2)==1: #only parameter ...
            par, sign = get_parsign(ds2[0], Hn)
            if sign==1:
                kind="max"
            else:
                kind="min"
            area="null"
        elif len(ds2)==2:
            par, sign = get_parsign(ds2[0], Hn)
            if sign==1:
                kind="max"
            else:
                kind="min"
            area=ds2[1]
        elif len(ds2)==3:
            par=ds2[0]
            kind=ds2[1]
            area=ds2[2]

        #Completing par
        self.par = par

        #Completing kind and associated value if relevant
        if kind[0]=="c":
            val=float(kind[2:len(kind)].replace("p","."))
            self.kind=kind[0:2]
            self.kindval=val
        else:
            self.kind=kind
            self.kindval=traject.missval

        #Completing area and associated value if relevant
        if area[0]=="s" or area[0]=="r":
            self.area=area[0]
            if len(area)>1:
                val=float(area[1:len(area)].replace("p","."))
                self.areaval=val
                area2=area
            elif area[0]=="s" and ss>0.0: 
                self.areaval=ss
                area2='s'+res2str(ss, 1)
            elif area[0]=="r" and rd>0.0: 
                self.areaval=rd
                area2='r'+res2str(rd, 1)
            else:
                print("Some information is missing to compute the area of "+par+" - Abort")
                exit()
        elif area=="p" or area=="o":
            self.area=area
            area2=area
        else:
            print("area not defined in Tools.diagdef() - WARNING")
            area2=""

        #self.strg = full string name of the diagnostic
        self.strg = par + "_" + kind + "_" + area2

        #Nice name and unit for plotting (to be completed)
        self.nicename = ""
        self.unit = ""
        self.plot_fct = 1.0
        if self.par[0:2]=="rr":
            dur=self.par[-1]
            t=self.par[2:-1]
            self.nicename = "rainfall "+t+dur
            self.unit = "mm"
            self.plot_unit = self.unit
        elif self.par=="mslp":
            self.nicename="mslp"
            self.unit = "Pa"
            self.plot_unit = "hPa"
            self.plot_fct = 0.01
        elif self.par=="btir":
            self.nicename="IR Brightness Temperature"
            self.unit = "K"
            self.plot_unit = self.unit
        else:
            print("you may declare nicename and unit in Tools.diagdef()")

        return

    def equals(self,diag):

        eq = self.par==diag.par and self.kind==diag.kind and self.area==diag.area

        return eq

    def compares(self,diag):
        #Checks if two diagnostics are comparable (same parameter, same kind
        # no matter of the area around which it has been computed)

        eq = self.par==diag.par and self.kind==diag.kind

        return eq

def guess_diag(strdiag,Hn):
    #Guesses the diagnostic from a string diagnostic
    #Input: string diagnostic value
    #Output: a diagdef object ... with possibly some fake values

    diag=diagdef(strdiag,Hn,-999.0,-999.0)    

    return diag

def make_diags(ldiag,obj,filin,inst,indf,domtraj,Hn,res,basetime,olon,olat,subnproc,**kwargs):
    #olon, olat: origin points (tracking parameter)

    for diag in ldiag:
        obj.diags.append(diag.strg)
        if "parfilt" in kwargs and "filtapply" in kwargs:
            if kwargs["filtapply"]==0:
                filtrad=0.0
            else:
                parfilt=kwargs["parfilt"]
                filtrad=parfilt[diag.par]*kwargs["filtapply"]
        else:
            filtrad=0.0
        fld=Inputs.extract_data(filin,inst,indf,diag.par,domtraj,res[diag.par],basetime,subnproc,filtrad=filtrad)
        if diag.area=="o":
            val= fld.getvalue_ll(olon,olat,interpolation="linear")
            setattr(obj,diag.strg,[olon, olat, val])
        elif diag.area=="p":
            val= fld.getvalue_ll(obj.lonc,obj.latc,interpolation="linear")
            setattr(obj,diag.strg,[obj.lonc, obj.latc, val])
        elif diag.area=="s": #On cherche valeur extreme dans carre de cote ss et on lisse
            if diag.kind=="max" or diag.kind=="min":
                if diag.kind=="max":
                    sign=1
                else:
                    sign=-1
                lon,lat,val,dist,gook2 = find_aextr(sign,fld,obj.lonc,obj.latc,ss=diag.areaval)
                setattr(obj,diag.strg,[lon, lat, val])
            elif diag.kind=="mean":
                print("option kind=mean A CODER DANS Tools.make_diags")
                exit()
            elif diag.kind=="cx" or diag.kind=="cn":
                if diag.kind=="cx":
                    sign=1
                else:
                    sign=-1
                lon, lat, val = list_abovebelow(sign,fld,obj.lonc,obj.latc,diag.kindval,diag.areaval)
                listg=[]
                for ivi in range(len(lon)):
                    listg.extend([lon[ivi],lat[ivi],val[ivi]])
                setattr(obj,diag.strg,listg)

        elif diag.area=="r": #On cherche valeur extreme dans un circle de rayon rd
            print("Option area=r A CODER DANS Tools.make_diags")
            exit()
        else:
            print("ABORT - Wrong option in diagnostic " + diag)
            exit()
    
    return

def res2str(res, d):
    #converts a float res to a string by replacing the '.' by 'p'
    #and avoiding a long chain (no more than d digits)

    eps = 10**(-d-1)

    if res-int(res)<eps: #we assume it is an integer then
        st=str(int(res))
    else:
        st2=str(res)
        st3=st2.split(".")
        st=st3[0]+"p"+st3[1][d-1]

    return st

#--------------------------------------------------------------------------#
#                      Some basic plotting routines                        #
#--------------------------------------------------------------------------#
def plot_F_LL(fld,tlon,tlat,fname,nl=10,vect=[]):
    #Plot field fld and points tlon,tlat - just to see how things are going
    #As an option, plots vector vect
    #Number of contours : nl
    #Output is fname.png

    proj = ccrs.PlateCarree()
    fig1=plt.figure(figsize=(30,15))
    ax1 = plt.axes(projection = proj)
    ax1.gridlines(crs=proj, linewidth=2, color='black',draw_labels=True, alpha=0.5, linestyle='--')
    ax1.coastlines()
    plon,plat=fld.geometry.get_lonlat_grid()
    ax1.contour(plon,plat,fld.data,nl)
    ax1.scatter(tlon,tlat,marker="x",c="k")
    if len(vect)==2*len(tlon):
        for ivi in range(len(tlon)):
            plt.arrow(tlon[ivi],tlat[ivi],vect[2*ivi]/10,vect[2*ivi+1]/10,head_width=0.3)
    fig1.savefig(fname+'.png')
    plt.close(fig1)

#--------------------------------------------------------------------------#
#--------------------------------------------------------------------------#


#--------------------------------------------------------------------------#
#                      Some field processing                               #
#--------------------------------------------------------------------------#

def max_diag(par, traj, indf, dom, res, mintime, maxtime, basetime,subnproc=1):
    #Computes the maximum value of a parameter inside a domain, between mintime and maxtime
    #Inputs:
    #traj: the trajectory corresponding to the data (to get some metadata if needed)
    #par: parameter name
    #indf: inputdef of the files
    #dom: domain (domtraj format), res
    #mintime, maxtime (time format)
    #basetime : used only for fc case (0 by default) - str format
    #Outputs:
    #val, lon, lat, time of the maximum value

    #print("computing max_diag for "+par+" between ",mintime," and ",maxtime)

    #Processing inputdef and get list of files and instants
    if isinstance(indf,str):
        indf2=Inputs.inputdef()
        indf2.read_input(indf)
    else:
        indf2=indf

    if "member" in indf2.__dict__.keys():
        indf2.member = traj.inputdef["member"]

    #Processing instants
    if indf2.origin=="fc":  
        terms = dict()
        terms["init"] = comp_difftime(basetime, mintime)
        terms["final"] = comp_difftime(basetime, maxtime)
        step=1
        terms["step"] = step
        lfile,linst = indf2.get_filinst(terms,basetime)
        while (not Inputs.check_file(lfile[step],indf2,par,par)) and (terms["init"]+step <= terms["final"]):
            step = step + 1
        terms["step"] = step
    else: #analysis or climate
        print("max_diag to be completed for an or cl ...")

    #Initialisations
    val = -999.0
    lon = -999.0
    lat = -999.0
    time = datetime.strptime("2000010100", Inputs.time_fmt)
    it = 0
    
    lfile,linst = indf2.get_filinst(terms,basetime)
    while it<=len(lfile)-1: # and Inputs.check_file(lfile[it],indf2,par,par):

        #Extraction of field (no filtering applied)
        flda = Inputs.extract_data(lfile[it],linst[it],indf2,par,dom,res,basetime,subnproc)
        lon0 = 0.5*(dom["lonmax"]+dom["lonmin"])
        lat0 = 0.5*(dom["latmax"]+dom["latmin"])
        #print(lon0,lat0)
        lon1, lat1, val1, dmin1, gook = find_amax(flda, lon0, lat0)
        if val1 > val:
            val = val1
            lon = lon1
            lat = lat1
            time = linst[it]
        it = it + 1

    return val, lon, lat, time

