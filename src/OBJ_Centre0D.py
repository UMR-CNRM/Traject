#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Traject/Centre0D object definitions

https://github.com/UMR-CNRM/Traject
Traject is disributed under Cecill-C license
(For rights and obligations, see LICENSE.txt)

UMR-CNRM, Meteo-France & CNRS
"""

# =============================================================================
#      GENERIC CYCLONE OBJECT (VALID FOR EXTRATROPICAL AND TROPICAL CYCLONE)
# =============================================================================

#------------------------------------------------------------------------------
                                # Importations
#------------------------------------------------------------------------------

import numpy as np
import json

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
from datetime import datetime
from traject import *
import Tools

#------------------------------------------------------------------------------
                        # Les classes relatives à l'objet
#------------------------------------------------------------------------------

                    ## Classe de l'object météorologique ##

''' Un objet météorologique est défini par sa position longitude, latitude pour
chaque temps. Des variables diagnostics peuvent être associées pour chaque
temps si nécessaire.'''


class ObjectM:

    def __init__(self, ldiag, ltrap,**kwargs):
    #Object constructor
    #If lonc,latc,time are given:
        #Generate object with appropriate conventions
        #Longitude has ]-180,180] range
        #Time has YYYYMMDDHH format

        # Définition
        if "lonc" in kwargs:
            self.lonc = kwargs["lonc"]
            if self.lonc>180.0:
                self.lonc = self.lonc - 360.0
        else:
            self.lonc = missval
        if "latc" in kwargs:
            self.latc = kwargs["latc"]
        else:
            self.latc = missval

        if "time" in kwargs:
            time2 = kwargs["time"]
            if isinstance(time2,str):
                self.time=time2
            else:
                self.time=datetime.strftime(time2,time_fmt)
        else:
            self.time = str(missval)

        self.nameobj=""

        # Diagnostics (optionnal)
        self.diags=ldiag

        # Parameters used by the tracking algorithm (optionnal)
        if len(ltrap)>0:
            self.traps={par:missval for par in ltrap}

    def copy(self):
        #Copy the object

        obj2 = ObjectM(vars(self.diags),vars(self.traps))

        obj2.lonc = self.lonc
        obj2.latc = self.latc
        obj2.time = self.time
        
        if obj2.diag:
            obj2.diags = self.diags
        if obj2.trap:
            obj2.traps = self.traps

        return obj2

    def write(self, outputfile, write_diag=True):
        '''Write the object in a JSON output file
           write_diag: True if the diagnostic values should be stored in the JSON file
'''
        # Initialisation of the dictionary that will be dumped in the json file
        dataout={}

        # Récuperation de l'objet
        t = self.time

        # Récuperation des attributs de définition
        if self.nameobj=="":
            dico = {'lonc': self.lonc, 'latc': self.latc}
        else:
            dico = {'lonc': self.lonc, 'latc': self.latc, 'nameobj':self.nameobj}
        if "diags" in vars(self)>0 and write_diag:
            dico["diags"]=self.diags
            for diag in self.diags:
                dico[diag]=obj.__dict__[diag]

        # Alimentation du dictionnaire par les valeurs souhaitées
        dataout["time:" + str(t)] = dico

        # L'ecriture du dictionnaire dans un fichier Json
        with open(outputfile +'.json', 'w') as outfile:
            json.dump(dataout, outfile, indent = 4)


    def read(self,inputfile):
        '''Reads the track from the json file (including diags)
        '''

        with open(inputfile+'.json') as json_file:
            datain = json.load(json_file)
        
        for dico in datain.items():
            diconame = dico[0]
            if diconame[0:4]=="time":
                t = diconame[5:]
                latc = dico[1]['latc']
                lonc = dico[1]['lonc']
                ldiag = dico[1]
                del ldiag['latc']
                del ldiag['lonc']
                obj = ObjectM(ldiag,[])
                obj.time=t
                obj.latc=latc
                obj.lonc=lonc
                if "nameobj" in dico:
                    obj.nameobj=dico[1]['nameobj']
                for diag in ldiag:
                    setattr(obj,diag,dico[1][diag])
        self = obj


    def check_dom(self,dom):
    #check if the object point is inside the dom domain

        #print(self.lonc,dom["lonmin"],dom["lonmax"])
        b = self.lonc>dom["lonmin"] and self.lonc<dom["lonmax"] and self.latc>dom["latmin"] and self.latc<dom["latmax"] 

        return b

    def advect(self,vtime,u,v,pos=""):
        #Creates a new object (valid at time vtime), after the advection of self by u and v (m/s)
        #If kind="o", we advect traps["olon"],traps["olat"] (raw tracking position)
        #Otherwise we advect lonc, latc (default)

        obj = ObjectM([],[],time=vtime)
        dt = (datetime.strptime(obj.time,time_fmt) - datetime.strptime(self.time,time_fmt)).seconds #Timedelta (in seconds)

        if pos=="o":
            lon1=self.traps["olon"]
            lat1=self.traps["olat"]
        else:
            lon1=self.lonc
            lat1=self.latc

        dlon = Tools.comp_length(lon1,lat1,lon1+1.0,lat1)*1000.0 #Distance (en m) d'un deg de longtitude
        dlat = Tools.comp_length(lon1,lat1,lon1,lat1+1.0)*1000.0 #Distance (en m) d'un deg de latitude

        obj.lonc = lon1 + u*dt/dlon
        obj.latc = lat1 + v*dt/dlat

        return obj

    def comp_mvt(self,obj,pos=""):
        #computes the speed movement between obj (origin) and self
        #output : u,v (m/s)

        dt = (datetime.strptime(obj.time,time_fmt) - datetime.strptime(self.time,time_fmt)).seconds #Timedelta (in seconds)

        if pos=="o":
            lon0=obj.traps["olon"]
            lat0=obj.traps["olat"]
            lon1=self.traps["olon"]
            lat1=self.traps["olat"]
        else:
            lon0 = obj.lonc
            lat0 = obj.latc
            lon1 = self.lonc
            lat1 = self.latc

        dlon = Tools.comp_length(self.lonc,self.latc,self.lonc+1.0,self.latc)*1000.0 #Distance (en m) d'un deg de longtitude
        dlat = Tools.comp_length(self.lonc,self.latc,self.lonc,self.latc+1.0)*1000.0 #Distance (en m) d'un deg de latitude

        if not dt==0.0:
            u_speed = (lon1 - lon0)*dlon/dt
            v_speed = (lat1 - lat0)*dlat/dt
        else:
            u_speed = 0.0
            v_speed = 0.0

        return u_speed, v_speed

    def plot(self, ax):

        ax.plot(self.lonc, self.latc, 'x', markersize=10)

    def tdist(self, obj):
        #computes the t-distance between self and obj
        #such as:
            # if self.time==obj.time, t-dist=distance (km) between self dans obj
            # if self.time==obj.time, t-dist=maxval (km)

        if self.time==obj.time:
            td = Tools.comp_length(self.lonc,self.latc,obj.lonc,obj.latc)
        else:
            td = Tools.maxval

        return td

    def mask(self,lons,lats,diag,diagthr,diagrad,centre):
        #returns the mask of the object, depending on different input parameters
        #Input:
            #lons, lats: longitudes and latitudes of the domain on which to compute the mask
            #diag: if "" then the mask is computed from the object properties (lonc, latc) and diagrad (radius, km),
            # otherwise, diag is a diagnostic name that may be used to compute the mask :
                #diagthr: threshold to apply to the diagnostic (only diag values is above/below diagthr are kept in the mask)
                #diagrad: distance to add to the diagnostic map (if cn or cx kind), or as a circle around the central point
                #centre: "obj" to mask around the object centre, "diag" to mask around the diagnostic centre (if max or min kind)
        #Output: mask: array of shape (len(lons), len(lats)), 1 if in the mask, 0 otherwise

        nlon = len(lons)
        nlat = len(lats)
        tab = np.zeros((nlat,nlon),dtype=int)

        #Initialisation
        if diag=="":
            dkind=""
            ldiag=[]
            lonm=obj.lonc
            latm=obj.latm
        else:
            Hn = (lats[0]>=0.0) #Northern Hemisphere
            dd = Tools.guess_diag(diag,Hn)
            par, parsign = Tools.get_parsign(dd.par,Hn)
            dkind=dd.kind

        #Different diag cases
        if dkind=="": #Circle around the object
            tab = Tools.maskrad(lons,lats,self.lonc,self.latc,diagrad)

        elif dkind[0]=="c": #all points that are below or above a value
            #Take all the values
            ldiag = self.__dict__[diag]
            ivi=0
            while ivi<len(ldiag):
                if parsign*ldiag[ivi+2] >= parsign*diagthr:
                    lon1=self.__dict__[diag][ivi]
                    lat1=self.__dict__[diag][ivi+1]
                    xi, yi, chk = Tools.ll_to_ij(lons,lats,lon1,lat1)
                    if chk:
                        tab[yi,xi] = 1
                ivi = ivi + 3
                #Add buffer around tab:
            if diagrad>0.0:
                tab = Tools.bufferrad(lons,lats,tab,diagrad)

        else: #circle around the object or the diag center
            ldiag = self.__dict__[diag]
            if parsign*ldiag[2] >= parsign*diagthr:
                if centre=="diag":
                    tab = Tools.maskrad(lons,lats,ldiag[0],ldiag[1],diagrad)
                elif centre=="obj":
                    tab = Tools.maskrad(lons,lats,self.lonc,self.latc,diagrad)
                else:
                    print("Wrong value of centre in OBJ_Centre0D.mask ; should be obj or diag - ABORT")
                    exit()

        return tab

            ## Some routines for object detection around a guess

    def search_core(self,fic,inst,param,Hn,indf,algo,domtraj,res,parfilt,filtapply,track_parameter,basetime,subnproc,ldiag,**kwargs):
        #Search core (local extremum) in field at a distance below dist_max and,
        #if successful, creates the associated object
        #Returns None if no valid detection
        #Parameters:
          #self: input object around which to search
          #fic: file where to read parameter
          #inst: instant
          #param:field parameter
          #Hn:Hemisphere
          #indf,algo:inputdef,algodef
          #domtraj:domain
          #res: resolutions
          #parfilt,filtapply: filtering parameters
          #track_parameter
          #basetime
          #subnproc
          #ldiag: if not [], computes the list of diagnostics in ldiag
          #kwargs: criterions for detection
                #max_dist:maximal distance where to search extremum
                #ss:size of the square (degrees)
                #thr_param:threshold on the value
                #pairing : if this is a pairing parameter or not
                #smooth: if smoothing is applied or not

        #initialisation
        trackpar, signtrack = Tools.get_parsign(param, Hn)
        if filtapply==0:
            filtrad=0.0
        else:
            filtrad=parfilt[trackpar]*filtapply
        chk_dist=False
        ss=1000
        thr_track=Tools.maxval
        max_dist=Tools.maxval
        if "max_dist" in kwargs:
            check_dist=True
            max_dist=kwargs["max_dist"]
        if "ss" in kwargs:
            ss=kwargs["ss"]
        if "thr_param" in kwargs:
            thr_track=kwargs["thr_param"]
        smooth=False
        pairing=False
        if "smooth" in kwargs:
            smooth=kwargs["smooth"]
        if "pairing" in kwargs:
            pairing=kwargs["pairing"]

        fld = Inputs.extract_data(fic,inst,indf,trackpar,domtraj,res[trackpar],basetime,subnproc,filtrad=filtrad)
        lon,lat,val,dist,gook = Tools.find_locextr(signtrack,fld,self.lonc,self.latc,ss=ss,thr=thr_track)

        if gook and ((not chk_dist) or (dist<=max_dist)):

            if smooth:
                lon, lat, val = Tools.smooth_extr(signtrack,fld,lon,lat)
            #CREATE OBJECT
            objectm = ObjectM([], track_parameter,lonc=lon,latc=lat,time=inst)
            if trackpar in track_parameter:
                objectm.traps[trackpar] = val

            #DIAGNOSTIC PARAMETERS
                #if ldiag is [], Tools.make_diag gets out quickly
            if pairing:
                Tools.make_diags(ldiag,objectm,fic,inst,indf,domtraj,Hn,res,basetime,self.lonc,self.latc,subnproc,parfilt=parfilt,filtapply=filtapply)
            else:
                Tools.make_diags(ldiag,objectm,fic,inst,indf,domtraj,Hn,res,basetime,lon,lat,subnproc,parfilt=parfilt,filtapply=filtapply)

        else:
            objectm=None

        return objectm

    def conditiontype(self,fic,inst,indf,algo,domtraj,res,parfilt,filtapply,basetime,subnproc,init=False):
        #Determines if the object satisfies a condition type
        #init : if it is the first obj in the track
        #Output : isok (True is the object is of the type)
        #(in that case nameobj takes algo.conditiontype.nameobj value)
        #Output : exclude, if the object has to be removed from the track

        isok=True
        exclude=False
        if algo.conditiontype is not None:
            if "casetype" in algo.conditiontype:
                if algo.conditiontype["casetype"]=="TCVitart97":
                    if (not algo.conditiontype["when"]=="init") or init==True:
                        #print("Test condition "+algo.conditiontype["casetype"])
                        #print("... a tropical cyclone - "+algo.conditiontype["nameobj"]+" - criterion based on t200-t300-t400-t500, z200 and z1000 fields")

                        #Pameters
                        if abs(domtraj["lonmax"]-self.lonc)>3.0:
                            distmax=Tools.comp_length(self.lonc,self.latc,self.lonc+2.0,self.latc) #Distance (en km) de 2deg de longtitude
                        else:
                            distmax=Tools.comp_length(self.lonc,self.latc,self.lonc-2.0,self.latc) #Distance (en km) de 2deg de longtitude
                        dt = -0.5 #derivative criterion on temperature
                        dz = -500 #derivative criterion on thickness (500 mgp)

                        #Extraction of fields
                        t200 = Inputs.extract_data(fic,inst,indf,"t200",domtraj,res["t200"],basetime,subnproc,filtrad=parfilt["t200"]*filtapply)
                        t300 = Inputs.extract_data(fic,inst,indf,"t300",domtraj,res["t300"],basetime,subnproc,filtrad=parfilt["t300"]*filtapply)
                        t400 = Inputs.extract_data(fic,inst,indf,"t400",domtraj,res["t400"],basetime,subnproc,filtrad=parfilt["t400"]*filtapply)
                        t500 = Inputs.extract_data(fic,inst,indf,"t500",domtraj,res["t500"],basetime,subnproc,filtrad=parfilt["t500"]*filtapply)

                        #Computation of warm core criterion
                        tmean=t200
                        tmean.operation('+',t300)
                        tmean.operation('+',t400)
                        tmean.operation('+',t500)
                        tmean.operation("/",4)

                        lon,lat,val,dist,gook = Tools.find_locextr(1,tmean,self.lonc,self.latc,ss=4,thr=0.0)
                        if gook and (dist<distmax):
                            #Derivative criterion
                            deriv = Tools.comp_deriv(tmean,self.lonc,self.latc,8,domtraj)
                            if deriv[0]<dt and deriv[1]<dt and deriv[2]<dt and deriv[3]<dt:
                                isok=True
                            else:
                                isok=False
                        else:
                            isok=False

                        #Computation of geopotential criterion (if isok)
                        if isok:
                            z200 = Inputs.extract_data(fic,inst,indf,"z200",domtraj,res["z200"],basetime,subnproc,filtrad=parfilt["z200"]*filtapply)
                            z1000 = Inputs.extract_data(fic,inst,indf,"z1000",domtraj,res["z1000"],basetime,subnproc,filtrad=parfilt["z1000"]*filtapply)
                            z200.operation('-',z1000) #thickness
                            lon,lat,val,dist,gook = Tools.find_locextr(1,z200,self.lonc,self.latc,ss=4,thr=0.0)
                            if gook and (dist<distmax):
                                #Derivative criterion
                                deriv = Tools.comp_deriv(z200,self.lonc,self.latc,8,domtraj)
                                if deriv[0]<dz and deriv[1]<dz and deriv[2]<dz and deriv[3]<dz:
                                    isok=True
                                else:
                                    isok=False
                            else:
                                isok=False
                    #Consequences of isok: exclude, self.nameobj
                    if isok:
                        self.nameobj=algo.conditiontype["nameobj"]
                    if algo.conditiontype["exclude"]=="True" and not isok:
                        exclude=True
                else:
                    print(algo.conditiontype["casetype"] + " is not a valid condition type - Not applied")
            else:
                print("condition type does not apply - declare casetype in algo")
                exit()

        return isok, exclude

            ## Classe de la trajectoire de l'objet météorologique ##

''' La trajectoire est définie par une liste d'objet. Ces objets correspondent
aux objets météorologiques crées tout au long de l'exécution de l'algorithme,
soit à chaque pas de temps. Ceci permet de visualiser la trajectoire d'un objet
météorologique.
Elle prend en compte les coordonnées du domaine d'étude, utiles pour renvoyer
une carte dont le domaine est celui choisi. '''


class Track:

    def __init__(self,name="",**kwargs):
        #A track is a list of objects
        
        # Definition
        self.classobj = "Centre0D"
        self.nameobj = "" #Type of object : cyclone or ...
        self.name = name #Name of the track (for instance, name of the cyclone if it has been declared by RSMC)
        self.nobj = 0
        self.traj = []
        if "algodef" in kwargs.keys():
            self.algodef = kwargs['algodef'] #Only for model tracks 
        if "inputdef" in kwargs.keys():
            self.inputdef = kwargs['inputdef'] #Only for model tracks 
        if "basetime" in kwargs.keys():
            bt=kwargs['basetime'] #Only for model tracks
            if isinstance(bt,str):
                self.basetime=bt
            else:
                self.basetime=datetime.strftime(bt,time_fmt)
        if "nameobj" in kwargs.keys():
            self.nameobj=kwargs['nameobj'] #Only for model tracks

    def add_obj(self,objectm):
        #Add objectm in the track
	#No sorting is done (by time or other)
	#The name of the track does not change

        self.traj.append(objectm)
        self.nobj = self.nobj + 1

    def find_inst(self,ftime):
	#finds the object in the track that corresponds to the instant ftime
        #ftime has a datetime format or a "YYYYMMDDHH" string format

        objectm=[]
        found=False
        for i in range(self.nobj):
            otime=self.traj[i].time
            if isinstance(ftime,str):
                ftime2=ftime
            else:
                ftime2=datetime.strftime(ftime,time_fmt)

            if otime==ftime2:
                found=True
                objectm.append(self.traj[i])

        return objectm, found

    def find_ind(self,ftime):
	#finds the index of the object in the track that corresponds to the instant ftime
        #ftime has a datetime format or a "YYYYMMDDHH" string format

        indm=[]
        found=False
        for i in range(self.nobj):
            otime=self.traj[i].time
            if isinstance(ftime,str):
                ftime2=ftime
            else:
                ftime2=datetime.strftime(ftime,time_fmt)

            if otime==ftime2:
                found=True
                indm.append(i)

        return indm, found

    def dt_traj(self,unit="h"):
        #Compute time step between two objects in the track
        #If the track is a single point, then the output is 0

        if self.nobj>1:
            diffs=(datetime.strptime(self.traj[1].time,time_fmt)-datetime.strptime(self.traj[0].time,time_fmt)).seconds
            diffd=(datetime.strptime(self.traj[1].time,time_fmt)-datetime.strptime(self.traj[0].time,time_fmt)).days

            if unit.lower()=="h":
                td=diffs/3600.0+diffd*24.0
            elif unit.lower()=="d":
                td=diffs/3600.0/24.0+diffd
            elif unit.lower()=="s":
                td=diffs+diffd*3600.0*24.0
            else:
                print("wrong unit in OBJ_CYC traj tlen routin")
                exit()
        else:
            td=0

        return td

    def tlen(self,unit="h"):
        #Compute length (in time unit) of the track

        diffs=(datetime.strptime(self.traj[-1].time,time_fmt)-datetime.strptime(self.traj[0].time,time_fmt)).seconds
        diffd=(datetime.strptime(self.traj[-1].time,time_fmt)-datetime.strptime(self.traj[0].time,time_fmt)).days

        if unit.lower()=="h":
            tl=diffs/3600.0+diffd*24.0
        elif unit.lower()=="d":
            tl=diffs/3600.0/24.0+diffd
        elif unit.lower()=="s":
            tl=diffs+diffd*3600.0*24.0
        else:
            print("wrong unit in OBJ_CYC traj tlen routin")
            exit()


        return tl

    def write(self, write_diag=True):
        '''Write the track in a dictionary
           write_diag: True if the diagnostic values should be stored
'''

        # Initialisation of the dictionary that will be dumped in the json file
        dataout={}

        # Ecriture dans le fichier
        for i in range(self.nobj):
            
            # Récuperation de l'objet
            obj = self.traj[i]
            #print(obj.__dict__)
            t = obj.time
            
            # Récuperation des attributs de définition
            if obj.nameobj=="":
                dico = {'lonc': obj.lonc, 'latc': obj.latc}
            else:
                dico = {'lonc': obj.lonc, 'latc': obj.latc, 'nameobj':obj.nameobj}
            if "diags" in vars(obj) and write_diag:
                dico["diags"]=obj.diags
                for diag in obj.diags:
                    dico[diag]=obj.__dict__[diag]

            # Alimentation du dictionnaire par les valeurs souhaitées
            dataout["time:" + str(t)] = dico

        return dataout


    def read(self,dictin):

        '''Reads the track from input dictionary (including diags)
        '''

        for dico in dictin.items():
            diconame = dico[0]
            if diconame[0:4]=="time":
                t0 = diconame[5:]
                lat0 = dico[1]['latc']
                lon0 = dico[1]['lonc']
                ldiag = dico[1]
                del ldiag['latc']
                del ldiag['lonc']
                obj = ObjectM(ldiag,[],lonc=lon0,latc=lat0,time=t0)
                if "nameobj" in dico:
                    obj.nameobj=dico[1]['nameobj']
                if len(ldiag)>0:
                    obj.diags=ldiag
                    for diag in ldiag:
                        setattr(obj,diag,dico[1][diag])
                self.add_obj(obj)


    def plot(self, ax):

        '''visualisation de la trajectoire'''

        plt.text(self.traj[0].lonc,self.traj[0].latc+1.0,self.name)
        for i in range(len(self.traj)-1):

            lat_i = self.traj[i].latc
            long_i = self.traj[i].lonc
            lat_i1 = self.traj[i+1].latc
            long_i1 = self.traj[i+1].lonc
            ax.plot(long_i, lat_i, '.')
            ax.plot(long_i1, lat_i1, '.')
            plt.plot([long_i, long_i1], [lat_i, lat_i1], color = 'red',
                     linestyle = 'solid')

    def tdist(self, tra):
        #computes the list of t-distances the objects in each track
        #Output : list of distance of objects that correspond to the same instants

        tld=[]

        for obj1 in self.traj:
            for obj2 in tra.traj:
                td=obj1.tdist(obj2)
                if td<Tools.maxval/2.0:
                    tld.append(td)

        return tld

    def tmax(self,dia):
        #Computes the maximum value of the diagnostic along the track

        val = missval #assumes missval is negative
        for obj in self.traj:
            vect=obj.__dict__[dia]
            nd=int(len(vect)/3)
            for ivi in range(nd):
                if vect[3*ivi+2] > val:
                    val = vect[3*ivi+2] 

        return val

    def tmin(self,dia):
        #Computes the minimum value of the diagnostic along the track

        val = 1e10 #assumes missval is negative
        for obj in self.traj:
            vect=obj.__dict__[dia]
            nd=int(len(vect)/3)
            for ivi in range(nd):
                if (vect[3*ivi+2]>missval) and vect[3*ivi+2] < val:
                    val = vect[3*ivi+2] 

        return val

#Some routines
def search_allcores(fic,inst,indf,trackpar,signtrack,domtraj,res,basetime,track_parameter,subnproc,filtrad=0.0,dist=Tools.maxval, thr=Tools.maxval):
    #Search for all objects in a domain
    #that are at a minimal distance of radmax one to the other
    #and above thr_core
    #Output : list of objects

    fld = Inputs.extract_data(fic,inst,indf,trackpar,domtraj,res,basetime,subnproc,filtrad=filtrad)
    tlon, tlat, tval = Tools.find_allextr(signtrack, fld, dist=dist, thr=thr)
    lobj=[]
    for ivi in range(len(tlat)):
        objectm = ObjectM([], track_parameter,lonc=tlon[ivi],latc=tlat[ivi],time=inst)
        if trackpar in track_parameter:
            objectm.traps[trackpar] = tval[ivi]
        lobj.append(objectm)

    return lobj


