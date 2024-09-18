# -*- coding: utf-8 -*-
"""
Traject main program

https://github.com/UMR-CNRM/Traject
Traject is disributed under Cecill-C license
(For rights and obligations, see LICENSE.txt)

UMR-CNRM, Meteo-France & CNRS

"""
import copy
import glob
import os, glob
import sys
import time as timen
from datetime import datetime,timedelta,time
import h5py
import json
import importlib
import numpy as np
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import Inputs, Tools
import epygram
import concurrent.futures

#Generic variables
traject_version = 0.92
missval = -9999.0
time_orig = np.datetime64("1970-01-01T00:00:00")
time_fmt = "%Y%m%d%H"

def DefObject(classobj, ldiag, ltrap, **kwargs):
    '''Function to import the class of meteorological object that is needed by the user
    and creates an object

    Inputs: classobj : object class
            **kwargs arguments are optional and depend on the object
    Ouputs : an object of ObjectM class.
    '''

    #L'importation du module convenable
    module = __import__("OBJ_"+classobj)

    #Définition des classes qu'on veut utiliser
    myobj_cls = getattr(module, "ObjectM") # This defines the class to be used

    return myobj_cls(ldiag,ltrap,**kwargs)

def DefTrack(classobj,**kwargs):
    '''Function to import the class of meteorological object that is needed by the user
    and creates a track

    Inputs: classobj : object class ;
    Ouputs : an object of ObjectM class, an objet of Track class.
    '''

    #L'importation du module convenable
    module = __import__("OBJ_"+classobj)

    #Définition des classes qu'on veut utiliser
    #myobj_cls = getattr(module, "ObjectM") # This defines the class to be used
    track_cls = getattr(module, "Track")

    return track_cls(**kwargs)

def Search_allcores(classobj,**kwargs):

    #L'importation du module convenable
    module = __import__("OBJ_"+classobj)

    search_allcores = getattr(module,"search_allcores")

    return search_allcores(**kwargs)

def PrepareFiles(algo,indf,dirkeys,timetraj,**kwargs):
    #Same input arguments as track (algo, indf, timetraj + termtraj if fc)
    #dirkeys: {"vortex"="", --> where data should be downloaded through vortex first (=$MTOOLDIR)
               # "extract"="", --> if the selection of fields should apply
               # "filt"="",  --> if subdomain and filter should apply
               # "dirout"="", --> where to write the resulting inputdef files (optional)}
    #If a dirkey is not specified, the action (vortex/extract/filt...) does not apply
    #As optional kwargs, we can add a loop on members (given as lists)
    #Output : None - the inputdef configuration objects are written for each selection, though

    #Processing inputdef
    if isinstance(indf,str):
        #indf is a string: we read the file
        indf2=Inputs.inputdef()
        indf2.read_input(indf)
    else:
        #indf is already a dictionary
        indf2=indf
    indf0=copy.deepcopy(indf2)

    #Processing algodef
    if isinstance(algo,str):
        #algo is a string: we read the file
        algo2=Inputs.algodef()
        algo2.read_input(algo)
    else:
        #algo is already a dictionary
        algo2=algo
    lvar=list(algo2.parfilt.keys())
    
    #Special case "rr" and decumulation
    rrdecum, rrbefore = Inputs.ifdecum(indf2)
    isparr=False
    ivi=0
    lvar2=[]
    for var in lvar:
        ivi=ivi+1
        if var[0:2]=="rr" and rrdecum:
            lvar2.append("rr")
        else:
            lvar2.append(var)
    lvar = [*set(lvar2)] #remove duplicates

    #subnproc
    subnproc=1
    if algo2.parallel is not None:
        if "subnproc" in algo2.parallel:
            subnproc=algo2.parallel["subnproc"]

    if "members" in kwargs and "member" in indf2.__dict__.keys():
        lmb=kwargs['members'] #should be a list of integers
        nmb=len(lmb)
    else:
        lmb=[]
        nmb=1

    #Loop on members
    for imb in range(nmb):

        #initializations
        lfile=[]
        linst=[]
        lbase=[]
        indf2=copy.deepcopy(indf0)
        
        ###
        # - READ INPUTFILE DEFINITION
        # - PERFORM VORTEX DOWNLOAD IF REQUIRED
        ###
        

        if len(lmb)>0: #Several members
            mbs=str(lmb[imb]).rjust(3,'0')

        if indf2.origin=="fc": #Forecast data is processed: algo is applied on term steps
            termtraj=kwargs["termtraj"]
        
            for basetime in Inputs.comptimes(timetraj): #Loop on timetraj
                #print(basetime)
                if len(lmb)>0: #Several members
                    indf2.member=mbs
                if "vortex" in dirkeys.keys():
                    indf2.get_vortex(termtraj,basetime)
                    if imb==0 and "dirout" in dirkeys.keys():
                        indf2.write_input(dirkeys["dirout"]+"indef_vortex.json")
                lfileb,linstb = indf2.get_filinst(termtraj,basetime)
                lbaseb=[basetime for ivi in range(len(lfileb))]
                lfile.extend(lfileb)
                linst.extend(linstb)
                lbase.extend(lbaseb)

        elif indf2.origin=="an" or self.origin=="cl": #algo is applied on instants
            if "vortex" in dirkeys.keys():
                indf2.get_vortex(timetraj)
                if imb==0 and "dirout" in dirkeys.keys():
                    indf2.write_input(dirkeys["dirout"]+"indef_vortex.json")
            lfileb,linstb = indf2.get_filinst(timetraj)
            lbaseb=Inputs.comptimes(timetraj)
            lfile.extend(lfileb)
            linst.extend(linstb)
            lbase.extend(lbaseb)

        else:
            print("ABORT - indf.origin not defined - Should be fc, an or cl")
            exit()

        #Here, lfile is the list of files to process ... we launch the processing
        #(parallelisation is possible)

        ###
        # - PERFORM FIELD EXTRACTION IF REQUIRED
        ###
        if "extract" in dirkeys.keys():
            print("Performing field extraction ...")
            for ivi in range(len(lfile)):
                fic=lfile[ivi]
                print(fic)
                inst=linst[ivi]
                binst=lbase[ivi]
                term=linst[ivi]-datetime.strptime(lbase[ivi],time_fmt)
                bterm=str(int(24*term.days + round(term.seconds)/3600.)).rjust(4,'0')
                rep, fname, ext = Inputs.split_filename(fic)
                if len(lmb)>0: #Several members
                    fname2 = dirkeys["extract"]+'XTRACT'+binst+'_'+bterm+'_m'+mbs+'.'+ext
                elif indf2.origin=="fc":
                    fname2 = dirkeys["extract"]+'XTRACT'+binst+'_'+bterm+'.'+ext
                else:
                    fname2 = dirkeys["extract"]+'XTRACT'+binst+'.'+ext
                #print("Extract " + fic + " to " + fname2)
                #Extraction parameter from the grib
                Inputs.extract_field(fic,inst,indf2,lvar,fname2)
                lfile[ivi]=fname2

            if imb==0: #first member loop only
                indf2.directory=dirkeys["extract"]
                if len(lmb)>0: #Several members
                    indf2.filename='XTRACT[YYYYMMDDHH]_[term:4]_m[member].'+ext
                elif indf2.origin=="fc":
                    indf2.filename='XTRACT[YYYYMMDDHH]_[term:4].'+ext
                else:
                    indf2.filename='XTRACT[YYYYMMDDHH].'+ext
                if "dirout" in dirkeys.keys():
                    indf2.write_input(dirkeys["dirout"]+"indef_xtract.json")

        ###
        # - PERFORM SUBDOMAIN AND FILTERING IF REQUIRED
        ###
        if "filter" in dirkeys.keys():
            print("Performing filtering and subdomain extraction ...")

            if imb==0: #first member loop only
                #Define domains and resolutions
                if algo2.domtraj is None:
                    domtraj, res, lons, lats = Inputs.extract_domain(lfile[0],linst[0],indf,lvar[0]) #domtraj is the total grid
                else:
                    domtraj=algo2.domtraj
                    parres = Tools.get_parres(algo2,indf2,lfile[0],linst[0])
                    #print(parres)

            for ivi in range(len(lfile)):
                fic=lfile[ivi]
                print(fic)
                inst=linst[ivi]
                binst=lbase[ivi]       
                term=linst[ivi]-datetime.strptime(lbase[ivi],time_fmt)
                bterm=str(int(24*term.days + round(term.seconds)/3600.)).rjust(4,'0')
                rep, fname, ext = Inputs.split_filename(fic)
                if len(lmb)>0: #Several members
                    fname2 = binst+'_'+bterm+'_m'+mbs+'.nc'
                elif indf2.origin=="fc":
                    fname2 = binst+'_'+bterm+'.nc'
                else:
                    fname2 = binst+'.nc'
                #print("Process " + fic + " to " + dirkeys["filter"]+"FILT$res..."+fname2)
                param_file, param_nc = Inputs.filter_field(fic,inst,indf2,algo2.parfilt,domtraj,parres,dirkeys["filter"],binst,fname2,subnproc)

            if imb==0: #first member loop only
                indf2.directory=dirkeys["filter"]
                indf2.nativefmt="nc"
                indf2.update_input({"filtered":algo2.parfilt})
                indf2.domain=algo2.domtraj
                indf2.update_input({"param_nc":param_nc})
                if "rv_av" in indf2.special_keys:
                    indf2.special_keys.remove("rv_av") #filter_field traduit le av en rv ...
                if "rr_before" in indf2.special_keys:
                    indf2.special_keys.remove("rr_before") #decumulation has been done ...
                if "rr_after" in indf2.special_keys:
                    indf2.special_keys.remove("rr_after") #decumulation has been done ...
                indf2.update_input({"param_file":param_file})
                if len(lmb)>0: #Several members
                    indf2.filename='FILT[param_file]-[YYYYMMDDHH]_[term:4]_m[member].nc'
                elif indf2.origin=="fc":
                    indf2.filename='FILT[param_file]-[YYYYMMDDHH]_[term:4].nc'
                else:
                    indf2.filename='FILT[param_file]-[YYYYMMDDHH].nc'

                if "dirout" in dirkeys.keys():
                    indf2.write_input(dirkeys["dirout"]+"indef_filter.json")


def read_nml(fic):
    #reads namelist file fic (json format)
    #checks the type of namelist (by look at some keys inside)
    #and outputs the right type of namelist object (inputdef or algodef)

    if isinstance(fic,str):
    #fic is a filename
        with open(fic,'r') as json_file:
            fic2 = json.load(json_file)
    else:
    #fic should already be an object
        fic2 = fic.__dict__

    #Detects the type of namelist from the content of fic2 dictionary
    if "origin" in fic2.keys():
        #inputdef file
        objdef=Inputs.inputdef()
        objdef.update_input(fic2)
    elif "varalgo" in fic2.keys():
        #algodef file
        objdef=Inputs.algodef()
        objdef.update_input(fic2)
    else:
        print("Wrong "+ fic + " file to be read ... abort")
        exit()

    return objdef


def track(algo,indf,timetraj,**kwargs):
    #Definition and launch of the algorithm
    #traj = DefAlgo(algo,**kwargs)
    #Compulsory input parameters:
        #algo: algodef namelist (can be json file or a dictionary - see Inputs.py)
        #indf: definition of the inputdata (can be json file or a dictionary - see Inputs.py)
        #timetraj: time on which tracking should be ran (basetime for (re-)forecasts,
            #validtime for analyses or climate) - format: {'start':"2020021500",'end':"2020021500",'step':"06"} (step in hours)
        #--compulsory only for forecasts-- termtraj: forecast terms along which to track
            #- format: {'init':0,'final':45,'step':3}, init=0 by default
    #Optional input parameters:
        #outfile : output file to store the tracks
    #Output : list of tracks
    
    #Processing inputdef and algodef
    indf2=read_nml(indf)
    algo2=read_nml(algo)
    
    #Definition of the tracking algorithm
    module = __import__("ALGO_"+algo2.name)
    myalgo_func = getattr(module,"track")

    #Initialisation of trajlist
    outf=[]
    outw=[]
    trajlist=[]
    indfw=[]

    if "members" in kwargs and "member" in indf2.__dict__.keys():
        lmb=kwargs['members'] #should be a list of integers
        nmb=len(lmb)
    else:
        lmb=[]
        nmb=1
    print("List members :", lmb)

    #Read reftraj (if any)
    if "reftraj" in kwargs:
        refnam = kwargs["reftraj"]
        nref = 1
        kwargs.pop("reftraj")
    else:
        lref = ["noref"]
        nref = 0

    #Parallelization options (only applies on fc)
    write_parall = False
    if "ntasks" in algo2.parallel:
        ntasks = algo2.parallel["ntasks"]
        print("Compute ",ntasks," tracks in parallel.")
        #executor = concurrent.futures.ThreadPoolExecutor(max_workers=ntasks)
        executor = concurrent.futures.ProcessPoolExecutor(max_workers=ntasks)
        if "write_parall" in algo2.parallel:
            if algo2.parallel["write_parall"] == "True":
                print("Writing output files will be done in parallel")
                write_parall = True
                execwrite = concurrent.futures.ProcessPoolExecutor(max_workers=ntasks)
    else:
        ntasks = 1

    if indf2.origin=="fc": #Forecast data is processed: algo is applied on term steps
        tsk = 0

        for basetime in Inputs.comptimes(timetraj): #Loop on timetraj
            lfile,linst2 = indf2.get_filinst(kwargs["termtraj"],basetime)
            if nref>0:
                lref = Tools.get_reftraj(linst2,algo2.domtraj,refnam)

            for reftraj in lref: #Loop on lref
                #Loop on members
                for imb in range(nmb):

                    if ntasks==1: #No parallelisation
                        print("No Upper-level Parallelisation")
                        outf.append(track_parallel_fc(myalgo_func,0,algo2,indf2,basetime,lmb,imb,reftraj,**kwargs))
                    else: #Parallelisation
                        tsk = tsk + 1
                        outf.append(executor.submit(track_parallel_fc,myalgo_func,tsk,algo2,indf2,basetime,lmb,imb,reftraj,**kwargs))

        tsk = 0 
        for out in outf:
            if ntasks==1:
                outtraj, indf3 = out
            else:
                tsk = tsk + 1
                try:
                    outtraj, indf3 = out.result()
                except Exception as exc:
                    print("Tracking - Exception raised by Task "+str(tsk)+":",exc)
                    outtraj=[]
                    indf3=[]

            #Writing outputfile file in parallel if outfile is specified
            if "outfile" in kwargs and write_parall:
                outfile = '.'.join(kwargs["outfile"].split('.')[0:-1])+"_"+str(tsk)+"."+kwargs["outfile"].split('.')[-1]
                print("Write in parallel file ",outfile)
                indfw = []
                for ivi in range(len(outtraj)):
                    indfw.append(copy.deepcopy(indf3))
                outw.append(execwrite.submit(write_fc,outtraj,outfile,algo2,indfw))
                for out1 in outw:
                    try:
                        ot = out1.result()
                    except Exception as exc:
                        print("Writing - Exception raised by Task "+str(tsk)+":",exc)
            else:
                trajlist.extend(outtraj) #Full list of tracks
                for ivi in range(len(outtraj)):
                    indfw.append(copy.deepcopy(indf3))

    elif indf2.origin=="an" or self.origin=="cl": #algo is applied on instants
        #No parallelization
        lfile,linst = indf2.get_filinst(timetraj)
        outtraj = myalgo_func(algo2,indf2,linst,lfile,**kwargs)
        trajlist.extend(outtraj)
        for ivi in range(len(outtraj)):
            indfw.append(copy.deepcopy(indf2))

    else:
        print("ABORT - indf.origin not defined - Should be fc, an or cl")
        exit()

    #Writing outputfile file if outfile is specified
    if "outfile" in kwargs and not write_parall:
        print("Write all tracks in single file")
        outfile = kwargs["outfile"]
        write_fc(trajlist,outfile,algo2,indfw)

    if "plotfile" in kwargs and len(trajlist)>0:
        if algo2.domtraj is None:
            varname=list(algo.parfilt.keys())
            domtraj, lons, lats = Inputs.extract_domain(lfile[0],linst[0],indf2,varname[0]) #domtraj is the total grid
        else:
            domtraj=algo2.domtraj

        if indf2.origin=="fc":
            rep, fn , ext = Inputs.split_filename(kwargs["plotfile"])
            plotfile = rep+fn
            print("Plot all tracks in "+ plotfile +"_[basetime]." + ext + " (one file per basetime)")
            for basetime in Inputs.comptimes(timetraj): #Loop on timetraj
                trajb=[tra for tra in trajlist if tra.basetime==basetime]
                Plot(trajb,domtraj,outputfile=plotfile+"_"+basetime+'.'+ext)
        else:
            plotfile = kwargs["plotfile"]
            print("Plot all tracks in "+ plotfile)
            Plot(trajlist,domtraj,outputfile=plotfile)

    return trajlist

def track_parallel_fc(func,b,algo,indf,basetime,lmb,imb,reftraj,**kwargs):

    print("Start PROCESS ",b," :",datetime.now().strftime("%H:%M:%S"))
    
    indf4=copy.deepcopy(indf)
    if len(lmb)>0: #Several members
        mbs=str(lmb[imb]).rjust(3,'0')
        indf4.member=mbs
        print("Tracking - member="+ str(lmb[imb]) + ", basetime="+basetime)

    lfile,linst = indf4.get_filinst(kwargs["termtraj"],basetime)
    outtraj = func(algo,indf4,linst,lfile,basetime=basetime,reftraj=[reftraj],**kwargs)

    print("End PROCESS ",b," :",datetime.now().strftime("%H:%M:%S")," - ",len(outtraj)," tracks have been found")

    return outtraj, indf4

def write_fc(ltraj,ofile,algo,lindf):

    if len(ltraj)>0:
        if os.path.exists(ofile):
            os.remove(ofile)
        for i in range(len(ltraj)):
            Write(ltraj[i],ofile,algo,lindf[i],mode="a")
    else:
        Write([],ofile)

def Plot(ltraj,dom,outputfile="out.png"):
    #Plots on a single figure the tracks that are
    #in the list ltraj in outputfile (default=out.png)
    #if a trajectory has a single object, then is it plotted as an object 

    #Object module call
    if len(ltraj)>0:
        __import__("OBJ_"+ltraj[0].classobj)

    #Make figure
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())
    ax.set_extent([dom["lonmin"],dom["lonmax"],dom["latmin"],dom["latmax"]])

    ax.add_feature(cfeature.OCEAN.with_scale('50m'))
    ax.add_feature(cfeature.COASTLINE.with_scale('50m'))
    ax.add_feature(cfeature.RIVERS.with_scale('50m'))
    ax.add_feature(cfeature.BORDERS.with_scale('50m'), linestyle=':')

    if len(ltraj)>0:
        for i in range(len(ltraj)):
            if ltraj[i].nobj>1:
                ltraj[i].plot(ax)
            elif ltraj[i].nobj==1:
                ltraj[i].traj[0].plot(ax)

    #Storage of the plot
    plt.savefig(outputfile)

def write_header(traj,algo,indf):
    #Input : a single trajectory and inputdef and algodef objects
    #If they are not given (ie : {}), then the information is taken from traj (if available)
    #Ouput : a header (dict format)

    adict={}
    
    if isinstance(indf,str):
        #indf is a string: we read the file
        indf2=Inputs.inputdef()
        indf2.read_input(indf)
    elif isinstance(indf,dict):
        #indf is a dictionary
        if indf=={} and "inputdef" in traj.__dict__:
            indf2=Inputs.inputdef()
            indf2.update_input(traj.inputdef)
        else:
            indf2=Inputs.inputdef()
            indf2.update_input(indf)
    else: #indf is an object
        indf2 = indf

    if isinstance(algo,str):
        #algo is a string: we read the file
        algo2=Inputs.algodef()
        algo2.read_input(algo)
    elif isinstance(algo,dict):
        #algo is a dictionary
        if algo=={} and "algodef" in traj.__dict__:
            algo2=Inputs.algodef()
            algo2.update_input(traj.algodef)
        else:
            algo2=Inputs.algodef()
            algo2.update_input(algo)
    else: #algo is an object
        algo2 = algo

    if not indf2=={}:
        if indf2.origin=="obs":
            idict = {"origin":"obs",
                "database":indf2.database,
                "domain":indf2.domain}
            if "nameobj" in traj.__dict__.keys():
                namobj=traj.nameobj
            else:
                namobj=traj.classobj
            tdict = {"classobj":traj.classobj,
                 "nameobj":namobj,  
                "name":traj.name,
                "basetime":traj.basetime,
                "nobj":traj.nobj}
        else:
            idict = {"origin":indf2.origin,
                "model": indf2.model,
                "domain": indf2.domain,
                "experiment": indf2.experiment}
            if "cutoff" in indf2.__dict__.keys():
                idict.update({"cutoff": indf2.cutoff})
            if "member" in indf2.__dict__.keys():
                idict.update({"member": indf2.member})
            if "scenario" in indf2.__dict__.keys():
                idict.update({"scenario": indf2.scenario})
            if not algo2=={}:
                adict = {"name":algo2.name,
                     "classobj":algo2.classobj,
                     "nameobj":algo2.nameobj,
                     "varalgo":algo2.varalgo,
                     "parfilt":algo2.parfilt,
                     "parres":algo2.parres,
                     "specfields":algo2.specfields,
                     "domtraj":algo2.domtraj}
            else:
                if "algodef" in traj.__dict__:
                    adict=traj.algodef
                else:
                    print("Warning : no algodef information could be found to write ...")
                    adict={}
            tdict = {"name":traj.name,
                "basetime":traj.basetime,
                "nobj":traj.nobj}
    else:
        if "inputdef" in traj.__dict__:
            idict=traj.inputdef
        else:
            print("Warning : no inputdef information could be found to write ...")
            idict={}

    return {"inputdef":idict, "algodef": adict, "trajdef":tdict}

def read_header(hdict):
    #Input: dictionary read from the json file
    #Output: empty track object with the header 

    #print('--------- keys ---------')
    #print(hdict.keys())
    
    if not "algodef" in hdict.keys(): #Observation
        if "basetime" in hdict.keys(): #version<=0.51
            bt=hdict["basetime"]
        elif "basetime" in hdict["trajdef"]: #version>0.51
            bt=hdict["trajdef"]["basetime"]
        else:
            bt=""
        if "inputdef" in hdict.keys():
            indf = hdict["inputdef"]
        else:
            indf={}
        if "nameobj" in hdict["trajdef"]:
            otraj = DefTrack(hdict["trajdef"]["classobj"],name=hdict["trajdef"]["name"],nameobj=hdict["trajdef"]["nameobj"],inputdef=indf, \
                        basetime=bt,nobj=0)
        else:
            otraj = DefTrack(hdict["trajdef"]["classobj"],name=hdict["trajdef"]["name"],basetime=bt,nobj=0)
    else:
        otraj = DefTrack(hdict["algodef"]["classobj"],name=hdict["trajdef"]["name"],algodef=hdict["algodef"],\
             trajdef=hdict["trajdef"],inputdef=hdict["inputdef"],basetime=hdict["trajdef"]["basetime"],nobj=0)

    return otraj

def Write(traj,outputfile,algo={},indf={},write_diag=True,mode='a'):
    #Write (append or write if mode='w') a single track traj in outputfile
    #Some indf and algo information can also be written in the file
    #Output: file (JSON format)

    if isinstance(traj,list):
        if len(traj)==0:
            with open(outputfile, 'w') as fout:
                json.dump([],fout)
        else:
            print("Write only empty list or single track - ABORT")
            exit()
    else:
        #Create data in a dictionary format
        __import__("OBJ_"+traj.classobj)
        wdict = write_header(traj,algo,indf)
        if wdict["algodef"]=={}:
            wdict.pop("algodef")
        wdict.update({"traj":traj.write(write_diag)})

        result = list()

        if mode=='w' or not os.path.exists(outputfile): #A new data is created
            result=[wdict]
        else: #json data is merged (as a list)
            with open(outputfile, 'r') as infile:
                result.extend(json.load(infile))
            result.extend([wdict])

        with open(outputfile, 'w') as fout:
            json.dump(result, fout, indent=4)


def Read(inputfile,select=["all"]):
    ''' Reads inputfile and generates a list of track
     according to the content of the file
     select (optional) is a dictionary that applies the Select() routine on the list of tracks
     (exemple : out=Read("in.json",select={"member":000})
     '''

    lname=inputfile.split('.')
    filefmt=lname[-1]
    filename=inputfile.replace('.'+filefmt,'')
    if filefmt.lower()=="json":
        json_file = open(inputfile)
        datain = json.load(json_file)

        ltraj=list()
        for i in range(len(datain)):
            traj=read_header(datain[i])
            __import__("OBJ_"+traj.classobj)
            traj.read(datain[i]["traj"])
            ltraj.append(traj)

    ltraj = Select(ltraj, select)

    return ltraj

def Select(ltraj,select):
    #selects the list of tracks the ones that correspond to the values given in the select list:
    # - a member ({"member"=}),
    # - a basetime ({"basetime"=...}),
    # - a track name ({"name"=...}),
    # - an interval of instants ({"time"={}), which are the valid times - in that case the tracks can be truncated
    # ... more possible to come

    #Apply select
    ltrajout=ltraj
    if select==["all"]:
        ltrajout=ltraj
    elif select==[]:
        ltrajout=[]
    else:
        for ff in select:

            if ff=="member":

                if isinstance(select[ff],str):
                    mbs=select[ff].rjust(3,"0")
                else:
                    mbs=str(select[ff]).rjust(3,"0")
                ltrajout=[traj for traj in ltrajout if traj.inputdef["member"]==mbs]

            elif ff=="basetime":
                bs=select[ff]
                ltrajout=[traj for traj in ltrajout if traj.basetime==bs]

            elif ff=="name":
                ns=select[ff]
                ltrajout=[traj for traj in ltrajout if traj.name==ns]

            elif ff=="time":
                ltrajout2=[]
                timetraj=select[ff]
                if not "step" in timetraj.keys():
                    timetraj["step"]=1
                linst=Inputs.comptimes(timetraj)
                for traj in ltrajout:
                    trajout2 = copy.deepcopy(traj)
                    trajout2.nobj = 0
                    trajout2.traj=[]
                    for obj in traj.traj:
                        if obj.time in linst:
                            trajout2.add_obj(obj)
                    if trajout2.nobj>0:
                        ltrajout2.append(trajout2)
                ltrajout = ltrajout2
            else:
                print("Wrong select key in traject.Select()")
    
    return ltrajout

def Clean_should(rep=".",force=False):
    #Clean the should* files generated by vortex in rep

    lfiles=glob.glob(rep+"/shouldfly-*")
    print("#############################")
    print("Directory " + rep + " - List of should-fly files to be deleted ... ")
    for name in lfiles:
        print(name)
    print("#############################")
    if force:
        print("Clean the shouldfly-* files ...")
        for name in lfiles:
            print("Removing " + name + " ...")
            os.remove(name)
    else:
        print("Clean the shouldfly-* files?... if ok press y ")
        k=input()
        if k=='y' or k=='Y':
            print("Removing all shouldfly- files...")
            for name in lfiles:
                print("Removing " + name + " ...")
                os.remove(name)

def Clean_inputfiles(indf,timetraj,**kwargs):
    #Removes all the files corresponding to inputdef indf and timetraj
    #In kwargs, we can declare termtraj (needed for forecasts)
        #and use Force=True

    if isinstance(indf,str):
        #indf is a string: we read the file
        indf2=Inputs.inputdef()
        indf2.read_input(indf)
    else:
        #indf is already a dictionary
        indf2=indf

    if "members" in kwargs and "member" in indf2.__dict__.keys():
        lmb=kwargs['members'] #should be a list of integers
        nmb=len(lmb)
    else:
        lmb=[]
        nmb=1

    force=False
    if "Force" in kwargs.keys():
        force=kwargs["Force"]

    #Loop on members
    lfiles=[]
    for imb in range(nmb):

        if len(lmb)>0: #Several members
            mbs=str(lmb[imb]).rjust(3,'0')

        if indf2.origin=="fc": #Forecast data is processed: algo is applied on term steps
            termtraj=kwargs["termtraj"]
        
            for basetime in Inputs.comptimes(timetraj): #Loop on timetraj
                if len(lmb)>0: #Several members
                    indf2.member=mbs
                lfile,linst = indf2.get_filinst(termtraj,basetime)
                lfiles.extend([name for name in lfile if os.path.exists(name)])

        elif indf2.origin=="an" or self.origin=="cl": #algo is applied on instants
            lfile,linst = indf2.get_filinst(timetraj)
            lfiles.extend([name for name in lfile if os.path.exists(name)])

    print("#############################")
    print(" List of files to be deleted ... ")
    for name in lfiles:
        print(name)
    print("#############################")
    if force:
        print("Clean the files ...")
        for name in lfiles:
            print("Removing " + name + " ...")
            os.remove(name)
        print("done")
    else:
        print("Clean the files?... if ok press y ")
        k=input()
        if k=='y' or k=='Y':
            print("Removing all files...")
            for name in lfiles:
                print("Removing " + name + " ...")
                os.remove(name)
        print("done")

    print("#############################")
    print(" Delete empty directories ... ")
    for ivi in range(4): #recursive loop
        lrep=[]
        print("Recursive deletion - Step ",ivi)
        for name in lfiles:
            ln=name.split('/')[0:-1]
            rep=''
            for sn in ln:
                if len(sn)>0:
                    rep=rep+'/'+sn
            #print(rep)
            if rep not in lrep:
                lrep.append(rep)

        for rep in lrep:
            for root,lrep2,files in os.walk(rep):
                for rep2 in lrep2:
                    rep3=os.path.join(root,rep2)
                    if len(os.listdir(rep3)) == 0:
                        print("Delete path ",rep3)
                        os.rmdir(rep3)
        lfiles=lrep

    print("#############################")

    return

def match_tracks(ltraj, tmax, dist, mininst, minmb, prefname):
    #Creates groups of tracks in ltraj that match (distance below dist for at least mininst instants)
    #Groups must contain at least minnb tracks
    #Two tracks from the same member and the same basetime cannot belong to the same group
    #Inputs: ltraj list of tracks
        #dist : maximum radius (km) of the circle that contains the tracks of a group
        #mininst : number of instants for which tracks must match
        #minmb : minimum number of tracks in a group
        #tmax : maximum instants from the beginning the tracks to be taken into account for matching (not applied if zero)
        #prefname : prefix name that will be used to name the tracks
    #Outputs: list of tracks, tracks that match have the same name (prefname-*)

    #Initializations
    nn=Tools.maxval
    if tmax==0:
        tmax2 = 10000 #extreme max value
    else:
        tmax2 = tmax
    numane=0
    lname=[]
    ltrajout = []
    ltraj2 = ltraj
    gook=True

    #Introduire while ici
    while gook:

        ntra=len(ltraj2)
        print("Number of tracks:", len(ltrajout))
        print("Number of tracks remaining:", ntra)
        mtc= np.zeros(shape=(ntra,ntra)) #array to check if tracks match
        tdist = np.zeros(shape=(ntra,ntra)) #Mean interdistance between two tracks
        tn = np.zeros(shape=(ntra,ntra)) #Number of common instances between two tracks
        totdist = [Tools.maxval for ivi in range(ntra)] #total distance of a track with its neighbours

        #Add member if not defined
        for tk in ltraj2:
            if "member" not in tk.inputdef:
                tk.inputdef["member"]=""
            
        idcd= [tk.basetime+tk.inputdef["member"] for tk in ltraj2] #identifier of each trajectory

        #
        print(" ----- Start match ------")
        
        ik=-1
        for tk in ltraj2:
            ik=ik+1
            #print("Track no", ik," : ", tk.__dict__)
            ij=-1
            for tj in ltraj2:
                ij=ij+1
                if not idcd[ik] == idcd[ij]: #Different basetime or member
                    dtkj=tk.tdist(tj,tmax=tmax2)
                    #print(ik, ij, dtkj)
                    dist2=[dtkj[ivi] for ivi in range(len(dtkj)) if dtkj[ivi]<=dist]
                    #print(dist2)
                    if len(dist2)>=mininst:
                        mtc[ik,ij]=1
                        tn[ik,ij]=len(dist2)
                        tdist[ik,ij]=np.mean(dist2)
                        #print("match!" + idcd[ik] + " - " + idcd[ij])
            #TEST TO FILTER OUT DOUBLE VALUES
            for ij in range(ntra):
                if mtc[ik,ij]==1:
                    sameij = [ii for ii in range(ntra) if (mtc[ik,ii]==1 and idcd[ii]==idcd[ij])]
                    #print("SAMEIJ:",sameij)
                    if len(sameij)>1: #there are multiple values: we keep the closest one to tk
                        jmax=sameij[0]
                        for ii in sameij:
                            if tn[ik,ii]>tn[ik,jmax]:
                                jmax=ii
                            elif tn[ik,ii]==tn[ik,jmax]:
                                if tdist[ik,ii]<tdist[ik,jmax]:
                                    jmax=ii
                        for ii in sameij:
                           if not ii==jmax:
                               #print("REMOVE ",ii)
                               mtc[ik,ii]=0

        #Checks the tracks that has the highest number of neighbours
        nn=0
        for ik in range(ntra):
            nk = len([ivi for ivi in range(ntra) if mtc[ik,ivi]==1])
            if nk>=nn:
                nn=nk
        #print("maximum number of neigbours found: ",nn)

        kmax=[]
        ik=-1
        for tk in ltraj2:
            ik=ik+1
            nk = len([ivi for ivi in range(ntra) if mtc[ik,ivi]==1])
            if nk==nn:
                kmax.append(ik)
        #print("kmax:",kmax)

        #If several tracks reach the maximum kmax, we choose the only one (maximum matches in time,
        #and if equal, minimum cumulated distance)
        ikm=-1
        if len(kmax)==1:
            ikm=kmax[0]
        elif len(kmax)>1:
            tnk=[]
            for ik in kmax: #parcours de tous les ik qui ont le même nombre de voisins
                tnk.append(np.mean([tn[ik,ii] for ii in range(ntra) if mtc[ik,ii]==1])) #number of common instants for each
            tmax=0
            for ivi in range(len(tnk)):
                if tnk[ivi]>tmax:
                    tmax=tnk[ivi]
            itmax=[]
            for ivi in range(len(tnk)):
                if tnk[ivi]==tmax:
                    itmax.append(kmax[ivi])
            if len(itmax)==1:
                ikm=itmax[0]
            else:
                distk=[]
                for it in itmax: #parcours de tous les it qui ont le même nombre de voisins et meme nombre d'instants
                    distk.append(np.mean([tdist[ik,ii] for ii in range(ntra) if mtc[ik,ii]==1])) #number of common instants for each
                dmin=Tools.maxval
                for ivi in range(len(distk)):
                    if distk[ivi]<dmin:
                        dmin=distk[ivi]
                for ivi in range(len(distk)): #We assume there is only one ...
                    if distk[ivi]==dmin:
                        ikm=itmax[ivi]
       
        #Finalization : all must have the same name ...
        if nn>=minmb and ikm>-1:
            print("Finally we keep kmax=",ikm,", which has ", nn, " close member" )
            lind=[]
            mtc[ikm,ikm]=1
            numane=numane+1
            nam0=prefname+"-"+str(numane)
            lname.append(nam0)
            ltraj2[ikm].name=nam0
            for ij in range(len(ltraj2)):
                #print(ikm,ij,mtc[ikm,ij])
                if mtc[ikm,ij]==1:
                    ltraj2[ij].name=nam0
                    lind.append(ij)
            #print("kmax ... we should keep ",lind," as ",nam0)
        #Update of ltraj2 for the next loop
            lind.reverse()
            #print("lind=",lind)
            ltraj3=copy.deepcopy(ltraj2)
            for ivi in lind:
                ltrajout.append(ltraj2[ivi])
                #print(len(ltraj3),ivi)
                ltraj3.pop(ivi)
            del ltraj2
            ltraj2 = ltraj3

        else:
            gook = False

        del mtc
        del tdist
        del tn
        del totdist

    return ltrajout, lname


def ConvertIBTRACS(filename, timetraj , domtraj, diags=[]):
    #Converts a IBTRaCS input file to a list of trajectories
    #Has been tested for last3years csv file : https://www.ncei.noaa.gov/data/international-best-track-archive-for-climate-stewardship-ibtracs/v04r00/access/csv/ibtracs.last3years.list.v04r00.csv
    #Inputs : IBTRaCS filename
    #timetraj : valid times of the output tracks (only 'start' and 'final' are used)
    #domtraj : domain of the output tracks (lonmin must be < lonmax)
    #diags: list of diagnostics
    #    the possible ones so far are "mslp_min", "ff10m_max"
    #The whole track of a system is kept if it is at least in the domain for 1 instant
    
    import pandas as pd
    
    ibtimefmt="%Y-%m-%d %H:%M:%S"
    
    #Instant boundaries from timetraj
    mintime = datetime.strptime(timetraj["start"],Inputs.time_fmt)
    maxtime = datetime.strptime(timetraj["final"],Inputs.time_fmt)

    #Diagnostics
    namediag={"ff10m_max":"USA_WIND","mslp_min":"USA_PRES"} #Correspondance des noms
    ldiag=[dg for dg in diags if dg in namediag.keys()]
    print("List of diagnostics: ", ldiag)
    convfct={"USA_WIND":1.0,"USA_PRES":100.0}
    
    ltraj=[]    
    data = pd.read_csv(filename)
    all_sid=list(data.SID.unique()) #all the systemID in the file
    all_sid=[all_sid[ivi] for ivi in range(len(all_sid)) if not (all_sid[ivi] == " " or all_sid[ivi] == "")]
    #print(all_sid)
    
    #Inputdef of the track
    indef = Inputs.inputdef()
    indef.update_input({"origin":"obs","database":"IBTRaCS","domain":domtraj})
    namobj="TC"
    
    for sid in all_sid: #Spans all SID
        #Step 1: check if the system is inside the domain and valid in time
        chktime=False
        chkdom=False
        datid = data.loc[data["SID"]==sid]
        for ivi in datid.index:

            if not (chktime and chkdom):
                inst = datetime.strptime(datid.loc[ivi]["ISO_TIME"],ibtimefmt)
                chktime = (inst >= mintime) and (inst <= maxtime)
                lon = float(datid.loc[ivi]["LON"])
                lat = float(datid.loc[ivi]["LAT"])
                chkdom = (lat >= domtraj["latmin"]) and (lat <= domtraj["latmax"])
                if domtraj["lonmax"]<180.0 and lon>180.0:
                    lon = lon-360.0
                chkdom = chkdom and (lon >= domtraj["lonmin"]) and (lon <= domtraj["lonmax"])

        #Step 2: Create SID track if in the domain and in time
        if chkdom and chktime:

            #Gets list of instants sorted in time
            linst=[]
            for ivi in datid.index:
                linst.append(datetime.strptime(datid.loc[ivi]["ISO_TIME"],ibtimefmt))
                linst.sort()

            #Create traj
            traj = DefTrack("Centre0D",basetime=linst[0],inputdef=indef)
            traj.name = sid
            traj.nameobj=namobj
                        
            #Add objects to traj
            for inst in linst:
                data_sub1 = datid.loc[datid['ISO_TIME'] == datetime.strftime(inst,ibtimefmt)]
                ind=list(data_sub1.index)[0]
                #print(data_sub1.loc[ind])
                lon = float(data_sub1.loc[ind]["LON"])
                if domtraj["lonmax"]<180.0 and lon>180.0:
                    lon = lon-360.0
                lat = float(data_sub1.loc[ind]["LAT"])
                
                objectm = DefObject("Centre0D", ldiag, [],lonc=lon,latc=lat,time=datetime.strftime(inst,Inputs.time_fmt))
                for dg in ldiag: #adding diagnostics
                    if len(data_sub1.loc[ind][namediag[dg]])>0 and not data_sub1.loc[ind][namediag[dg]]==' ':
                        setattr(objectm,dg,[missval, missval, float(data_sub1.loc[ind][namediag[dg]])*convfct[namediag[dg]]])
                        if dg=="mslp_min":
                            objectm.mslp_min[0]=objectm.lonc
                            objectm.mslp_min[1]=objectm.latc
                    else:
                        setattr(objectm,dg,[missval, missval, missval])

                if traj.name==sid and (not (data_sub1.loc[ind]["NAME"]=="NOT_NAMED")):
                    traj.name=data_sub1.loc[ind]["NAME"]
                
                traj.add_obj(objectm)

            traj.nobj=len(linst)                
            ltraj.append(traj)    

    return ltraj
