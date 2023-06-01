#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Traject/Interfaces for input data

https://github.com/UMR-CNRM/Traject
Traject is disributed under Cecill-C license
(For rights and obligations, see LICENSE.txt)

UMR-CNRM, Meteo-France & CNRS
"""
# =============================================================================
#                       INTERFACE FOR INPUT MET & CLIM DATA 
#                (somehow inspired by VORTEX, but still independent)
# =============================================================================

#------------------------------------------------------------------------------
                                # Importations
#------------------------------------------------------------------------------
import os
import epygram, json
import time
import datetime as dt
import numpy as np
import Tools

time_fmt = "%Y%m%d%H"
grib_default = "grib2mf" #Default value when "grib" is called

#os.environ["ECCODES_SAMPLES_PATH"] = ("/home/common/epygram/ext/eccodes/"
#                                      "share/eccodes/samples/")

#______________________________________________________________________________

#                  DEFINITION AND MANAGEMENT OF GRIDS
#______________________________________________________________________________

#The domain for tracking is a regular lat/lon grid defined by the parameters
#domtrack, restrack, that are provided by the user

#We use the epygram geometries and their computation methods

#______________________________________________________________________________

#                  DEFINITION AND MANAGEMENT OF INPUT DATA
#______________________________________________________________________________

class inputdef():
    #Generic parameters that define input data (for models)

    def __init__(self,namelist=""): #Constructor

        #runtype
        self.origin=None #can be "an" (analysis), "fc" (forecast) or "cl" (climate simulation) or "obs" (observation)

        if self.origin=="obs":
            self.database=None #Input database (IbTracks, manual, etc)
            self.domain=None #domain (of the input data)

        else:
            if self.origin=="an" or "fc":
                self.cutoff=None #can be "assim", "production" (for NWP or for reanalysis/reforecast) (optional)
            self.model=None #model (arpege, arome, ...)
            self.domain=None #domain (of the input data)
            self.experiment=None #Experiment name that is introduced by the user

        self.directory=None #Input directory
        self.filename=None #filename (that depends on the parameters)
        self.special_keys = {}

        if not namelist=="": #Builds from the namelist if it is provided
            self.read_input(namelist)

    def update_input(self,dico):
        #Updates the self objects with the keys and values in the dictionary dico

        if "member" in dico.keys():
            self.member = None
        if "scenario" in dico.keys():
            self.scenario = None
        if "experiment" in dico.keys():
            self.experiment = None
        if "param_nc" in dico.keys():
            self.param_nc = None
        if "param_file" in dico.keys():
            self.param_file = None
        if "filtered" in dico.keys():
            self.filtered = None

        self.__dict__.update(dico)

    def read_input(self,namelist):
        #Reads the namelist that defines the input data
        #(including the optional parameters)
        
        with open(namelist,'r') as json_file:
            datain = json.load(json_file)
        
        self.update_input(datain)

    def write_input(self,namelist):
        #Writes self in the namelist file
        #(including the optional parameters)
        
        with open(namelist,'w') as json_file:
            json.dump(self.__dict__,json_file,indent=4)
        
    def get_filename(self,basictime,term):
        #returns the filename that corresponds to:
        #basictime as validtime for climate, reanalysis or analysis
        #basictime and forecast term for forecasts
        #basictime shoud be given in string format "YYYYMMDDHH" (and written like this in the file)
        #term should be given in hours (integer value) (the number of digits is given in [term:3]

        fname = self.directory + self.filename

        dicts=self.__dict__

        for it in dicts.items(): #Replaces each [keys] by their values
            if isinstance(it[1],str):
                fname = fname.replace('['+str(it[0])+']',it[1])
        
        fname = fname.replace('[YYYY]',basictime[0:4])
        fname = fname.replace('[YYYYMM]',basictime[0:6])
        fname = fname.replace('[MM]',basictime[4:6])
        fname = fname.replace('[YYYYMMDD]',basictime[0:8])
        fname = fname.replace('[DD]',basictime[6:8])
        fname = fname.replace('[YYYYMMDDHH]',basictime[0:10])
        fname = fname.replace('[HH]',basictime[8:10])

        if fname.find('[term')>-0.5:
            ndigits=int(fname[fname.find('[term')+6])
            fname = fname.replace('[term:'+str(ndigits)+']',str(term).rjust(ndigits,'0'))

        return fname

    def get_ecdico(self,param):
        #eccode/epygram dictionary to extract the relevant field in grib format
        #different grib codes are implemented here
        #The grib format should be declared in inputdef.nativefmt

        gribfmt = self.nativefmt.lower()
        if gribfmt == "grib":
            gribfmt = grib_default

        #Selection of meteorological parameters
        if gribfmt=="grib2mf": #Cas du grib2 a Meteo-France

            #Some conventions of parameter names:
            #lower case letters
            #name in English
            #if a vertical level is specified, then the name is
                #parameter+"level" (if hPa)
                #parameter+"level"+m (if meter)
            if str(param) == 'mslp': #Mean sea level pressure
                dict = {'parameterCategory': 3, 'parameterNumber': 1,
                        'typeOfFirstFixedSurface': 101, 'typeOfLevel': 0}
            elif str(param) == 't2m': #temperature at 2m
                dict = {'parameterCategory': 0, 'parameterNumber': 0, 'name':
                    'Temperature', 'level': 2}
            elif str(param) == 'av850': #Absolute vorticity at 850hPa
                dict = {'parameterCategory': 2, 'parameterNumber': 10,
                        'typeOfFirstFixedSurface':100, 'level': 850}
            elif str(param) == 'rv850': #Relative vorticity at 850hPa
                dict = {'parameterCategory': 2, 'parameterNumber': 12,
                        'typeOfFirstFixedSurface':100, 'level': 850}
            elif str(param) == 'u500': #1st Horizontal component of wind at 500hPa
                dict = {'parameterCategory': 2, 'parameterNumber': 2,
                        'typeOfFirstFixedSurface':100, 'level': 500}
            elif str(param) == 'u700':
                dict = {'parameterCategory': 2, 'parameterNumber': 2,
                        'typeOfFirstFixedSurface':100, 'level': 700}
            elif str(param) == 'u850':
                dict = {'parameterCategory': 2, 'parameterNumber': 2,
                        'typeOfFirstFixedSurface':100, 'level': 850}
            elif str(param) == 'v500': #2nd Horizontal component of wind at 550hPa
                dict = {'parameterCategory': 2, 'parameterNumber': 3,
                        'typeOfFirstFixedSurface':100, 'level': 500}
            elif str(param) == 'v700':
                dict = {'parameterCategory': 2, 'parameterNumber': 3,
                        'typeOfFirstFixedSurface':100, 'level': 700}
            elif str(param) == 'v850':
                dict = {'parameterCategory': 2, 'parameterNumber': 3,
                        'typeOfFirstFixedSurface':100, 'level': 850}
            elif str(param) == 'u10m':
                dict = {'parameterCategory':2,'parameterNumber':2,'typeOfFirstFixedSurface': 103,
                        'typeOfSecondFixedSurface':255,'level':10}
            elif str(param) == 'v10m':
                dict = {'parameterCategory':2,'parameterNumber':3,'typeOfFirstFixedSurface': 103,
                        'typeOfSecondFixedSurface':255,'level':10}
            elif str(param) == 'btir':
                dict = {'parameterCategory': 5, 'parameterNumber': 7,
                        'scaleFactorOfCentralWaveNumber': 2,'scaledValueOfCentralWaveNumber': 9259259}
            elif str(param) == 'rr':
                dict = {'discipline':0, 'parameterCategory' : 1, 'parameterNumber' : 65}
                #dict = {'shortName':'rprate'}
            elif str(param) == 'rr1h':
                dict = {'discipline':0, 'parameterCategory' : 1, 'parameterNumber' : 65}
                #dict = {'shortName':'rprate'}
            elif str(param) == 'rr12h':
                dict = {'discipline':0, 'parameterCategory' : 1, 'parameterNumber' : 65}
            else:
                print("Parameter " + param + " not implemented in Inputs.py/inputdef/get_ecdico under grib2mf format")

        elif gribfmt=="grib1mf": #Cas du grib (version 1) a Meteo-France
            if str(param) == 'btir':
                dict = {"indicatorOfParameter":1,'indicatorOfTypeOfLevel':100,'level': 108}
            elif str(param) == 'rr':
                dict = {"indicatorOfParameter":150}
            elif str(param) == 'rr1h':
                dict = {"indicatorOfParameter":150}
            elif str(param) == 'rr12h':
                dict = {"indicatorOfParameter":150}
            else:
                print("Parameter " + param + " not implemented in Inputs.py/inputdef/get_ecdico under grib1mf format")
        else:
            print("Format " + self.nativefmt + " not implemented in Inputs.py/inputdef/get_ecdico")

        return dict

    def get_filinst(self,instraj,basetime=0):
    #Get list of files and instants that a single tracking should process
    #(depends on whether tracking is applied on forecasts --> instraj are terms and basetime is given,
    #or analyses --> instraj=timetraj are instants and basetime is not given (0))

        #Initialisations
        linst=[] #List of instants
        lfile=[] #List of files to process

        if self.origin=="fc": #Forecast data is processed: algo is applied on term steps
            #print("forecast data...")
            termtraj = instraj

            init=0
            if 'init' in termtraj:
                init=termtraj['init']

            for t in range(init,termtraj['final']+termtraj['step'],termtraj['step']):
                linst.append(dt.datetime.strptime(basetime,time_fmt)+dt.timedelta(hours=t))
                lfile.append(self.get_filename(basetime,t))

        elif self.origin=="an" or self.origin=="cl": #algo is applied on instants
            #print("climate or analyses...")
            timetraj=instraj

            t=0
            while dt.datetime.strptime(timetraj["start"],time_fmt)+dt.timedelta(hours=t) <= dt.datetime.strptime(timetraj['final'],time_fmt):
                ctime = dt.datetime.strptime(timetraj["start"],time_fmt)+dt.timedelta(hours=t)
                linst.append(ctime)
                lfile.append(self.get_filename(dt.datetime.strftime(ctime,time_fmt),-1))
                t = t + int(timetraj['step'])
        else:
            print("ABORT - indf.origin not defined - Should be fc, an or cl")
            exit()

        return lfile, linst

    def get_vortex(self,instraj,basetime=0):
        #Get the file using vortex, corresponding to the description given by:
        #self resource, list of instants instraj
        #the directory in self is not used and is changed to match the vortex environment

        import vortex, olive 

        #Generate the dictionary of the vortex resource
        dicin={'namespace':'vortex.multi.fr', #Default
                    'block':'forecast',
                    'vapp':self.model,
                    'vconf':'4dvarfr',
                    'cutoff':'assim',
                    'origin':'historic',
                    'kind':"gridpoint",
                    'nativefmt':self.nativefmt,
                    'model':self.model,
                    'suite':self.experiment.lower(),
                    'geometry':self.domain,
                    'date':'',
                    'term':'00',
                    'shouldfly':True,
                    'fatal':True
                      }

        if self.model=="arome":
            dicin['vconf'] = '3dvarfr'
        if self.origin=="fc":
            dicin['cutoff'] = 'production'
        if not (self.experiment.lower()=="oper" or self.experiment.lower()=="dble"): 
            dicin['suite']='olive'
            dicin['experiment']=self.experiment.lower()

        if self.nativefmt in ["grib1mf","grib2mf"]:
            dicin["nativefmt"] = "grib"

        dicts=self.__dict__

        if 'member' in dicts.keys():
            dicin['member'] = dicts['member']
            if self.model=="arpege":
                dicin['vconf'] = 'pearp'
            elif self.model=="arome":
                dicin['vconf'] = 'pefrance'

        if self.origin=="an":
            for inst in comptimes(instraj):
                dicin['date']=inst+'00'
                rhs=vortex.toolbox.input(dicin)
                rhs[0].get()
                self.directory=os.environ["MTOOLDIR"]+'cache/vortex/'+dicin['model']+'/'+dicin['vconf']+ \
                              '/'+self.experiment.upper()+'/[YYYYMMDD]T[HH]00A/forecast/'
                self.filename='grid.'+dicin['model'].lower()+'-forecast.'+dicin['geometry'].lower()+'+0000:00.grib'
        else:
            dicin['date']=basetime+'00'
            for inst in compterms(instraj):
                dicin['term']=str(inst).rjust(3,'0')

                rhs=vortex.toolbox.input(dicin)
                rhs[0].get()
                #print(os.environ["MTOOLDIR"]+'cache/vortex'+rhs[0].uridata['path'])
                if 'member' in dicin.keys():
                    self.directory=os.environ["MTOOLDIR"]+'cache/vortex/'+dicin['model']+'/'+dicin['vconf']+ \
                              '/'+self.experiment.upper()+'/[YYYYMMDD]T[HH]00P/mb[member]/forecast/'
                else:
                    self.directory=os.environ["MTOOLDIR"]+'cache/vortex/'+dicin['model']+'/'+dicin['vconf']+ \
                              '/'+self.experiment.upper()+'/[YYYYMMDD]T[HH]00P/forecast/'
                self.filename='grid.'+dicin['model'].lower()+'-forecast.'+dicin['geometry'].lower()+'+[term:4]:00.grib'

        return

def compterms(termtraj):
    #Outputs the list of fc terms corresponds to termtraj:{'start','final','step'}

    lterm=[]
    init=0
    if 'init' in termtraj:
        init=termtraj['init']

    for t in range(init,termtraj['final']+termtraj['step'],termtraj['step']):
        lterm.append(t)
    
    return lterm

def comptimes(timetraj):
    #Outputs the list of instants that corresponds to timetraj:{'start','final','step'}

    ltime=[]
    t=0
    tstep=max(1,int(timetraj['step']))#Minimum step=1h (sinon ca boucle dans le vide) !
    while dt.datetime.strptime(timetraj["start"],time_fmt)+dt.timedelta(hours=t) <= dt.datetime.strptime(timetraj['final'],time_fmt):
        ctime = dt.datetime.strptime(timetraj["start"],time_fmt)+dt.timedelta(hours=t)
        ltime.append(dt.datetime.strftime(ctime,format=time_fmt))
        t = t + tstep
        #print(ctime,t,dt.datetime.strptime(timetraj['final'],time_fmt))
    
    return ltime


class algodef():
    #Generic parameters that define the algorithm parameters

    def __init__(self,namelist=""): #Constructor

        #Name of the algorithm (corresponds to the .py filename)
        self.name=None

        #Object class
        self.classobj=""

        #Default object name (could be changed by the algorithm if relevant)
        self.nameobj=""

        #Algorithm parameters (depends on the algorithm)
        self.varalgo={}

        #Domain definition for tracking
        self.domtraj=None

        #Filter values (in km of the input fields)
        self.parfilt={}

        #Change of resolution of the input fields
        self.parres={}

        #Diagnostic parameters (asked by the user)
        self.diag_parameter={}

        #Parallel options (optional)
        self.parallel={}

        #Special field options (optional, depends on the algorithm)
        self.specfields=None

        if not namelist=="": #Builds from the namelist if it is provided
            self.read_input(namelist)

    def update_input(self,dico):
        #Updates the self objects with the keys and values in the dictionary dico

        self.__dict__.update(dico)

    def read_input(self,namelist):
        #Reads the namelist that defines the algorithm parameters
        
        with open(namelist,'r') as json_file:
            datain = json.load(json_file)
        
        self.__dict__.update(datain)


    #------------------------------------------------------------------------------
                            # Traitement des fichiers
    #------------------------------------------------------------------------------

def get_filetype(inputfile):
    '''
    Gets the type of inputfile
    The output is:
    ftype=a string in lower_case that can be grib, nc, json, hdf5, etc
    epytype=a string that corresponds to the epygram fileformat (if relevant), "" otherwise
    '''

    ftype=""
    epytype=""

    lname=inputfile.split('.')
    filefmt=lname[-1].lower()

    if filefmt in ['grib','grib2']:
        ftype="grib"
        epytype='GRIB'
    elif filefmt in ['nc','netcdf']:
        filetype='nc'
        epytype='netCDF'
    elif filefmt in ['json']:
        filetype='json'
    elif filefmt in ['hdf5','hdf','h5']:
        filetype='hdf5'

    return ftype, epytype

def check_file(filename,indf,param,param0):
    #computes the final name of the file (substitution of param_file name)
    #and checks is the file exists
    #and outputs the format named as for epygram

    #Substition of parameter in filename (if needed)
    if '[param_file]' in filename:
        param2=param
        #print(indf.__dict__)
        for par in indf.param_file.items():
            #print(par,param,param0)
            if par[0]==param:
                param2=par[1]
            elif par[0]==param0:
                param2=par[1]
        filename2 = filename.replace('[param_file]',param2)
    else:
        filename2=filename

    gook = os.path.exists(filename2)

    ftype, epyfmt = get_filetype(filename2)

    return gook, filename2, epyfmt

def get_paramnc(indf,param,lev):
    #computes the parameter name to be used in the netcdf file

    param2=param
    for par in indf.param_nc.items():
        if par[0]==param+lev:
            param2=par[1]
        elif par[0]==param:
            param2=par[1]
    #print("get_paramnc",param,param2,lev)

    return param2

def split_filename(fic):
    #Separates filename directory, name and extension
    
    lfic=fic.split('/')
    fic2=lfic[-1]
    rep=fic[0:len(fic)-len(fic2)]

    fics=fic2.split('.')
    ext=fics[-1]
    fname=fic2[0:len(fic2)-len(ext)-1]

    return rep, fname, ext

def split_param(param):
    #Separates parameter name and level from param
    #Can be used only for Traject parameters (defined in Inputs.get_ecdico())

    if param[0:2]=="rv":
        par="rv"
        lev=param[2:5]
    elif param[0:1]=="u":
        par="u"
        lev=param[1:4]
    elif param[0:1]=="v":
        par="v"
        lev=param[1:4]
    elif param[0:5]=="fgust":
        par="fgust"
        lev=param[5:7]
    elif param[0:2]=="ff":
        par="ff"
        lev=param[2:5]
    elif param[0:1]=="t":
        par="t"
        lev=param[1:4]
    elif param in ["mslp","btir"] or param[0:2]=="rr":
        par=param
        lev=""
    else:
        print("Error in Inputs.split_param : the parameter name "+param+" should be processed here")
        exit()

    return par, lev

def read_nc(f1,inst,parnc,levnc):
    #reads the ncfile netcdf file to get field of parameter parnc at level levnc at instant inst
    #The output should be a single epygram field

    ll=f1.listfields() #=f1.find_fields_in_resource(dict(level=850))

    single_level=True
    single_time=True

    if not levnc=="":
        #Find list of levels
        levvar=["level","levels","lev"] #possible names for levels in nc file...
        levname=[]
        for il in levvar:
            for jl in ll:
                if il==jl:
                    levname=il
        if not levname==[]:
            single_level=False
            levels=f1.readfield(levname).data.tolist()
            #print("--Levels--")
            #print(single_level)
            #print(levels,levnc)
        else:
            single_level=True
            #Infer the nc name in the file-->
            parnc2=parnc
            if parnc in ll:
                parnc2=parnc
            elif parnc+levnc in ll:
                parnc2=parnc+levnc
            parnc=parnc2
            print("Warning - Could not find vertical levels in nc file")
            print("Warning - Available fields are: ", ll)
            print("Warning - So I take the single field that corresponds to --> " + parnc)

    #Find list of times
    timevar=["time","times"] #possible names for times in nc file...
    timename=[]
    for il in timevar:
        for jl in ll:
            if il==jl:
                timename=il
    if not timename==[]:
        single_time=False
        #print(timename)
        times=f1.readfield(timename).validity.get()
        #print("--Times--")
        #print(type(times))
        #print(times)
        #print(inst)
        if isinstance(times,dt.date):
            single_time=True
    else:
        single_time=True
        print("Warning - Could not find times in nc file")
        print("Warning - Available fields are: ", ll)
        print("Warning - So I take the single field available...")

    if not single_time:
        indt=[x for x in range(len(times)) if times[x]==inst]
    else:
        indt=[-1]

    if not single_level:
        indl=[x for x in range(len(levels)) if float(levels[x])==float(levnc)]
    else:
        indl=[-1]

    #check if the latitudes are reversed or not...
    adhb={}
    latvar=["lat","latitude","latitudes"] #possible names for latitudes in nc file...
    latname=[]
    for il in latvar:
        for jl in ll:
            if il==jl:
                latname=il
    if not latname==[]:
        try :
            lat1=f1.readfield(latname,only={latname:1}).data
        except:
            print("Warning ... could not read latitude dimension ... we assume Yorder is correct")
        else:
            #lat est bien un vecteur et on peut le lire
            lat2=f1.readfield(latname,only={latname:2}).data
            if lat2<lat1:
                adhb={'reverse_Yaxis':True}

    if len(indt)==1 and len(indl)==1:
        if not single_level and not single_time:
            print("Extraction of data at time " + str(times[indt[0]]) + " and level "+ str(levels[indl[0]]))
            fs = f1.readfield(parnc,only={timename:indt[0],levname:indl[0]},adhoc_behaviour=adhb)
        elif not single_time:
            print("Extraction of data at time " + str(times[indt[0]]))
            fs = f1.readfield(parnc,only={timename:indt[0]},adhoc_behaviour=adhb)
        elif not single_level:
            print("Extraction of data at level " + str(levels[indl[0]]))
            fs = f1.readfield(parnc,only={levname:indl[0]},adhoc_behaviour=adhb)
        else:
            print("Extraction of data")
            fs = f1.readfield(parnc)

    else:
        print("Several instants and levels found in the file...")
        print("abort")
        exit()

    #fs.what()

    return fs

def open_field(filename,inst,indf,param):
    #Open the epygram field given in filename, at instant inst,
    #described by indf and param
    #Input file can be grib or netcdf
    #Output : a single epygram field

    #Definition of input parameter
    par, lev = split_param(param)
    param0 = param
    par0 = par
    par1 = ""
    #print("split :",par, lev)

    #Processing of some special cases using indf special keys
    if par=="rv" and "rv_av" in indf.special_keys: #relative vorticity is computed from absolute vorticity
        par0="av"
        param0=param.replace("rv","av")
    if par=="ff" and "ff_uv" in indf.special_keys: #wind module (ff) is computed from u and v
        par0="u"
        param0=param.replace("ff","u")

    #check of input file and possible substitution of param_file in name
    gook, filename2, epyfmt = check_file(filename,indf,par0,param0)

    if gook and epyfmt=="GRIB":
        f1 = epygram.formats.resource(filename=filename2,openmode="r",fmt = epyfmt)
        f1.open()
        fs = f1.readfield(indf.get_ecdico(par0+lev))
        #We assume here that there is only one instant in the grib file ...
        #if it is not the case, some adaptation is required

    elif gook and epyfmt=="netCDF":

        f1 = epygram.formats.resource(filename=filename2,openmode="r")

        #Definition of nc keys from par0
        par1 = get_paramnc(indf,par0,lev)
        #print('param_nc',par0,par0+lev,par1)
        fs = read_nc(f1, inst, par1, lev)

    else:
        print("File "+filename2+" does not exist ... ABORT")
        exit()

    f1.close()

    return fs, par, lev

def extract_field(filename,inst,indf,lfields,filout):
    #extracts the fields in lfields from filename at instant inst and writes it in filout
    #Manages grib and nc format
    #The output file has the same format as the input file

    ftype, epytype = get_filetype(filename)

    rout=epygram.formats.resource(filout, "w", fmt=epytype)

    for var in lfields:
       fld, par, lev=open_field(filename,inst,indf,var)
       rout.writefield(fld)

    rout.close()

    return

def filter_field(filename,inst,indf,filt_var,domtraj,dres,repout,basetime,filesuf,subnproc):
    #filter the fields at scale (km) given in filt_var and extracts on domtraj
    #at resolution given by dres (for every parameter)
    #the list_var variables from filename at instant inst
    #Input file can be grib or nc
    #output file is in netcdf format (rep+FILT$res+filesuf)

    param_file={}
    param_nc={}
    fileres=[] #list of resolution files
    rout=[]
    ivi=0

    for var, filteff in filt_var.items():
        if dres[var] not in fileres:
            filout=repout+"/FILT"+str(dres[var]).replace('.','p')+"-"+filesuf
            rout.append(epygram.formats.resource(filout, "w", fmt="netCDF"))
            fileres.append(dres[var])
            ivi=ivi+1

        ic=fileres.index(dres[var])
        fld = extract_data(filename,inst,indf,var,domtraj,dres[var],basetime,subnproc,filteff)
        fld.fid['netCDF'] = var
        rout[ic].writefield(fld)
        param_file[var]=str(dres[var]).replace('.','p')
        param_nc[var]=var

    for ic in range(len(fileres)):
        rout[ic].close()

    return param_file, param_nc


def extract_data(filename,inst,indf,param,domtraj,res,basetime,subnproc,filtrad=0):
    '''
    Extracts data from filename that corresponds to:
    - instant inst,
    - parameter param,
    - on the grid given by domtraj (boundaries) and res (resolution),
    given the definition of input data as in indf.

    subnproc: number of cores used for pyresample
    If filtrad is declared and if it is >0, then smoothing is applied
    at the equivalent resolution of filtrad (in km).

    The output is a single epygram fieldset.
    '''

    #Opens the file and gets the fieldset
    fs, par, lev = open_field(filename,inst,indf,param)

    #Gets the original grid
    (lons, lats) = fs.geometry.get_lonlat_grid()

    # Extraction of domain characteristics
    lonmin=np.min(lons)
    lonmax=np.max(lons)
    latmin=np.min(lats)
    latmax=np.max(lats)

    #Test longitudes
    if lonmax>180.0:
        fs.global_shift_center(-180)
        (lons, lats) = fs.geometry.get_lonlat_grid()
        lonmin=np.min(lons)
        lonmax=np.max(lons)
    res0=lons[0][1]-lons[0][0]

    #Special case: relative vorticity is computed from absolute vorticity
    if par=="rv" and "rv_av" in indf.special_keys:
        fs = comp_rv(fs)

    #Special case: wind module ff is computed from u and v
    if par=="ff" and "ff_uv" in indf.special_keys:
        fs = comp_ff(par,lev,indf,inst,basetime,domtraj,res,filtrad,subnproc)

    #Special case: decumulate rain rate
    rrdecum, rrbefore = ifdecum(indf)
    if rrdecum and par[0:2]=="rr" and len(par)>2:
        fs = comp_rr(par,indf,inst,basetime,domtraj,res,subnproc,filtrad)

    #If there is no change of resolution nor smoothing --> extract subdomain
    if filtrad<1e-6 and abs(res0-res)<1e-6:
        fld=fs.resample_on_regularll({"lonmin":domtraj["lonmin"],"latmin":domtraj["latmin"],
        "lonmax":domtraj["lonmax"],"latmax":domtraj["latmax"]},res, nprocs=subnproc)
    else:
    #There is some change of resolution or smoothing --> resampling
        if filtrad>1e-6:
            filtrad2=filtrad
            ngb = int((res/res0)*(res/res0)*(filtrad/(res*100))*(filtrad/(res*100))*4)
        else:
            filtrad2 = res*100.0 #we apply an equivalent radius to the new grid
            ngb = int((res/res0)*(res/res0)*(filtrad2/(res*100))*(filtrad2/(res*100))*4)+1
        fld=fs.resample_on_regularll({"lonmin":domtraj["lonmin"],"latmin":domtraj["latmin"],
        "lonmax":domtraj["lonmax"],"latmax":domtraj["latmax"]},res,weighting='gauss',
                radius_of_influence=filtrad2*1e3,neighbours=ngb,sigma=(filtrad2/2.0)*1e3, reduce_data=True, nprocs=subnproc)

    return fld

def extract_domain(filename,inst,indf,param):

    '''Extracts the total domain from the initial grid
    Input : filename, instant, indf and a parameter that is present in the file
    Output : domtraj that covers the whole grid, and
            associated lon et lat arrays.
    Output longitudes are set to be in the range ]-180,180].

    Useful in case domtraj is not given by the user - Then tracks are computed
    on the total domain (can be costly!)'''

    #Open file and get field
    fs, par, lev = open_field(filename,inst,indf,param)

    # Extraction des lon, lat et données du domaine global
    (lons, lats) = fs.geometry.get_lonlat_grid()
    lonmin=np.min(lons)
    lonmax=np.max(lons)
    latmin=np.min(lats)
    latmax=np.max(lats)
    res=int((lons[0][1]-lons[0][0])*1e6)/1e6

    #Test longitudes
    if lonmax>180.0:
        fs.global_shift_center(-180)
        (lons, lats) = fs.geometry.get_lonlat_grid()
        lonmin=np.min(lons)
        lonmax=np.max(lons)

    domtraj={"lonmin":lonmin,"lonmax":lonmax,"latmin":latmin,"latmax":latmax}

    return domtraj, res, lons[0,:].tolist(), lats[:,0].tolist() #for lons and lats, we return list

def comp_rv(f1):
    #Compute relative vorticity from f1 absolute vorticity
    #from epygram f1 fieldset

    from Tools import omega

    (lons, lats) = f1.geometry.get_lonlat_grid()
    sinlat = 2*omega*np.sin(lats*np.pi/180.0)
    sinlat2=f1.deepcopy()
    sinlat2.data=sinlat
    f1.operation("-",sinlat2)

    return f1

def ifdecum(indf):
    #Tests indef to verify is decumulation should be done or not
    #Decumulation is done only if rr_before or rr_after is in special_keys.

    rrdecum=False
    rrbefore=False

    #Sets if accumulation is done before inst (rr_before=True) or after
    if "special_keys" in indf.__dict__:
        if "rr_before" in indf.special_keys:
            rrdecum=True
            rrbefore=True
        elif "rr_after" in indf.special_keys:
            rrdecum=True
            rrbefore=False

    return rrdecum, rrbefore

def comp_rr(parname,indf,inst,basetime,domtraj,res,subnproc,filtrad):
    #Decumulates rainfall according to parname and indf data at instant inst
    #if rr_after in indf.special_keys then the difference is done after inst,
    #else it is done before inst (default option)
    #Decumulation is done only if rr_before or rr_after is in special_keys.

    btime=dt.datetime.strptime(basetime,time_fmt)

    rrdecum, rrbefore = ifdecum(indf)

    if rrdecum:

        #Gets the accumulation frequency (parname should be rr1h or rr6h or ...)
        freq=int(parname.replace('rr','').replace('h',''))

        #print("comp_rr")
        #print("basetime,inst,freq : ",basetime,inst,freq)
        if indf.origin=='fc':
            #Defines the first and second instants & read the fields
            if rrbefore:
                inst1=inst-dt.timedelta(hours=freq)
                inst2=inst
            else:
                inst1=inst
                inst2=inst+dt.timedelta(hours=freq)
            t1 = int((inst1-btime).days*24 + ((inst1-btime).seconds)/3600)
            t2 = int((inst2-btime).days*24 + ((inst2-btime).seconds)/3600)
            #print(basetime)

            if t1>0:
                fname1=indf.get_filename(basetime,t1)
                fname2=indf.get_filename(basetime,t2)   
                #print(fname1, inst1, t1)
                #print(fname2, inst2, t2)
                f1=extract_data(fname1,inst1,indf,'rr',domtraj,res,subnproc,filtrad)
                f2=extract_data(fname2,inst2,indf,'rr',domtraj,res,subnproc,filtrad)
                f2.operation("-",f1)
            else:
                fname2=indf.get_filename(basetime,t2)   
                f2=extract_data(fname2,inst2,indf,'rr',domtraj,res,subnproc,filtrad)
        else:
            print("Fatal error in Inputs.comp_rr(): rain decumulation works only on forecasts")
            print("You may remove "+ indf.special_keys +" key to process without decumulation")
            exit()
        #print("Minimum/maximum rainfall values : ",f2.min(),"/",f2.max())

    else: #no decumulation: we take the value at instant inst in the file
        print("Fatal error in Inputs.comp_rr(): please add rr_before or rr_after in inputdef/special_keys")
        print("to process rain decumulation")
        exit()
        #t2 = int(((inst-btime).seconds)/3600)
        #fname2=indf.get_filename(basetime,t2)   
        #f2=extract_data(fname2,inst,indf,'rr',domtraj,res,subnproc,filtrad)

    return f2

def comp_ff(parname,lev,indf,inst,basetime,domtraj,res,filtrad,subnproc):
    #Computes wind module (ff) from u and v
    #if ff_uv in indf.special_keys

    term = Tools.comp_difftime(basetime,inst)
    #print("Compute ff ", basetime, '+', term, "Valid time=",inst)
    #print("Filter= ", filtrad, " Resolution: ",res)

    fname1=indf.get_filename(basetime,term)
    uf=extract_data(fname1,inst,indf,'u'+lev,domtraj,res,subnproc,filtrad)
    vf=extract_data(fname1,inst,indf,'v'+lev,domtraj,res,subnproc,filtrad)

    wind = epygram.fields.make_vector_field(uf,vf)
    ff = wind.to_module()

    return ff

