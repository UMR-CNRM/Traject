# -*- coding: utf-8 -*-
"""
Some Visualisation functions associated to Traject

https://github.com/UMR-CNRM/Traject
Traject is disributed under Cecill-C license
(For rights and obligations, see LICENSE.txt)

UMR-CNRM, Meteo-France & CNRS


"""

import os
from traject import *
import pandas as pd
import matplotlib as mpl
import matplotlib.dates as mdates
from datetime import datetime
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import cartopy.feature as cfeature
import epygram
import copy

#---------------------------------------------------------
# Some simple plots to plot tracks and diagnostics
#---------------------------------------------------------

def plot_diag(trakin,diagname,ficout,**kwargs):
    #Plot diagnostic "diagname" in file ficout
    #Options in kwargs:
    #leg : legend of each track (list)
    #cum (boolean) : plot cumulation instead of actual value

    # create figure
    fig, ax = plt.subplots()

    if isinstance(trakin,str):
        ltraj=Read(trakin,**kwargs)
    else:
        ltraj=trakin
    lname=[]
    
    #conversion factor
    fct=1.0
    if diagname=="mslp":
        fct=0.01 #(conversion to hPa)

    for ivi in range(len(ltraj)):
        plot_single_diag(ltraj[ivi],diagname,factor=fct,**kwargs)
        #print(ltraj[ivi].name)
        lname.append(ltraj[ivi].name)        
    if "leg" in kwargs:
        lname=kwargs["leg"]

    # set graphical parameters
    myFmt = mdates.DateFormatter('%d/%H')
    ax.xaxis.set_major_formatter(myFmt)

    plt.xlabel('Time')
    plt.ylabel(diagname)
    plt.legend(lname,ncol=3)

    plt.savefig(ficout)

def plot_track(trakin,ficout,**kwargs):
    #Plot tracks on a map in file ficout
    #Options in kwargs:
    #dom : other domain than domtraj
    #leg : legend of each track (list)

    if isinstance(trakin,str):
        ltraj=Read(trakin,**kwargs)
    else:
        ltraj=trakin

    # create figure
    ax = plt.axes(projection=ccrs.PlateCarree())

    lname=[]

    for ivi in range(len(ltraj)):
        #print(ltraj[ivi].name)
        #print(ltraj[ivi].__dict__)
        plot_single_track(ltraj[ivi])
        lname.append(ltraj[ivi].name)        
    if "leg" in kwargs:
        lname=kwargs["leg"]


    # set graphical parameters
    ax.coastlines(resolution = '50m')
    gl = ax.gridlines(draw_labels=True, linestyle = '--', linewidth = 0.5)
    if "dom" in kwargs:
        domtraj=kwargs["dom"]
    else:
        domtraj=ltraj[0].algodef["domtraj"]
    ax.set_extent([domtraj["lonmin"],domtraj["lonmax"],domtraj["latmin"],domtraj["latmax"]])
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    plt.legend(lname,title='Tracks',ncol=3)

    plt.savefig(ficout)


#-----------------------------------------------------------#
#####   Routines to process a single track (used by other routines)
#-----------------------------------------------------------#
def plot_single_diag(tra,diagname,factor=1.0,col='',**kwargs): #Plot diagnostic "diagname"

    # create empty lists
    time = []
    diag = []
    iscum = False #Should we plot cumulation ?
    if "cum" in kwargs:
        if kwargs["cum"]:
            iscum = True
    cum = []

    #print(tra.inputdef["member"],tra.name)

    # read variables
    for ivi in range(tra.nobj):
        #print(tra.traj[ivi].__dict__['diags'], diagname)
        if diagname in tra.traj[ivi].__dict__['diags']:
            time.append(tra.traj[ivi].__dict__['time'])
            diag.append(tra.traj[ivi].__dict__[diagname][2])
            if iscum:
                if len(cum)==0:
                    cum.append(tra.traj[ivi].__dict__[diagname][2])
                else:
                    cum.append(cum[-1]+tra.traj[ivi].__dict__[diagname][2])

    # convert for plot
    timestamp = pd.to_datetime(time, format='%Y%m%d%H')
    dconv = np.array(diag) * factor
    
    if iscum:
        dcum = np.array(cum) * factor
    

    #print(tra.nobj)
    #print(timestamp)
    #print(time)

    # plot time evolution
    if len(timestamp)>1:
        if not col=='':
            if not iscum:
                plt.plot(timestamp,dconv,color=col)
            else:
                plt.plot(timestamp,dcum,color=col,linestyle='--')
        else:
            if not iscum:
                plt.plot(timestamp,dconv)
            else:
                plt.plot(timestamp,dcum,linestyle='--')
    else:
        plt.plot(timestamp,dconv,'o')

#-----------------------------------------------------------#
def plot_single_track(tra,col=''): #Plot track
    # create empty lists
    latc = []
    lonc = []
    # read variables
    for ivi in range(tra.nobj):
        latc.append(tra.traj[ivi].__dict__['latc'])
        lonc.append(tra.traj[ivi].__dict__['lonc'])
    # plot map
    if len(lonc)>1:
        if not col=='':
            plt.plot(lonc,latc,color=col)
        else:
            plt.plot(lonc,latc)
    else: # track contains one point only
        plt.plot(lonc,latc,'o')


#---------------------------------------------------------
# Plots for probabilistics diagnostics
#---------------------------------------------------------

def map_plot(ltraj, diag, typeplot="line", leg=[], opt = "std", fig={}, filename="", colormap="viridis", colorlevels=[], diagthr=0.0, diagrad=0.0, centre="obj",dom=[], **kwargs):
    #ltraj, #list of tracks to take into account in the time plot
    #diag, #diagnostic name to plot (default is nam=””: the tracked object is plotted)
    #typeplot:, #can be “line” (default),  “dot” or “strike”
    #fig = figure name and axis on which to plot (optional input, default is {})
    #leg=[], #legend for every track (if default, use member number) if relevant (“line” or “dot”),
    #colorlevels:[] # thresholds for different colors (default is []) if typeplot=“dot”, or scale of probabilistic values if typeplot="“strike”
    #colormap=, #colormap (used for leg, default is “viridis”)
    #diagthr=0.0, #threshold on diag to take into account the object in the plot,
    #diagrad=0.0, #radius (km) to add around the diagnostic
    #centre="obj" #if we plot the object centre or the diagnostic centre (if kind=max or min)
    #“opt”: # list of options, can be “max” (for maximum value along the track), "std" for whole domain (default), "zoom": zoom on the track features, "user": for user-defined domain (dom)
    #filename, #optional, output filename
    #kwargs #arguments passed to plot() routine

    #Plot definition
    if fig=={}:
        #Create new figure and 1 axis
        fig, ax0 = start_figure_map(**kwargs)
        #Set up domain limits
        lonmin,lonmax,latmin,latmax,res = set_dom_limits(ltraj,opt,diag=diag,diagrad=diagrad,dom=dom)
        ax0.set_extent([lonmin,lonmax,latmin,latmax])

    lax = fig.get_axes()
    cax = lax[0]

    #Plotting
    if typeplot=="line":
        plot_map_line(ltraj, diag, cax, leg, opt, colormap, **kwargs)
    elif typeplot=="dot":
        plot_map_dot(ltraj, diag, fig, cax, leg, opt, colormap, colorlevels, **kwargs)
    elif typeplot=="strike":
        plot_map_strike(ltraj, diag, fig, diagrad, diagthr, centre, cax, leg, opt, colormap, colorlevels, dom, **kwargs)
    else:
        print("Wrong typeplot option in time_plot() : typeplot=",typeplot)
        exit()

    #save and close figure
    if len(filename) > 0:
        fig.savefig(filename)

    return fig
    
def plot_map_line(ltraj, nam, cax, leg, opt, cmap, **kwargs):        
    #Plot individual members in ltraj as a function of time on cax
    #(called by map_plot)
    
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    
    #read list of colors
    ntraj=len(ltraj)
    lcol=get_col(cmap,ntraj)
    
    lname = []
    for ivi in range(ntraj):
        latc=[]
        lonc=[]
        tra=ltraj[ivi]
        for ivj in range(tra.nobj):
            #print(tra.traj[ivi].__dict__['diags'], diag)
            if nam=="":
                latc.append(tra.traj[ivj].__dict__['latc'])
                lonc.append(tra.traj[ivj].__dict__['lonc'])
            elif nam in tra.traj[ivj].__dict__['diags']:
                dia=Tools.guess_diag(nam,True)
                if not dia.kind[0]=='c':
                    lonc.append(tra.traj[ivj].__dict__[nam][0])
                    latc.append(tra.traj[ivj].__dict__[nam][1])
        
        if len(lonc)>1:
            cax.plot(lonc,latc,color=lcol[ivi])
        elif len(lonc)==1:
            cax.plot(lonc,latc,color=lcol[ivi],marker='x')
        else:
            print("WARNING in plot_map_line: nothing to plot...")
        if len(lonc)>0:
            if "member" in tra.inputdef:
                lname.append(tra.inputdef["member"][1:3])
            else:
                lname.append(tra.name)
        
    #legend
    if len(leg) > 0:
        lname=leg

    cax.legend(lname,title=nam+' tracks',ncol=3,loc="best")
    #Traits legende noires a cause de epygram_departments ... a corriger ... voir comment faire un 2e axe?
    
    return

def plot_map_dot(ltraj, nam, fig, cax, leg, opt, cmap, lev, **kwargs):  
    #Plot individual members in ltraj as dots on cax
    #if opt=="max", only the maximum value along the forecast is computed
    #(called by map_plot)
    
    #read list of colors
    if lev==[]: #a color different by member of ltraj
        ntraj=len(ltraj)
        lcol=get_col(cmap,ntraj)
        cbar=False
    elif len(lev)>1: # color depends on the diag level
        nlev=len(lev)
        gmap = mpl.cm.get_cmap(name=cmap,lut=nlev)
        lcol = [mpl.colors.to_hex(gmap(ivi), keep_alpha=True) for ivi in range(nlev)]
        #print(lcol)
        cbar = True
    else:
        print("Please include at least 2 colorlevels - ABORT")
        exit()
        
    #print(lcol)
        
    lname = []
    for ivi in range(len(ltraj)):
        lonc=[]
        latc=[]
        pdiag=[]
        if opt=="max":
            pdiag.append(-999.0)
            lonc.append(-999.0)
            latc.append(-999.0)
        tra=ltraj[ivi]
        for ivj in range(tra.nobj):
            if nam=="":
                lon0=tra.traj[ivj].__dict__['lonc']
                lat0=tra.traj[ivj].__dict__['latc']
                val=-999.0
            elif nam in tra.traj[ivj].__dict__['diags']:
                dia=Tools.guess_diag(nam,True)
                if not dia.kind[0]=='c':
                    lon0=tra.traj[ivj].__dict__[nam][0]
                    lat0=tra.traj[ivj].__dict__[nam][1]
                    val=tra.traj[ivj].__dict__[nam][2]
            if opt=="max":
                if val > pdiag[0]:
                    pdiag[0] = val
                    lonc[0] = lon0
                    latc[0] = lat0
            else:
                pdiag.append(val)
                lonc.append(lon0)
                latc.append(lat0)

        if lev==[]:
            #print(lonc,latc)
            cax.plot(lonc,latc,'o', color=lcol[ivi]) #, **kwargs)
        else:
            #print(pdiag)
            #print(lev)
            for ivc in range(len(pdiag)):
                col=lcol[np.digitize(pdiag[ivc],lev)-1]
                cc=cax.plot(lonc[ivc],latc[ivc],'o', color=col)
        #print(ltraj[ivi].inputdef["member"])
        if "member" in ltraj[ivi].inputdef:
            lname.append(ltraj[ivi].inputdef["member"][1:3])
        else:
            lname.append("")
        
    #legend
    if len(leg)>0:
        lname=leg

    if nam=="":
        dia = Tools.guess_diag(ltraj[0].algodef["specfields"]["track"],True)
    else:
        dia = Tools.guess_diag(nam,True)
    if cbar:
        add_colorbar(fig,gmap,lev,dia)
        #cba = fig.add_axes([0.93, 0.2, 0.02, 0.6])
        #cb1 = mpl.colorbar.ColorbarBase(cba,cmap=gmap,orientation='vertical')
        #cb1.set_ticks([float(ivi)/(nlev-1) for ivi in range(nlev)])
        #cb1.set_ticklabels(lev)
        #cb1.set_label(dia.nicename+" ("+dia.unit+") : dots")
    else:
        plt.legend(lname,ncol=3)
        cax.set_ylabel(dia.nicename+" ("+dia.unit+") : dots")    
    return

def plot_map_strike(ltraj, diag, fig, diagrad, diagthr, centre, cax, leg, opt, colormap, colorlevels, dom, **kwargs):
    #Plots strike probabilities 
    #(called by map_plot())

    #computation of strike (use of object masks)
    lonmin,lonmax,latmin,latmax,res = set_dom_limits(ltraj,opt,diag=diag,diagrad=diagrad,dom=dom)
    images, N_obj, lons, lats = init_strike(ltraj,[lonmin,lonmax,latmin,latmax],res)
    ivi=0
    for traj in ltraj:
        for obj in traj.traj:
            images[ivi] = obj.mask(lons,lats,diag,diagthr,diagrad,centre)
            ivi = ivi +1
    images[images > 0] = 1
    
    summed_img = np.sum(images, axis=0)
    imax = np.max(summed_img) #Normalisation 1

    if colorlevels==[]:
        maxc=1.0
    else:
        maxc=max(colorlevels)

    nlev=10
    col=[ivi/nlev for ivi in range(nlev+1)]

    if imax>0:
        summed_img = summed_img / imax * maxc
    
    gmap = mpl.cm.get_cmap(name=colormap,lut=nlev)
    plt.contourf(lons,lats,summed_img,col,cmap=gmap,transform=ccrs.PlateCarree())
    cba = fig.add_axes([0.85, 0.2, 0.02, 0.6])
    cb1 = mpl.colorbar.ColorbarBase(cba,cmap=gmap,orientation='vertical')
    cb1.set_ticks(col)
    cb1.set_ticklabels([str(col[ivi]*100.0)+"%" for ivi in range(len(col))])
    if not diag=="":
        cb1.set_label("Strike probability - " + Tools.guess_diag(diag,True).nicename+" > "+str(int(diagthr))+Tools.guess_diag(diag,True).unit+" + Neighbour "+str(int(diagrad))+"km")

    return

def init_strike(ltraj,dom,res):
    #Initialize strike probability array
    #dom: domain [lonmin,lonmax,latmin,latmax]
    #res : resolution
    
    lonmin=dom[0]
    lonmax=dom[1]
    latmin=dom[2]
    latmax=dom[3]
    
    #Scan tracks for objects
    N_obj = 0 
    for traj in ltraj:
        N_obj = N_obj + traj.nobj
    
    #Define lists of longitudes and latitudes
    if lonmax<lonmin:
        lonmin=lonmin-180.0
    lons = []
    ivi = 0
    while lonmin+ivi*res<=lonmax:
        lons.append(lonmin+ivi*res)
        ivi=ivi+1
    lats = []
    ivi = 0
    while latmin+ivi*res<=latmax:
        lats.append(latmin+ivi*res)
        ivi=ivi+1
    
    height, width = len(lats), len(lons)
    image = np.full((height, width), 0)
    images = np.repeat(np.expand_dims(image, axis=0),
                                 N_obj,
                                 axis=0)
    
    return images, N_obj, lons, lats
    
def start_figure_map(**kwargs):
    #make a new map figure

    #add options according to the content of kwargs...

    #Create figures 
    fig = plt.figure(figsize=(16, 9))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # set graphical parameters
    # set land and sea
    land = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                        edgecolor='face',
                                        facecolor=cfeature.COLORS['land'])
    sea = cfeature.NaturalEarthFeature('physical', 'ocean', '10m',
                                        edgecolor='none',
                                        facecolor=cfeature.COLORS['water'])
                                        
    # set coastlines
    coastlines = cfeature.NaturalEarthFeature('physical','coastline', scale="10m",
                                              edgecolor="black", facecolor='none')

    #set country borders
    states_provinces = cfeature.NaturalEarthFeature(
        category='cultural',
        name='admin_1_states_provinces_lines',
        scale='50m',
        facecolor='none')
        
    #ax.add_feature(land)
    #ax.add_feature(sea)
    ax.add_feature(coastlines)
    #ax.add_feature(states_provinces, edgecolor = 'grey')
    if "epygram_departments" in kwargs.keys():
        draw_departments(ax,kwargs["epygram_departments"])  

    #resolution = '10m'
    #category = 'cultural'
    #name = 'admin_0_countries'
    #shpfilename = shapereader.natural_earth(resolution, category, name)
    #df = geopandas.read_file(shpfilename)
    #poly = df.loc[df['ADMIN'] == 'France']['geometry'].values[0]
    #ax.coastlines(resolution = '50m')


    # set latlon grid
    gl = ax.gridlines(draw_labels=True, linestyle = '--', linewidth = 0.5)
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    return fig, ax

def draw_departments(ax,epygram_departments):
    #based on epygram
    #(see ~/.epygram/src/epygram/_plugins/with_cartopy/H2DField.py)

    projection = ccrs.PlateCarree() #Attention! should be consistent with the map
    
    if not isinstance(epygram_departments, dict):
       epygram_departments = dict(color='k')
    import json
    with open(epygram.config.installdir + '/data/french_departments.json', 'r') as dp:
        depts = json.load(dp)[1]
    for d in list(range(len(depts)))[:]:
        for part in range(len(depts[d])):
            dlon = np.array(depts[d][part][0])
            dlat = np.array(depts[d][part][1])
            xyz = projection.transform_points(ccrs.PlateCarree(), dlon, dlat)
            x = xyz[..., 0]
            y = xyz[..., 1]
            ax.plot(x, y, **epygram_departments)

def set_dom_limits(ltraj,opt,diag="",diagrad=0.0,dom=[]):

    #Computes resolution
    print(ltraj[0].__dict__)
    if ltraj[0].inputdef["origin"]=="obs":
        leng = ltraj[0].inputdef["domain"]["lonmax"] - ltraj[0].inputdef["domain"]["lonmin"]
        if leng>100.0:
            res = 1.0
        elif leng>10.0:
            res = 0.25
        else:
            res = 0.1
    elif diag=="":
        par = ltraj[0].algodef["specfields"]["track"]
        if par in ltraj[0].algodef["parres"].keys():
            res = ltraj[0].algodef["parres"][par]
        else:
            res = ltraj[0].algodef["parres"]["all"]
    else:
        dd = Tools.guess_diag(diag,True)
        if dd.par in ltraj[0].algodef["parres"].keys():
            res = ltraj[0].algodef["parres"][dd.par]
        else:
            res = ltraj[0].algodef["parres"]["all"]

    if opt=="std":
        if ltraj[0].inputdef["origin"]=="obs":
            domtraj=ltraj[0].inputdef["domain"]
        else:
            domtraj=ltraj[0].algodef["domtraj"]
        lonmin=domtraj["lonmin"]
        lonmax=domtraj["lonmax"]
        latmin=domtraj["latmin"]
        latmax=domtraj["latmax"]
    elif opt=="zoom":
        print("zoom")
        llon=[]
        llat=[]
        if diag=="":
            for traj in ltraj:
                for obj in traj.traj:
                    llon.append(obj.lonc)
                    llat.append(obj.latc)
        elif dd.kind[0]=='c':
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
        if len(llon)>=2 and len(llat)>=2:
        #Compute length of grid and grid boundaries
            #print(opt,llon, llat)
            dyy = Tools.comp_length(np.mean(llon),np.mean(llat),np.mean(llon),np.mean(llat)+res)
            lonmin = res * int( (min(llon)/res - diagrad/dyy ) )
            lonmax = res * int( (max(llon)/res + diagrad/dyy ) ) + res
            latmin = res * int( (min(llat)/res - diagrad/dyy ) )
            latmax = res * int( (max(llat)/res + diagrad/dyy ) ) + res
        else:
            domtraj=ltraj[0].algodef["domtraj"]
            lonmin=domtraj["lonmin"]
            lonmax=domtraj["lonmax"]
            latmin=domtraj["latmin"]
            latmax=domtraj["latmax"]
    elif opt=="user" and len(dom)>0:
        lonmin=dom[0]
        lonmax=dom[1]
        latmin=dom[2]
        latmax=dom[3]  
    else:
        print("bad opt in set_dom_limits -- "+opt)
        exit()

    return lonmin,lonmax,latmin,latmax, res

def time_plot(ltraj, diag, typeplot="line", leg=[], opt = "", fig={}, ax=0, filename="", colormap="viridis", colorlevels=[], **kwargs):
    #ltraj, #list of tracks to take into account in the time plot
    #diag, #diagnostic name to plot
    #typeplot:, #can be “line” (default), “whisker”, “bar” or “dot”
    #fig, ax = figure name and axis on which to plot (optional input, default is 0)
    #leg=[], #legend for every track (if default, use member number) if relevant (“line” or “dot”),
    #colorlevels:[] # thresholds for different colors (default is []), valid if typeplot=“dot” or “bar”
    #colormap=, #colormap (used for leg, default is “viridis”)
    #“opt”: # optional, can be “max” (for maximum value along the forecast)
    #filename, #optional, output filename
    #kwargs #arguments passed to plot() routine
    
    #Plot definition
    if fig=={}:
        #Create new figure and 1 axis
        fig, ax0 = start_figure_time(**kwargs)

    lax = fig.get_axes()

    #print("LEN0",len(lax),len(fig.get_axes()))
    
    if ax==0:
        cax = lax[0]
        #print("LEN1",len(lax),len(fig.get_axes()))
    else:
        while len(lax)<ax+1:
            #print("LEN1",len(lax),len(fig.get_axes()))
            ax2 = lax[0].twinx()
            lax = fig.get_axes()
            #print("LEN2",len(lax),len(fig.get_axes()))            
        cax = lax[ax]
    #print("LEN",len(lax),len(fig.get_axes()))
    
    #Plotting
    if typeplot=="line":
        plot_time_line(ltraj, diag, cax, leg, opt, colormap, **kwargs)
    elif typeplot=="dot":
        plot_time_dot(ltraj, diag, fig, cax, leg, opt, colormap, colorlevels, **kwargs)
    elif typeplot=="whiskers":
        plot_time_whiskers(ltraj, diag, cax, leg, opt, colormap, **kwargs)
    elif typeplot=="bar":
        plot_time_bar(ltraj, diag, fig, cax, leg, opt, colormap, colorlevels, **kwargs)
    else:
        print("Wrong typeplot option in time_plot() : typeplot=",typeplot)
        exit()
    
    #ax.set_ylim([0,100])
    #ax2.set_ylim([0,100])
    plt.xlabel('Time')
    
    if len(filename) > 0:
        plt.savefig(filename)
    
    return fig

def start_figure_time(**kwargs):
    #make a new time plot
    
    # create figure
    fig = plt.figure(figsize=(16, 9))
    ax = plt.axes()
    
    # set graphical parameters
    ax.xaxis_date()
    myFmt = mdates.DateFormatter('%d/%H')
    ax.xaxis.set_major_formatter(myFmt)

    return fig, ax

def plot_time_line(ltraj, nam, cax, leg, opt, cmap, **kwargs):        
    #Plot individual members in ltraj as a function of time on cax
    #(called by time_plot)
    
    from matplotlib.colors import ListedColormap, LinearSegmentedColormap
    
    #read list of colors
    ntraj=len(ltraj)
    lcol=get_col(cmap,ntraj)
    
    lname = []
    for ivi in range(ntraj):
        time=[]
        pdiag=[]
        tra=ltraj[ivi]
        for ivj in range(tra.nobj):
            #print(tra.traj[ivi].__dict__['diags'], diag)
            if nam in tra.traj[ivj].__dict__['diags']:
                time.append(tra.traj[ivj].__dict__['time'])
                pdiag.append(tra.traj[ivj].__dict__[nam][2])
        
        timestamp = pd.to_datetime(time, format='%Y%m%d%H')
        cax.plot(timestamp,pdiag, color=lcol[ivi], **kwargs)
        if "member" in ltraj[ivi].inputdef.keys():
            lname.append(ltraj[ivi].inputdef["member"][1:3])
        else:
            lname.append(ltraj[ivi].basetime)
        
    #legend
    if len(leg) > 0:
        lname=leg

    if "ylim" in kwargs:
        cax.set_ylim(kwargs["ylim"])

    cax.set_ylabel(Tools.guess_diag(nam,True).nicename+" ("+Tools.guess_diag(nam,True).unit+") : lines")
    plt.legend(lname,ncol=3)
    
    return
    
def plot_time_dot(ltraj, nam, fig, cax, leg, opt, cmap, lev, **kwargs):  
    #Plot individual members in ltraj as dots on cax
    #if opt=="max", only the maximum value along the forecast is computed
    #(called by time_plot)
    
    #read list of colors
    if lev==[]: #a color different by member of ltraj
        ntraj=len(ltraj)
        lcol=get_col(cmap,ntraj)
        cbar=False
    elif len(lev)>1: # color depends on the diag level
        nlev=len(lev)
        gmap = mpl.cm.get_cmap(name=cmap,lut=nlev)
        lcol = [mpl.colors.to_hex(gmap(ivi), keep_alpha=True) for ivi in range(nlev)]
        #print(lcol)
        cbar = True
    else:
        print("Please include at least 2 colorlevels - ABORT")
        exit()
        
    #print(lcol)
        
    lname = []
    for ivi in range(len(ltraj)):
        time=[]
        pdiag=[]
        if opt=="max":
            pdiag.append(-999.0)
            time.append("")
        tra=ltraj[ivi]
        #print("Nobj : ",tra.nobj)
        for ivj in range(tra.nobj):
            #print(tra.traj[ivj].__dict__['diags'], nam)
            if nam in tra.traj[ivj].__dict__['diags']:
                val = tra.traj[ivj].__dict__[nam][2]
                tim = tra.traj[ivj].__dict__['time']
                if opt=="max":
                    if val >= pdiag[0]:
                        pdiag[0] = val
                        time[0] = tim
                else:
                    pdiag.append(val)
                    time.append(tim)
        #print(opt,pdiag,time)
        timestamp = pd.to_datetime(time, format='%Y%m%d%H')
        if lev==[]:
            cax.plot(timestamp,pdiag,'o', color=lcol[ivi]) #, **kwargs)
        else:
            #print(pdiag)
            #print(lev)
            for ivc in range(len(pdiag)):
                col=lcol[np.digitize(pdiag[ivc],lev)-1]
                #print(col)
                cc=cax.plot(timestamp,pdiag[ivc],'o', color=col)
        #print(ltraj[ivi].inputdef["member"])
        lname.append(ltraj[ivi].inputdef["member"][1:3])
        
    if "ylim" in kwargs:
        cax.set_ylim(kwargs["ylim"])

    #legend
    if len(leg)>0:
        lname=leg

    if cbar:
        add_colorbar(fig,gmap,lev,Tools.guess_diag(nam,True))
        #cba = fig.add_axes([0.93, 0.2, 0.02, 0.6])
        #cb1 = mpl.colorbar.ColorbarBase(cba,cmap=gmap,orientation='vertical')
        #cb1.set_ticks([float(ivi)/(nlev-1) for ivi in range(nlev)])
        #cb1.set_ticklabels(lev)
        #cb1.set_label(Tools.guess_diag(nam,True).nicename+" ("+Tools.guess_diag(nam,True).unit+") : dots")
    else:
        plt.legend(lname,ncol=3)
        cax.set_ylabel(Tools.guess_diag(nam,True).nicename+" ("+Tools.guess_diag(nam,True).unit+") : dots")    
    return

def plot_time_whiskers(ltraj, nam, cax, leg, opt, cmap, **kwargs):
    #Plot ensemble in ltraj as whiskers on cax
    #(called by time_plot)
    
    #Gets the list of instants to draw
    time=[]
    for ivi in range(len(ltraj)):
        tra=ltraj[ivi]
        for ivj in range(tra.nobj):
            if nam in tra.traj[ivj].__dict__['diags']:
                tim = tra.traj[ivj].__dict__['time']
                if tim not in time:
                    time.append(tim)
    timestamp = pd.to_datetime(time, format='%Y%m%d%H').sort_values()
    data = pd.DataFrame(index=timestamp,dtype="float32")
    
    ivj=0
    for tim in timestamp:
        tims=datetime.strftime(tim,time_fmt)
        #print(tims,tim)
        data.append([])
        #print(data)
        nobj=0
        for ivi in range(len(ltraj)):
            obj,found=ltraj[ivi].find_inst(tim)
            if found:
                data.loc[tim,nobj]=obj[0].__dict__[nam][2]
                nobj=nobj+1
                #print(ivj,nobj)
                #data[ivj].append(obj[0].__dict__[nam][2])
        ivj=ivj+1
    

    data.T.boxplot(ax=cax,positions=mdates.date2num(data.index),grid=False,widths=0.03)
    cax.set_xlim([data.index[0],data.index[-1]])
    #print(data.index.strftime("%d/%H"))
    cax.set_xticklabels( data.index.strftime("%d/%H") )
    if "ylim" in kwargs:
        cax.set_ylim(kwargs["ylim"])
    cax.set_ylabel(Tools.guess_diag(nam,True).nicename+" ("+Tools.guess_diag(nam,True).unit+") : whiskers")
    
    return

def plot_time_bar(ltraj, nam, fig, cax, leg, opt, cmap, lev, **kwargs):
    #Plot ensemble in ltraj as bars on cax
    #(called by time_plot)
    
    nlev=len(lev)
    
    #Gets the list of instants to draw
    time=[]
    for ivi in range(len(ltraj)):
        tra=ltraj[ivi]
        for ivj in range(tra.nobj):
            if nam in tra.traj[ivj].__dict__['diags']:
                tim = tra.traj[ivj].__dict__['time']
                if tim not in time:
                    time.append(tim)
    timestamp = pd.to_datetime(time, format='%Y%m%d%H').sort_values()
    data = pd.DataFrame(index=timestamp,dtype="float32")
    
    ivj=0
    for tim in timestamp:
        tims=datetime.strftime(tim,time_fmt)
        #print(tims,tim)
        data.append([])
        #print(data)
        nobj=0
        for ivi in range(len(ltraj)):
            obj,found=ltraj[ivi].find_inst(tim)
            if found:
                #data.loc[tim,nobj]=np.digitize(obj[0].__dict__[nam][2],colorlevels)-1
                data.loc[tim,nobj]=obj[0].__dict__[nam][2]
                nobj=nobj+1
        ivj=ivj+1
    
    #print(data)
    #print(data.index)

    il=0
    y=[]
    for ivi in range(len(data.index)):
        d=data.iloc[ivi] #list of values
        y.append(np.count_nonzero([d[ivj]>=lev[il] and d[ivj]<lev[il+1] for ivj in range(len(d))]))
    y0=[0 for ivj in range(len(y))]
    cax.bar(x=mdates.date2num(data.index),height=y,width=0.03,color=cmap(il))
        
    for il in range(1,nlev-1):
        for ivi in range(len(y)):
            y0[ivi]=y[ivi]+y0[ivi]
            y[ivi]=0
        #print(y)
        #print(y0)
        #exit()
        for ivi in range(len(data.index)):
            d=data.iloc[ivi] #list of values
            y[ivi]=np.count_nonzero([d[ivj]>=lev[il] and d[ivj]<lev[il+1] for ivj in range(len(d))])
        cax.bar(x=mdates.date2num(data.index),height=y,width=0.03,bottom=y0,color=cmap(il))

    if "ylim" in kwargs:
        cax.set_ylim(kwargs["ylim"])

    cba = fig.add_axes([0.93, 0.2, 0.02, 0.6])
    cb1 = mpl.colorbar.ColorbarBase(cba,cmap=cmap,orientation='vertical')
    cb1.set_ticks([float(ivi)/(nlev-1) for ivi in range(nlev)])
    cb1.set_ticklabels(lev)
    cb1.set_label("Number of members - " + Tools.guess_diag(nam,True).nicename+" ("+Tools.guess_diag(nam,True).unit+") : barplots")
    
    return


def get_col(cmap,n):
#returns list of n colors from the colormap cmap

    gmap = mpl.cm.get_cmap(name=cmap, lut = n)
    lcol=[]
    if n>1:
        for ivi in range(n):
            lcol.append(mpl.colors.to_hex(gmap(ivi/(n-1)), keep_alpha=True))
    elif n>0:
        lcol=[mpl.colors.to_hex(gmap(0), keep_alpha=True)]
        
    return lcol
    
def add_colorbar(fig,gmap,lev,dia):
#add colorbar to a plot

    nlev=len(lev)
    cba = fig.add_axes([0.93, 0.2, 0.02, 0.6])
    cb1 = mpl.colorbar.ColorbarBase(cba,cmap=gmap,orientation='vertical')
    cb1.set_ticks([float(ivi)/(nlev-1) for ivi in range(nlev)])
    cb1.set_ticklabels(lev)
    cb1.set_label(dia.nicename+" ("+dia.unit+") : dots")
    
    
