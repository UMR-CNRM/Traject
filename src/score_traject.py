# -*- coding: utf-8 -*-
"""
Some Score functions associated to Traject

https://github.com/UMR-CNRM/Traject
Traject is disributed under Cecill-C license
(For rights and obligations, see LICENSE.txt)

UMR-CNRM, Meteo-France & CNRS

"""

import os
import numpy as np
from traject import *
import Tools
import pandas as pd
import matplotlib as mpl
import matplotlib.dates as mdates
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import pyproj


#---------------------------------------------------------
def score_comp_diff(lref, ltraj, basetime, ldiag, nmb, fname, proj="merc"):
    #Computes differences between forecast tracks and of reference tracks
    #We assume that the names in lref and in ltraj correspond to the same tracks
    #lref : list of reference tracks
    #ltraj : list of forecast tracks
    #basetime : list of basetimes
    #ldiag : list of diagnostics computed (ate, tte, cte or other name in the forecast tracks)
    #        If several diagnostic names in traj or ref correspond to a ldiag name, only the first one is taken into account !!
    #nmb : number of members
    #proj : projection used for computation of tte, cte, ate (default=mercator)
    #output : diagnostics are written in fname file (csv format)

    #Structure of the DataFrame:
    #to be completed

    hline = ["name","basetime","fctime","diag","obs"]
    for im in range(nmb):
        hline.append("mb"+str(im).rjust(3,"0"))
    ncol=len(hline)
    print(hline)
    data = pd.DataFrame(columns=hline,dtype="float32")
    print(data)
    #print(data.loc(['name']))

    #Correspondance of diagnostics names between ltraj and lref
    namediagr={}
    namediagt={}
    diagr=lref[0].traj[0].diags
    diagt=ltraj[0].traj[0].diags    
    Hn=lref[0].traj[0].latc > 0.0
    notdiag=['ATE','ate',"CTE","cte","tte","TTE","lonc","latc"]
    ldiag2=[dg for dg in ldiag if not dg in notdiag]
    for dg in ldiag2:
        d1=Tools.guess_diag(dg,Hn)
        for dr in diagr:
            d2=Tools.guess_diag(dr,Hn)
            chkis=False
            if d2.compares(d1) and not chkis:
                namediagr[dg]=dr
                chkis=True
        for dt in diagt:
            d2=Tools.guess_diag(dt,Hn)
            chkis=False
            if d2.compares(d1) and not chkis:
                namediagt[dg]=dt
                chkis=True
    print(namediagr)
    print(d1.__dict__)
    print(namediagt)
    print(d2.__dict__)
                        
    #Loop on tracknames (in reference)
    for tref in lref:
        tname = tref.name

        #Loop on basetime
        for bt in Inputs.comptimes(basetime):

            #Loop on forecast tracks
            for traj in Select(ltraj,{"basetime":bt,"name":tname}):

                #Loop on forecast lead time of traj
                for it in range(traj.nobj):
                    term = Tools.comp_difftime(bt,traj.traj[it].time)

                    #Search for reference object that corresponds to this lead time
                    itref, chk = tref.find_ind(traj.traj[it].time)

                    if chk:
                        itr = itref[0]

                        for diag in ldiag:
                            dg = diag.lower()
                            ind = get_index(tname,bt,term,dg)
                            if ind not in data.index: #Create line
                                dline=[tname,bt,term,dg]
                                dline.append(np.nan) #obs
                                for im in range(nmb):
                                    dline.append(np.nan)
                                data.loc[get_index(tname,bt,term,dg)] = dline
                            #print(data.index)

                            mb=0 #member
                            if "member" in traj.inputdef:
                                mb=int(traj.inputdef["member"])

                            if dg=="tte":
                                print("Compute ", dg)
                                #itref, chk0 = tref.find_ind(traj.traj[it].time)
                                ldist = comp_tte_cte_ate(tref,itr,traj.traj[it],proj)
                                data.loc[ind,"mb"+str(mb).rjust(3,"0")] = ldist[dg]
                                data.loc[ind,"obs"] = 0.0
                                #dline.append(ldist[dg.lower()])
                            elif dg=="ate" or dg=="cte":
                                print("Compute ", dg)
                                if it>0:
                                    itrefm1, chk1 = tref.find_ind(traj.traj[it-1].time)
                                else:
                                    chk1 = False
                                if it<traj.nobj-1:
                                    itrefp1, chk2 = tref.find_ind(traj.traj[it+1].time)
                                else:
                                    chk2 = False
                                if chk1 and chk2:
                                    ldist = comp_tte_cte_ate(tref,itr,traj.traj[it],proj,dt=[itrefm1[0],itrefp1[0]])
                                    data.loc[ind,"mb"+str(mb).rjust(3,"0")] = ldist[dg]
                                    data.loc[ind,"obs"] = 0.0
                                #dline.append(ldist[dg.lower()])
                            elif dg=="lonc":
                                data.loc[ind,"mb"+str(mb).rjust(3,"0")] = traj.traj[it].lonc - tref.traj[itr].lonc
                                data.loc[ind,"obs"] = tref.traj[itr].lonc
                            elif dg=="latc":
                                data.loc[ind,"mb"+str(mb).rjust(3,"0")] = traj.traj[it].latc - tref.traj[itr].latc
                                data.loc[ind,"obs"] = tref.traj[itr].latc
                            else: #other diagnostic
                                print(traj.traj[it].__dict__[namediagt[dg]], tref.traj[itr].__dict__[namediagr[dg]])
                                lvalt=traj.traj[it].__dict__[namediagt[dg]]
                                lvalr=tref.traj[itr].__dict__[namediagr[dg]]
                                if isinstance(lvalt,float):
                                    valt=lvalt
                                else:
                                    valt=lvalt[2]
                                if isinstance(lvalr,float):
                                    valr=lvalr
                                else:
                                    valr=lvalr[2]                                
                                if not (valt==missval or valr==missval):
                                    data.loc[ind,"mb"+str(mb).rjust(3,"0")] = valt - valr
                                    data.loc[ind,"obs"] = valr

    # Write output file
    data.to_csv(fname, sep=',')
    
    return data

def get_index(name,bt,term,diag):
    #Computes the unique index associated to name, basetime, term, and diagnostic

   return diag.upper()+"-"+name+"-"+bt+"-"+str(term).rjust(3,"0")

def LonLat_To_XY_merc(proj,Lon,Lat):

    P_proj = pyproj.Proj(proj=proj, ellps='WGS84') #,  units='m')

    return P_proj(Lon,Lat)

def XY_To_LonLat_merc(proj,x,y):

    P_proj = pyproj.Proj(proj=proj, ellps='WGS84') #,  units='m')

    return P_proj(x,y,inverse=True)

def distance_plan(x1,y1,x2,y2):
    return np.sqrt((x2-x1)**2+(y2-y1)**2)

def sinus (A, B, dA, dB) :
    return (A[0]*B[1]-A[1]*B[0])/(dA*dB)


def comp_tte_cte_ate(tref,itref,obj,proj,dt=[]):
    #Compute tte, cte, ate of object obj
    #from the tref reference track (index itref) ; direction is deduced from (itref-1) to (itref)
    #proj : projection used for plan distance ("merc" for Mercator)
    #dt: direction for computation of cte and ate (not used for tte)
    #Routine by C. Labadie and colleagues

    # calcul TTE
    TTE=Tools.comp_length(tref.traj[itref].lonc,tref.traj[itref].latc,obj.lonc,obj.latc)
    
    #Calcul CTE et ATE
    if len(dt)==0:
        
        ATE=missval * 1000.0
        CTE=missval * 1000.0
        
    else:
    
        # définition vecteurs AP et trajectoire + leur norme
        xA, yA = LonLat_To_XY_merc(proj,tref.traj[itref].lonc, tref.traj[itref].latc)
        xP, yP = LonLat_To_XY_merc(proj,obj.lonc, obj.latc)

        AP=np.array([xP-xA,yP-yA])
        dAP=distance_plan(xA,yA,xP,yP)

        print(dt)
        x_A_t0, y_A_t0 = LonLat_To_XY_merc(proj,tref.traj[dt[0]].lonc, tref.traj[dt[0]].latc)
        x_A_t2, y_A_t2 = LonLat_To_XY_merc(proj,tref.traj[dt[1]].lonc, tref.traj[dt[1]].latc)
        trajA=np.array([x_A_t2-x_A_t0, y_A_t2-y_A_t0])
        dtrajA=distance_plan(x_A_t0,y_A_t0, x_A_t2, y_A_t2)

        print("dtraj : ", dtrajA)

        # test valeur de TTE (pour éviter que ATE et CTE = nan qd TTE est nul)
        if TTE > 1e-6:

            print("TTE=",TTE)  

            # calcul cosinus et sinus de l'angle entre la trajectoire analysée et le vecteur AP

            cos_alpha = float(np.vdot(AP,trajA)/(dAP*dtrajA))
   
            sin_alpha = float(sinus(trajA, AP, dtrajA, dAP))
            #print("sin alpha : ", sin_alpha, type(sin_alpha))

            # calcul ATE, CTE
            ATE=cos_alpha*TTE
            CTE=sin_alpha*TTE
 
        else :
            ATE = TTE
            CTE = TTE
  
    return {"tte":TTE,"ate":ATE,"cte":CTE}

#######################
#### Plotting routines
#######################

def start_figure_score():
    #make a new score plot
    
    # create figure
    fig = plt.figure(figsize=(16, 9))
    ax = plt.axes()
    
    return fig, ax
    
def plot_boxplot_score(fig, df, diag, pdt, echmax, qmin, qmax, **kwargs) :
  #Boxplot of diag (can be "TTE" or "ATE" or "CTE" or other diagnostic that is present in df)
  #df : dataframe generated by score_comp_diff()
  
    if isinstance(df,str):
        data = pd.read_csv(df)
    else:
        data = df

    #List of members in the dataframe
    lcol = list(data.columns)
    lmembers = [col for col in lcol if col[0:2]=="mb"]
    
    #Filtering diagnostic name
    datid = data.loc[data["diag"]==diag.lower()]

    fct=1
    if diag.lower() not in ["tte","cte","ate"]:
        fct = Tools.guess_diag(diag,True).plot_fct
    
    # Boucle sur les différentes échéances
    l_TE = []
    l_nb_parech = []
    l_ech = []
    ech = 0
    while ech<=echmax:
        datech = datid.loc[datid["fctime"]==int(ech)]
        vals = []
        for mb in lmembers:
            vals.extend([x*fct for x in datech[mb] if not np.isnan(x)])
        nb_parech=len(vals)/len(lmembers)
        l_TE.append(vals)
        l_ech.append(ech)
        l_nb_parech.append(nb_parech)
        ech = ech + pdt

    #Plot properties
    plt.rcParams['font.size'] = 12
    lwidth = 2
    boxprops = dict(linestyle = '-', linewidth = lwidth, color = 'black')
    whiskerprops = dict(linestyle = '-', linewidth = lwidth, color = 'black')
    flierprops = dict(marker = 'D', markerfacecolor = 'k', markersize = 5,
                  markeredgecolor = 'k')
    medianprops = dict(linestyle = '-', linewidth = lwidth, color = 'tab:red')
    meanpointprops = dict(marker = 'o', markeredgecolor = 'tab:blue',
                      markerfacecolor = 'tab:blue', markersize = 5)

    #Make plot
    #plt.boxplot(l_TE, medianprops= {'color' : 'b'},whis=[int(qmin),int(qmax)], showfliers=False)
    plt.boxplot(l_TE,
            whiskerprops = whiskerprops,
            boxprops = boxprops,
            medianprops = medianprops,
            meanprops = meanpointprops,
            whis=[int(qmin),int(qmax)],
            showfliers=False,
            showmeans=True,
            zorder=2)

    ax=plt.gca()
    ax.set_xticklabels(labels=l_ech,rotation=45)
    ax.set_xlabel("Echéance (en h)")#, ha='left')
    if "ylim" in kwargs:
        ax.set_ylim(kwargs["ylim"])
    if diag.lower() in ["tte","cte","ate"]:
        ylab=diag+" (km)"
    else:
        ylab=Tools.guess_diag(diag,True).nicename+" ("+Tools.guess_diag(diag,True).plot_unit+") : lines"
    ax.grid(axis = 'y')
    ax.set_ylabel(ylab)

    # ajout de la ligne 0
    ax.plot([0]*(len(l_TE)+2),'--', color='lightgrey')

    # Ajout de la courbe : nb de cyclones par ech
    ax1=ax.twinx()
    x_decal = list(range(1,len(l_ech)+1)) #decalage car boxplot commence a 1 ...
    ax1.plot(x_decal,l_nb_parech,'-.',color='grey',zorder=1)#, label="nb cyclone")
    ax1.set_ylim([0,max(l_nb_parech)+10])
    ax1.set_ylabel("Number of tracks")
    
    return


def plot_line_score(fig, df, diag, pdt, echmax, lmetric, **kwargs) :
  #Linear plots of diag for several experiments (can be "TTE" or "ATE" or "CTE" or other diagnostic that is present in df)
  #ldf : list dataframe generated by score_comp_diff() (the different experiments)
  #lmetric : list of values that are plotted (can be "nb" (number of forecasts per term), "mean", "max", "min" or a quantile (q50, q90...))
  #by default, lmark is of different line options

    if fig=={}:
        #Create new figure and 1 axis
        fig, ax = start_figure_score()
    else:
        ax=plt.axes()

    #read input data
    if isinstance(df,str):
        data = pd.read_csv(df)
    else:
        data = df

    if isinstance(lmetric,str):
        lmetric1=[lmetric]
    else:
        lmetric1=lmetric

    #read plot options
    if "lmark" in kwargs:
        lmark=kwargs["lmark"]
    else:
        lmark=["-","o","."]

    #List of members in the dataframe
    lcol = list(data.columns)
    lmembers = [col for col in lcol if col[0:2]=="mb"]

    #Filtering diagnostic name
    datid = data.loc[data["diag"]==diag.lower()]

    fct=1
    if diag.lower() not in ["tte","cte","ate"]:
        fct = Tools.guess_diag(diag,True).plot_fct

    # Boucle sur les différentes échéances
    l_TE = []
    l_nb_parech = []
    l_ech = []
    ech = 0
    while ech<=echmax:
        datech = datid.loc[datid["fctime"]==int(ech)]
        vals = []
        for mb in lmembers:
            vals.extend([x*fct for x in datech[mb] if not np.isnan(x)])
        nb_parech=len(vals)/len(lmembers)
        l_TE.append(vals)
        l_ech.append(ech)
        l_nb_parech.append(nb_parech)
        ech = ech + pdt

   #Plot
    for met in lmetric1:

        plval=[]
        if met=="nb": #Plot nombre de forecasts parech
            plval=l_nb_parech
        elif met=="mean":
            for ivi in range(len(l_ech)):
                if len(l_TE[ivi])>0:
                    plval.append(np.mean(l_TE[ivi]))
                else:
                    plval.append(np.nan)
        elif met=="min":
            for ivi in range(len(l_ech)):
                if len(l_TE[ivi])>0:
                    plval.append(np.min(l_TE[ivi]))
                else:
                    plval.append(np.nan)
        elif met=="max":
            for ivi in range(len(l_ech)):
                if len(l_TE[ivi])>0:
                    plval.append(np.max(l_TE[ivi]))
        elif met[0]=="q":
            for ivi in range(len(l_ech)):
                if len(l_TE[ivi])>1:
                    plval.append(np.quantile(l_TE[ivi],float(met.replace("q",""))/100.0))
                else:
                    plval.append(np.nan)
        else:
            print("Bad metric in plot_line_score()")
            plval = [0.0 for ivi in range(len(l_ech))]

        ax.plot(l_ech,plval,**kwargs)

        ax=plt.gca()
        ax.set_xticks(l_ech)
        ax.set_xticklabels(labels=l_ech,rotation=45)
        ax.set_xlabel("Echéance (en h)")#, ha='left')

    return


