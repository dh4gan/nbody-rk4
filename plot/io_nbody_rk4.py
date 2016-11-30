# -*- coding: utf-8 -*-
# Written 25/11/16 by dh4gan
# Module contains methods for I/O relating to nbody_rk4 
# Includes: methods to read in output files (single body, snapshot, log)
#        Future methods to include ability to write setup files


import numpy as np
import matplotlib.pyplot as plt
import filefinder as ff
from orbitalelements import orbitalElements
from vector import Vector3D

# Declare some useful dictionaries for plotting here

# single particle files 
single_keys = ['t', 'mass', 'x', 'y', 'z', 'vx', 'vy', 'vz','ax','ay','az','semimaj','ecc','inc','longascend','argper','trueanom','ekin','epot','etot','angmomx','angmomy','angmomz']
single_labs = ['Time (yr)', 'Mass ($M_{\odot}$)','x (AU)', 'y  (AU)', 'z  (AU)', '$v_x$  (AU/yr)', '$v_y$ (AU/yr)', '$v_z$ (AU/yr)','$a_x$ (AU/yr)','$a_y$ (AU/yr)','$a_z$ (AU/yr)','$a$ (AU)','$e$','$i$','$\Omega$','Argument of Periapsis','True Anomaly','$E_{kin}$','$E_{pot}$','$E_{tot}$','$L_x$','$L_y$','$L_z$']
single_cols = range(len(single_keys))

singlecol = dict(zip(single_keys,single_cols))
singlelabel = dict(zip(single_keys,single_labs))

# snapshot files

snap_keys = ['mass','x', 'y', 'z', 'vx', 'vy', 'vz','ax','ay','az','semimaj','ecc','inc','longascend','argper','trueanom','ekin','epot','etot','angmomx','angmomy','angmomz']
snap_labs = ['Mass ($M_{\odot}$)','x (AU)', 'y  (AU)', 'z  (AU)', '$v_x$  (AU/yr)', '$v_y$ (AU/yr)', '$v_z$ (AU/yr)','$a_x$ (AU/yr)','$a_y$ (AU/yr)','$a_z$ (AU/yr)','$a$ (AU)','$e$','$i$','$\Omega$','Argument of Periapsis','True Anomaly','$E_{kin}$','$E_{pot}$','$E_{tot}$','$L_x$','$L_y$','$L_z$']
snap_cols = range(len(snap_keys))

snapcol = dict(zip(snap_keys,snap_cols))
snaplabel = dict(zip(snap_keys,snap_labs))

# log files

log_keys = ['t','dt','error', 'E', 'dE', 'L', 'dL']
log_labs = ['Time (yr)','Timestep (yr)','Error (% of tolerance)', 'Total Energy', 'dE', 'Angular Momentum', 'dL']
log_cols = range(len(log_keys))

logcol = dict(zip(log_keys,log_cols))
loglabel = dict(zip(log_keys,log_labs))

# G in units (AU,days, Msol)

Gmau = (2.0*np.pi/365.24)*(2.0*np.pi/365.24)

# One Earth Radius in AU
earthradiitoAU = 4.258e-5 

def getFileLineCount(filename):
    
    with open(filename) as f_obj:
        count = sum(1 for line in f_obj)
    return count

def select_variables(filetype='snapshots'):        

    if(filetype=='single'):
        variablekeys = single_keys
        coldict = singlecol
        labeldict= singlelabel
    elif(filetype=='log'):
        variablekeys = log_keys
        coldict = logcol
        labeldict = loglabel
    else:
        variablekeys = snap_keys
        coldict = snapcol
        labeldict= snaplabel
        
    print "Which variable for the x-axis?"

    for i in range(len(variablekeys)):
        print variablekeys[i],":\t \t \t", labeldict[variablekeys[i]]

    keyword = raw_input("Enter appropriate keyword:   ")

    xkey = keyword
    ix = coldict[keyword]
    xlabel = labeldict[keyword]

    print "Which variable for the y-axis?"

    for i in range(len(variablekeys)):
        if(i!=ix): print variablekeys[i],":\t \t", labeldict[variablekeys[i]]

    keyword = raw_input("Enter appropriate keyword:   ")

    ykey = keyword
    iy = coldict[keyword]
    ylabel = labeldict[keyword]
        
    return xkey, ykey, ix,iy,xlabel,ylabel

def read_body_file(filename):
    
    print 'Reading single particle data from file ',filename
    return np.genfromtxt(filename)

def read_log_file(filename):

    print 'Reading log data from file ',filename        
    return np.genfromtxt(filename)

def read_snapshot_file(filename):
    print 'Reading snapshot data from file ',filename
    
    f_obj = open(filename,'r')
    time = float(f_obj.readline())
    
    data = np.genfromtxt(filename,skiprows=1)
    return time, data

def read_all_bodyfiles(prefix):
        
    filenames = ff.find_sorted_local_input_fileset(prefix)
    
    # Ignore log files
    for i in range(len(filenames)):
        if '.log' in filenames[i]:
            print 'Removing log file ',filenames[i] 
            filenames.remove(filenames[i])
    
    nbodies = len(filenames)
    
    print 'Reading data for ',nbodies, ' bodies'        
              
    data = []      
    for i in range(nbodies):
        filename = filenames[i]
        print 'Reading ',filename        
        data.append(read_body_file(filename))
    
    return filenames, data

def plot_body(filename, showPlot=False):
    
    filetype = 'single'    
    xkey, ykey, ix,iy, xlabel,ylabel = select_variables(filetype)    
    data = read_body_file(filename)
    
    outputfile = ykey+'_vs_'+xkey+'_'+filename+'.png'

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
            
    ax1.plot(data[:,ix], data[:,iy])
    ax1.set_xlabel(xlabel, fontsize = 22)
    ax1.set_ylabel(ylabel, fontsize = 22)
    
    if(showPlot): 
        plt.show()
    plt.savefig(outputfile)


def plot_body_allfiles(prefix,showPlot=False):
    
    filetype = 'single'    
    xkey, ykey, ix,iy,xlabel,ylabel = select_variables(filetype)
    
    outputfile = 'multiparticle_'+ykey+'_vs_'+xkey+'.png'
    
    filenames,alldata = read_all_bodyfiles(prefix)
        
    nbodies = len(alldata)
    
    print nbodies
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    
    ax1.set_xlabel(xlabel, fontsize = 22)
    ax1.set_ylabel(ylabel, fontsize = 22)
        
    for i in range(nbodies):            
                
        try:
            ax1.plot(alldata[i][:,ix], alldata[i][:,iy], label=filenames[i])
        except IndexError:
            print "File ",filenames[i], " has no entries, skipping"
    
    ax1.legend()
    if(showPlot): 
        plt.show()
    plt.savefig(outputfile)
    
    
def plot_snapshot(filename,showPlot=False):
    filetype = 'snapshot'   
    
    xkey, ykey, ix,iy,xlabel,ylabel = select_variables(filetype)
    
    time, data = read_snapshot_file(filename)
    
    outputfile = ykey+'_vs_'+xkey+'_'+filename+'.png'

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
            
    ax1.scatter(data[:,ix], data[:,iy], marker = '.', size = np.sqrt(data[:,snapcol['mass']]))
    ax1.set_xlabel(xlabel, fontsize = 22)
    ax1.set_ylabel(ylabel, fontsize = 22)
    
    ax1.text(0.9, 0.9,'t={:.3E} yr'.format(time), bbox=dict(edgecolor='black',facecolor='none'), horizontalalignment='center', verticalalignment='center',transform = ax1.transAxes)
    
    if(showPlot): 
        plt.show()
    plt.savefig(outputfile)

    
def animate_bodies(prefix):    
    
    filetype = 'single'
    xkey, ykey, ix,iy,xlabel,ylabel = select_variables(filetype)
    
    imass = singlecol['mass']
    itime = singlecol['t']
    
    
    filenames, alldata = read_all_bodyfiles(prefix)   
        
    nbodies = len(filenames)
    
    plt.ion()
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
                
    
    nsteps = len(alldata[0])
 
    xmax = -1.0e30
    xmin = -xmax
    
    ymin = 1.0e30
    ymax = -ymin
       
    for i in range(nbodies):
        body_xmax = np.amax(alldata[i][:,ix])
        body_ymax = np.amax(alldata[i][:,iy])
        
        body_xmin = np.amin(alldata[i][:,ix])
        body_ymin = np.amin(alldata[i][:,iy])
        
        if(body_xmax>xmax): xmax = body_xmax
        if(body_ymax>ymax): ymax = body_ymax

        if(body_xmin<xmin): xmin = body_xmin        
        if(body_ymin<ymin): ymin = body_ymin
    
    
    print 'Plot range: '
    print xlabel,' : ',xmin, xmax
    print ylabel,' : ',ymin, ymax
    
                 
    iskip = np.zeros(nbodies)
    for j in range(nsteps):
       
        for i in range(nbodies):
            
            if(iskip[i]==1):
                continue
            
            try:
                t = alldata[i][j,itime]
            except IndexError:
                print 'No further data for ',filenames[i], ': skipping'
                iskip[i]=1
                continue
                                   
            if(alldata[i][j,imass]>0.0):
                
                ax1.scatter(alldata[i][j,ix], alldata[i][j,iy], s = np.sqrt(1.0e8*alldata[i][j,imass]))
            else:
                ax1.scatter(alldata[i][j,ix], alldata[i][j,iy], s = 5.0, color='k')
                
                        
        plt.savefig('animation_0'+str(j)+'.png', format='png')
        ax1.text(0.9, 0.9,'t={:.3E} yr'.format(t), bbox=dict(edgecolor='black',facecolor='none'), horizontalalignment='center', verticalalignment='center',transform = ax1.transAxes)
        
        ax1.set_xlabel(xlabel, fontsize = 22)
        ax1.set_ylabel(ylabel, fontsize = 22)
        ax1.set_ylim(ymin,ymax)
        ax1.set_xlim(xmin,xmax)
        #ax1.set_xscale('log')
        #ax1.set_yscale('log')
        plt.draw()
    
        ax1.clear()


def calc_orbit_ellipse(orbitdata,itime):
    
    a = orbitdata[itime,singlecol['a']]
    ecc = orbitdata[itime,singlecol['e']]
    inc = orbitdata[itime,singlecol['i']]
    argper = orbitdata[itime,singlecol['peri']]
    longascend = orbitdata[itime,singlecol['node']]
        
    # Find total mass of the system TODO
    
    totalmass = 1.0
    
    orbit = orbitalElements(a,ecc,inc,argper,longascend, 0.0, Vector3D(0.0,0.0,0.0),Vector3D(0.0,0.0,0.0),1.0,totalmass)
    
    npoints = 100
    x,y,z = orbit.calcOrbitTrack(npoints)
    
    return x,y,z
    
    