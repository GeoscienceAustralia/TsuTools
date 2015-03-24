#!/usr/bin/python

"""
Script to plot surface deformation from a gmt grd file

INPUT: 
grd_file

Nick Horspool
Geoscience Australia
Feb 2012
"""

import subprocess, sys, os
from subprocess import call, Popen, PIPE, STDOUT

def usage():
    return 'python plot_deformation.py <grd_file>'

def plot_deformation(filename):
    """ Call GMT commands to plot surface deformation
    """

    call(["gmtset","BASEMAP_TYPE","fancy","LABEL_FONT_SIZE","8p","ANNOT_FONT_SIZE","8p","HEADER_FONT_SIZE","10p","D_FORMAT","%lg"])    
        
    p = subprocess.Popen(["grdinfo","-I-",filename], stdout=PIPE)
    for line in p.stdout:
        region = line.strip()
#    region = '-R91.7/94.1/11.4/13.5'
    proj="-JM15c"     
    plot = open(filename + '.ps','w')   

    cmd = ('grdinfo -C %s > grdinfo.tmp' % filename)
    os.system(cmd)
    
    #p = subprocess.Popen(['''awk '{print "-"$7"/"$7"/"0.01}' grdinfo.tmp'''], stdout=PIPE)
    p = subprocess.Popen([r"awk","{print -$7,$7,0.05}","grdinfo.tmp"], stdout=PIPE)
    for line in p.stdout:
        z_range = line.strip()
    print z_range
    
    t1=z_range.split()[0]
    t2=z_range.split()[1]
    tz=z_range.split()[2]
    
    #Force range
    t1='-5.0'
    t2='5.0'
    tz = '0.5'
    T_range="-T"+t1+"/"+t2+"/"+tz
    print 'Warning forcing range and increment manually!!!'
    print T_range
    
    cpt_file = open("defm.cpt",'w')
    call(["makecpt","-D","-Cpolar",T_range,"-Z"],stdout=cpt_file)

    call(["grdimage","-Ba2f1",filename,region,proj,"-Cdefm.cpt","-P","-K","-Y5c"],stdout=plot)

    call(["pscoast",region,proj,"-Dh","-W0.5p/0","-O","-K","-P"],stdout=plot)

    call(["psscale","-Ba1.0f0.5:Co-seismic Deformation (m):","-D7.5/-1/8/0.5h","-Cdefm.cpt","-O","-P"],stdout=plot)

    call(["ps2raster",filename+".ps","-Tg","-A"])

if __name__ == '__main__':

    # Input checks
    if len(sys.argv) < 2:
        print usage()
        sys.exit()
    filename = sys.argv[1]  
    plot_deformation(filename)
