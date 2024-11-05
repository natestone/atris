import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import sys
import os

# I like colored printing:
STA = "\x1B["
STO = STA + "0m"
COL = ["31;40m", "36;40m", "32;40m", "35;40m"]

def aprint(text, color):
    """ custom colored printing using ascii escape sequences """
    print STA+COL[color] + "### AtRIS:" + text +STO


def extract_header(datapath):
    """ Simply extract the usefull data from the header """
    thefile = open(datapath,'r')
    bonusheader = ""
    status = 0
    nrrows = 0
    for k in thefile:
        if k == "\n":
            status = 0
            nrrows += 1
            continue
        if k[:6] == "1     ":
            nrrows+=1
            aprint("\tLast headerrow found: "+ str(nrrows),3)
            return nrrows
        kk = k.split()
        if kk[0] == "Selected":
            aprint("\tHeader found, now extracting column names",2)
            status = 1
        if kk[0] == "Year=" or "Latitude=" or "Prof.":
            bonusheader += k
        if status == 1:
            nr = kk[0]
            name = kk[1][:-1]
            unit = kk[2]
            #print nr, name, unit
        nrrows += 1



def nrml2psf(datapath):
    """ build the psf file
    Pressure is not part of the nrlmsise00 output. We therefore calculate it
    assuming ideal gas conditions. This is support by the relatively low
    densities which negate most deviations from ideal gas law. We also use
    Daltons law of partial pressures to calculate the net pressure. 
    """
    rows = extract_header(datapath)
    aprint("\tStarting processing... " + datapath, 2)
    data = np.loadtxt(datapath,skiprows=rows)
    alt  = data[:,0]
    temp = data[:,5]
    rho  = data[:,4]
    nN2  = data[:,2] + 0.5 * data[:,9] # molecular + atomic
    nO2  = data[:,3] + 0.5 * data[:,1] # molecular + atomic
    nHe  = data[:,6]
    nAr  = data[:,7]
    nH   = data[:,8]

    #Prepare output array
    nrSheets    = data.shape[0]           #one line per sheet
    output      = np.zeros((nrSheets,16)) #Array of zeros
    earthRadius = 6371000               #in meters
    thickness   = (data[1:,0]-data[:-1,0]) #uniform sheet thickness
    thickness  *= 1000                  #convert to meters
    print thickness
    thickness   = np.append(thickness,thickness[-1])
    aprint("\tAtmospheric thickness: "+str(thickness), 3)
    botAltitudes= earthRadius + alt*1000
    topAltitudes= botAltitudes + thickness

    #Fill out the pressures
    kb   = 1.38064852e-23 #m^2kgs^-2K^-1 - 
    pO2  = nO2 * kb * 1000000 * temp#convert to Pascals
    pN2  = nN2 * kb * 1000000 * temp#convert to Pascals
    pAr  = nAr * kb * 1000000 * temp#convert to Pascals
    pHe  = nHe * kb * 1000000 * temp#convert to Pascals
    pH   = nH  * kb * 1000000 * temp#convert to Pascals
    pTotal    = pO2+ pN2 + pAr + pHe + pH

    #Here are the molecular weights
    wO2  = 2*15.9994   
    wN2  = 2*14.0067   
    wAr  = 39.948      
    wHe  = 4.0026
    wH   = 1.008 
    u         = 1.66054e-24 #atomic mass unit in g
    mO2  = nO2*wO2*u #mass of oxygen in 1 cm3
    mN2  = nN2*wN2*u 
    mAr  = nAr*wAr*u 
    mHe  = nHe*wHe*u
    mH   = nH*wH*u
    total = mO2+mN2+mAr+mHe+mH
    mfO2 = mO2/total
    mfN2 = mN2/total
    mfAr = mAr/total
    mfHe = mHe/total
    mfH   = mH/total
    #for l in range(nrSheets):
    #    print "%10e ::: %.10f :: %.10f :: %.10f ::::: %.10f :: %.10f :: %.10f"%(nOxygen[l], mfOxygen[l],mfNitrogen[l],mfArgon[l],total[l], density[l],total[l]/ density[l])

    # we now fill out the output file
    output[::,0]  = botAltitudes / 1000.
    output[::,1]  = topAltitudes / 1000.
    output[::,2]  = 0 #dummy
    output[::,3] += 2*np.pi #unit is pi!
    output[::,4]  = 0 #dummy
    output[::,5] += 1*np.pi #unit is pi!
    output[::,6]  = np.arange(2,nrSheets+2)
    output[::,7]  = temp
    output[::,8]  = rho
    output[::,9]  = pTotal
    output[::,10]  = 0 #dummy
    output[::,11] = mfN2
    output[::,12] = mfO2
    output[::,13] = mfAr
    output[::,14] = mfHe
    output[::,15] = mfH
    
    aprint("\tCOMPLETE: writing out the gdml lines.",2)
    formatstring = "%e %e %e %e %e %e %i %e %e %e %e %e %e %e %e %e"
    #formatstring = "%i %i %.2f %.2f %.2f %.2f %i %.1f %e %e %.5f %.5f %.5f %.5f %.5f %.5f"
    np.savetxt(arg1+".psf",output,fmt=formatstring)
    # we now read the saved file an add the header, the core, crust and loss
    # volumes. 
    f      = open(arg1+".psf","r")
    ff     = f.readlines()
    f.close()

    cr = 20/1000. #thickness in kilometers
    er = earthRadius/1000.
    top= (topAltitudes[-1]+1)/1000.
    t  = 293.15 #temperature in Kelvins
    p  = 101325 #pressure in Pascals
    d  = 2.7    #density in g/cm3
    f  = open(arg1+".psf","w")
    f.write("#O478Si284Al84Fe52Ca37Na28Mg15 N O Ar He H\n") # Adding the header file
    f.write(formatstring%(0,er-cr,0,2*np.pi,0,1*np.pi,0,t,d,p,1.0,0,0,0,0,0))
    f.write("\n")
    f.write(formatstring%(er-cr,er,0,2*np.pi,0,1*np.pi,1,t,d,p,1.0,0,0,0,0,0))
    f.write("\n")
    for k in ff: f.write(k)
    f.write(formatstring%(top,top+cr,0,2*np.pi,0,1*np.pi,nrSheets+2,t,d,p,1.0,0,0,0,0,0))
    f.write("\n")
    f.close()




if __name__ == "__main__":
    os.system("clear")
    aprint("\t################################################################",1)
    aprint("\t$$                  NRLMSISE00 2 PSF                          $$",1)
    aprint("\t################################################################",1)
    aprint("\tThis is a simple interface to the NRLMSISE00 model data.",2)
    aprint("\tIt generates psf files from the ascii file that can be obtained", 2)
    aprint("\ton https://ccmc.gsfc.nasa.gov/modelweb/models/nrlmsise00.php.", 2)
    aprint("\tBecause of the unusual order of the nrlmsise data file, it is",2)
    aprint("\trather complicated to build a do-it-all interace. Use this as",2)
    aprint("\tan example. Two screenshots sett1.png and sett2.png show what",2)
    aprint("\tsettings have been used to generate compatible data. Usage",2)
    aprint("\t",2)
    aprint("\t\tpython nrlm2psf.py filenameprefix",2)
    aprint("\t",2)
    aprint("\twhere filenameprefix is the filename without the extension .nrlm",2)
    aprint("\t################################################################",1)
    aprint("\t################################################################",1)

    # load data, count the abundance and save the data:
    arg1 = sys.argv[1]
    if os.path.isfile(arg1+'.nrlm'):
        aprint('\tNow processing the ascii file: '+arg1+".nrlm",3)
        nrml2psf(arg1+'.nrlm')
    else:
        aprint('\tThe file ' + arg1 + '.nrlm does not exist.',0)
