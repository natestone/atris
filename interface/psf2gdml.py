import numpy as np
import matplotlib.pyplot as plt
import re
import os
import periodic
import sys 

# We load "nist.exist", which holds information about which materials already
# exist in the Geant4 material database. This information is useful when 
# creating sheet materials. 
f      = open('nist.exist','r')
nistdb = f.readline()
# Appending own materials. Here, I should add all the project-related molecules
nistdb = nistdb[:-1]+",'G4_VACUUM']"
f.close()

                                                                                                                                                                                                                                                            
# I like colored printing:
STA = "\x1B["
STO = STA + "0m"
COL = ["31;40m", "36;40m", "32;40m", "35;40m"]

def aprint(text, color):
    """ custom colored printing using ascii escape sequences """
    print STA+COL[color] + "### AtRIS:" + text +STO


def parse_solid(formula):
    """
    parse_solid is ment to deal with the first element in the header.
    In the PSF format spec, we have described in detail the representation
    of the header. To summarize, the first element is a series of chemical
    elements followed by their mass fraction in per milles. For example,
    O762C111H101N26 corresponds to:
        : 762 parts oxygen
        : 111 parts carbon
        : 101 parts hydrogen
        :  26 parts nitrogen
        ====================
        :1000 parts
    Notes:
        - Sum of per milles should add up to 1000
        - Order according to decreasing abundance
        - This is only for the header! The GDML material will be designated
          as a "solid"
        - The macroscopic properties of the core elemnts still have to be 
          provided in the first line of the PSF file. 
    """
    # convert formula using regular expressions in to a list of tuples
    first = re.findall(r'([A-Z][a-z]*)(\d*)', formula)
    for k in first:
        element = k[0]
    if element not in nistdb: 
        print "element not found in the nist database. ABORTING!!"
        return 0
    g4mat, totalmass = {}, 0.
    for k in first:
        material    = k[0]
        massfrac    = float(k[1]) / 1000.
        g4mat.update({material:{"mass fraction":massfrac, "material":material}})
        totalmass   += massfrac        
    if np.abs(totalmass - 1.0) >0.0001:
        aprint("\ttotal mass is not equal to 1!!\n\t\t\tformula: %s\n\t\t\ttotalmas: %.10f. Rescaling."%(formula,totalmass),0)
        new_g4mat = {}
        new_totalmass = 0.
        factor = 1.0/totalmass
        for k in first:
            material    = k[0]
            massfrac    = float(k[1]) / 1000. * factor
            new_g4mat.update({material:{"mass fraction":massfrac, "material":material}})
            new_totalmass   += massfrac        
        if np.abs(new_totalmass - 1.0) >0.0001:
            aprint("\ttotal mass is AGAIN not equal to 1!! formula: %s totalmas: %.10f. Serious problem encuntered."%(formula,new_totalmass),0)
        else:
            aprint("\ttotal mass is now fine. Problem solved.",2)
        return new_g4mat
    return g4mat


def gdml_material(mat, name, temperature, density, pressure,state="solid"):
    """
    # units for T, rho, P are K, g/ccm, Pa
    This function is specially tailored for the solid materials, which are
    to be used for the core and the crust, but it can be also used to make
    gaseous materials. 
    """
    header = '\t<material name="%s" state="%s">\n'%(name,state)
    footer = '\t</material>\n'
    output = ""
    output += header
    output += '\t\t<D value="%e" unit="g/cm3" />\n'%density
    output += '\t\t<T value="%e" unit="K" />\n'%temperature
    output += '\t\t<P value="%e" unit="Pa" />\n'%pressure
    for key in mat:
        frac, name = mat[key]["mass fraction"],mat[key]["material"]
        if frac != 0.:
            current = '\t\t<fraction n="%e" ref="G4_%s" />\n'%(frac,name)
            output += current
    output += footer
    return output



def parse_header(filename):
    thefile = open(filename,'r')
    header = thefile.readline()
    thefile.close()
    header     = header.split()
    nrformulas = len(header)
    
    #Generate the header gdml text
    coremat  = parse_solid(header[0])
    
    #generate template dictionary for other sheets
    template   = {}
    idx = 0
    for material in header[1:]:
        template.update({material:{"mass fraction":0, "material":material, "colidx":idx}})
        idx += 1
    return coremat, template

def line_processor(filename):
    core,template = parse_header(filename)
    columnnames   = {"rmin":0,"rmax":1,"phimin":2,"phimax":3,"thetamin":4,
                     "thetamax":5,"SDID":6,"temperature":7,"density":8,"pressure":9,"core":10}
    #print "template:", template
    for mat in template:
        idx  = template[mat]["colidx"] + 11
        name = template[mat]["material"]
        columnnames.update({name:idx})
    # now we generate the gdml material blocks
    data          = np.loadtxt(filename)
    nrsheets      = data.shape[0]
    assert data.shape[1] == len(columnnames), "data shape and columnnames missmatch: (%ix%i) vs %i"%(data.shape[0],data.shape[1],len(columnnames))
    
    #We generate the materials and geometry entries
    materials, geometry, struct = [], [], []
    for aLine in data:
        custom = aLine[columnnames["core"]] == 1.0
        sdid = int(aLine[columnnames["SDID"]])
        name = "material" + str(sdid)
        temp = aLine[columnnames["temperature"]]
        dens = aLine[columnnames["density"]]
        pres = aLine[columnnames["pressure"]]
        if custom:
            gdmlmater = gdml_material(core, name, temp, dens, pres,"solid")
            materials.append(gdmlmater)
        else: 
            atmomat = template.copy()
            for mat in template:
                nname = template[mat]["material"]
                idx  = columnnames[mat]
                atmomat[nname]["mass fraction"] = aLine[idx]
            gdmlmater = gdml_material(atmomat, name, temp, dens, pres,"gas")
            materials.append(gdmlmater)
        #now the geometry  entry
        solname  = "solid" + str(int(aLine[columnnames["SDID"]]))
        rmin     = aLine[columnnames["rmin"]]
        rmax     = aLine[columnnames["rmax"]]
        phimin   = aLine[columnnames["phimin"]]
        phimax   = aLine[columnnames["phimax"]]
        phimax   = phimax - phimin
        if np.abs(phimax -2.0*np.pi) < 0.01: phimax = "TWOPI"
        thetamin = aLine[columnnames["thetamin"]]
        thetamax = aLine[columnnames["thetamax"]]
        thetamax = thetamax - thetamin
        if np.abs(thetamax -np.pi) < 0.01: thetamax = "PI"
        if thetamax == 2.0: thetamax = "PI"
        sphere   = '\t<sphere name="%s" rmin="%f" rmax="%f" startphi="%.4f" deltaphi="%s" starttheta="%.4f" deltatheta="%s" lunit="km" aunit="rad"/>\n'%(solname,rmin,rmax,phimin,str(phimax),thetamin,str(thetamax))
        geometry.append(sphere)
        struct.append(['volume'+str(sdid),name,solname,sdid])
    wrld = data[-1,1]
    return materials, geometry, struct, wrld

def gdml_compile(materials, world, geometry, structure, filename):
    finalfile = ""
    header = '<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">\n\n\n\n'
    footer = '\n<setup name="Exoplanet" version="1.0">\n\t<world ref="World"/>\n</setup>\n</gdml>'
    finalfile += header
    
    #now the <define> part
    #here we need to modify the offset of the world box. 
    define = '<define>\n\t<constant name="PI" value="1.*pi"/>\n\t<constant name="TWOPI" value="2.*pi"/>\n\t<position name="center"/>\n\t<rotation name="identity"/>\n</define>\n'
    finalfile += define
    
    #Now the materials
    finalfile += "<materials>\n"
    #We also add the AtRIS material extensions
    f = open("atris.materials","r")
    ff = f.readlines()
    for fff in ff:
        finalfile += fff
    f.close()
    for mat in materials:
        finalfile += mat
    finalfile += "</materials>\n"
    #Now the geometry
    finalfile += "<solids>\n"
    #add the world
    wd         = world * 2
    add        = '\t<box name="WorldBox" x="%f" y="%f" z="%f" lunit="km"/>\n'%(wd,wd,wd)
    finalfile += add
    
    for vol in geometry:
        finalfile += vol
    finalfile += "</solids>\n"
    
    #now the structure
    finalfile += "<structure>\n"
    for st in structure:
        add = '\t<volume name="%s">\n\t\t<materialref ref="%s"/>\n\t\t<solidref ref="%s"/>\n\t\t<auxiliary auxtype="SensDet" auxvalue="%i"/>\n\t</volume>\n'%tuple(st)
        finalfile += add
    #now we need to add the world and place all the other volumes
    worldpart = '\t<volume name="World">\n\t\t<materialref ref="Vacuum"/>\n\t\t<solidref ref="WorldBox"/>\n'
    finalfile += worldpart
    for st in structure:
        vv = st[0]
        add = '\t\t<physvol>\n\t\t\t<volumeref ref="%s"/>\n\t\t\t<positionref ref="center"/>\n\t\t</physvol>\n'%vv
        finalfile += add
    finalfile += "\t</volume>\n"
    finalfile += "</structure>\n"
    #Now the footer
    finalfile += footer
    f = open(filename,"w")
    f.write(finalfile)
    f.close()

if __name__ == "__main__":
    os.system("clear")
    aprint("\t################################################################",1)
    aprint("\t$$                  PSF 2 GDML                                $$",1)
    aprint("\t################################################################",1)
    aprint("\tThis is a script which converts a psf file to a gdml file.",2)
    aprint("\tFor more information check the documentation or visit the wiki.",2)
    aprint("\tUse:",2)
    aprint("\t",2)
    aprint("\t\tpython psf2gdml.py filenameprefix",2)
    aprint("\t",2)
    aprint("\twhere filenameprefix is the filename without the extension .gdml",2)
    aprint("\t################################################################",1)
    aprint("\t################################################################",1)
    # load data, count the abundance and save the data:
    arg1 = sys.argv[1]
    if os.path.isfile(arg1+'.psf'):
        mat, geom, struct,wrld = line_processor(arg1+'.psf')
        wrld += 100
        gdml_compile(mat,wrld,geom,struct,arg1+'.gdml')
    else:
        aprint('\tThe file ' + arg1 + '.psf does not exist.',0)

