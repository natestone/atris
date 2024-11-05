# This batcher iterates over a collection of physics lists.
# Energy is constant
from threading import Thread
import random
import numpy as np
import subprocess
import os
import sys


def genMacro(nr1,nr2,N,name):
    lines =  ""
    lines +=  "/run/initialize\n"
    lines += "/run/verbose 2\n"
    lines += "/gps/particle proton\n"
    lines += "/gps/ene/type User\n"
    lines += "/gps/hist/type energy\n"
    bins   = np.loadtxt(name+'.bins')
    lines += "/gps/hist/point %.9f %.9f\n"%(bins[0],0.0)
    for thebin in bins[1:]:
        lines += "/gps/hist/point %.9f %.9f\n"%(thebin,1.0)
    lines += "/gps/pos/type Surface\n"
    lines += "/gps/pos/shape Sphere\n"
    lines += "/gps/pos/centre 0 0 0 m\n"
    lines += "/gps/pos/radius 6434.5241 km\n" #THIS HAS TO BE CONFIGURED MANUALLY!
    lines += "/gps/ang/type cos\n"
    lines += "/gps/ang/mintheta 0 deg\n"
    lines += "/gps/ang/maxtheta 90 deg\n"
    lines += "/gps/number 1\n"
    lines += "/run/setCut 100. m\n"           #CONFIGURE
    lines += "/run/particle/applyCuts true e-\n"
    lines += "/run/particle/applyCuts true e+\n"
    lines += "/run/particle/applyCuts true gamma\n"
    lines += "/run/particle/dumpCutValues\n"
    lines += "/random/setSeeds %i %i\n"%(nr1,nr2)
    lines += "/run/beamOn %i\n"%N
    return lines

def manager(prefix,nr1,nr2,N):
    name = str(nr1) + prefix
    os.system('cp '+prefix+'.gdml ' + name + '.gdml')
    os.system('cp '+prefix+'.bins ' + name + '.bins')
    lines = genMacro(nr1,nr2,N,name)
    filehandle = open(name + '.mac', 'w')
    filehandle.writelines(lines)
    filehandle.close()
    command = './AtRIS QGSP_BERT_HP '+name+' 0 > ' + name + ".out"
    print command
    subprocess.Popen(command,shell=True)



if __name__ == "__main__":
    randoms = []
    name = sys.argv[1]
    number = int(sys.argv[2])
    thrd = int(sys.argv[3])
    for k in range(1,thrd+1):
        a = random.randint(31,999999)
        b = random.randint(31,999999)
        manager(name,a,b,number)
        randoms.append(a)

    o = open('randoms','a')
    np.savetxt(o,randoms,fmt='%i')
    os.system('watch -n 1 ./monitor.sh')

