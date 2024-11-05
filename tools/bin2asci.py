import numpy as np
import sys


if __name__ == "__main__":
    docstr  = "Welcome to bin2ascii, a script which converts your 1.binary and 2.binary data\n"
    docstr += "For data analysis we recomend working with the binary data. Please examine the\n"
    docstr += "code of this script to se how. Binary data reads in approx 30 times faster and\n"
    docstr += "consumes only 1/3 of the storage space\n"
    docstr += "use:' python binary2asci.py prefix', where prefix corresponds to the name part\n"
    docstr += "that comes before '.binary'\n"
    print docstr
    dt      = np.dtype("i4,u2,u2, (2)f4")
    prefix  = sys.argv[1]
    print "### AtRIS: reading binary data....",
    d = np.fromfile(prefix+".binary", dtype = dt)
    print ".........done."
    pdg = d['f0']
    det = d['f1']
    ang = d['f2']
    flt = d['f3']
    d = np.column_stack((pdg,det,ang,flt))
    print "### AtRIS: now saving the data....",
    fmt = "%i %i %i %f %f"
    np.savetxt(prefix+".asci", d,fmt=fmt)
    print ".........done. Thank you for using this script"
    print "additional debuging"
    print "col 1:", np.min(d[:,0]), np.max(d[:,0])
    print "col 1:", np.min(d[:,1]), np.max(d[:,1])
    print "col 1:", np.min(d[:,2]), np.max(d[:,2]), "downward ", np.sum(d[:,2]>90), " upward", np.sum(d[:,2]<90)
    print "col 1:", np.min(d[:,3]), np.max(d[:,3])
    print "col 1:", np.min(d[:,4]), np.max(d[:,4])

