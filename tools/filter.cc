#include <vector>
#include <string>
#include <sys/stat.h>
#include <iostream>
#include <fstream>
using namespace std;

/* Author: S. Banjac
 * Leibnitzstr. 11
 * LS11-R308 
 * 0431 880 1191
 * banjac@physik.uni-kiel.de
 */

/* Description:
 *====================================================================================
 * The AtRIS simulation toolkit creates a lot of binary data, especially in mode 0. 
 *
 *      [PDG, next SDID, angle, ekin, esource]
 *
 *    * PDG specifies the particle species
 *        "(signed) int" ::  4 bytes
 *    * next SDID specifies the Sensitive Detector ID of the 
 *      detector in to which the particle is entering. This is:
 *        "usigned short" :: 2 bytes
 *    * Angle is the angle between the position vector and
 *      the direction vector. To save storage, this is an:
 *        "unsigned char"::  2 bytes
 *    * Kinetic energy of the particle at the interface. Type:
 *       "float" ::          4 bytes
 *    * Kinetic energy of the primary particle. Type:
 *       "float" ::          4 bytes
 *
 * Usage:
 *      filter simname PDG IID
 *
 * Examples:
 *    * simname    :: adleo
 *    * pdg        :: 2212    (neutrons)
 *    * iid        :: 0       (all interfaces)
 *
 * In this program, we want to implement:
 *    * find out which integer numbers can be found in the
 *      binary data names within the current folder. Remember
 *      AtRIS batch scripts pick random integers to initialize
 *      the PRNG and then use the first integer to generate the
 *      names which are used for the output.
 *    * Here is how processing works:
 *          i. create String sufix = simname + ".binary"
 *         ii. check in existing folder for all files 
 *             ending with "sufix"
 *        iii. open each file
 *         iv. process each line
 *          v. if pdg and iid match, write data out.
 *         vi. use sufix (without initial integer seeds)
 *             as output.
 *    * calculate from SDID and angle the (interface) IID
 *    * filterint of particles with a given PDG
 *    -> we set the IID argument to 666
 *    * filtering of particles with a given IID
 *    -> we set the PDG argument to 666
 *    * filtering of particles with a given IID and PDG
 *    * joining of all binary files in to a single file:
 *    -> we set both the PDG and the IID arguments to 666
 *
 * Compilation:
 *   * we need to use the -std=c++11 flag for g++ since we are 
 *     using stoi. Compile with:
 *     ./g++ -std=c++11 filter.cc
 *
 *====================================================================================
 */

/* We define a structure Line and declare a function which is going to process it */
typedef unsigned short ushort;
typedef unsigned char uchar;
typedef unsigned int uint;

int main(int argc, char *argv[])
{
    int intsize = sizeof(int);
    int floatsize = sizeof(float);
    int ucharsize = sizeof(uchar);
    int ushortsize = sizeof(ushort);
    int len = intsize + ushortsize +ucharsize + floatsize + floatsize;

    if (argc != 4)
    {
        /* If no argument is passed to main, remaind user how to use the programm*/
        printf( "Usage: %s simname PDG IID\n", argv[0]);
	printf( "set PDG or IID to 666 if you want to ignore them\n");
	printf( "set both PDG or IID to 666 if you only wish to join the binaries\n");
	return 1;
    }
    string  sim  = argv[1];   
    string sufix = sim + ".binary";
    int    pdgc  = stoi(argv[2]);
    int    iidc  = stoi(argv[3]);
    string output= sim + "." + to_string(pdgc)+"_"+to_string(iidc);
    cout << "#################################################################" << endl;
    cout << "#################################################################" << endl;
    cout << "### AtRIS filter script #########################################" << endl;
    cout << "#################################################################" << endl;
    cout << "#################################################################" << endl;
    cout << "### AtRIS: simname \t\t::\t\t " << sim   << endl;
    cout << "### AtRIS: sufix   \t\t::\t\t " << sufix << endl;
    cout << "### AtRIS: output  \t\t::\t\t " << output<< endl;
    cout << "### AtRIS: PDG code\t\t::\t\t " << pdgc  << endl;
    cout << "### AtRIS: IID code\t\t::\t\t " << iidc  << endl;

    // reading out the integers used to initialize the PRNG and name the files
    cout << "### AtRIS: now parsing 'randoms' file to read seed ints..." << endl;
    string randoms = "randoms";
    int    files   = 0;
    int    curran;
    vector<int> rands;
    vector<string> filenames;
    ifstream is(randoms);
    while (is >> curran) {
	    rands.push_back(curran);
	    string newname = to_string(curran) + sufix;
	    filenames.push_back(newname);
	    files++;
    }
    cout << "### AtRIS: must process " << files << " files:" << endl;
    int existing = 0;
    for (int k=0; k<files; k++) {
	   // We check if the files exist
           struct stat buffer;
	   bool exists = (stat (filenames[k].c_str(), &buffer) == 0);
	   cout << "\t\t" << k << ". " << filenames[k] << " exists: " << exists << endl;
	   existing++;
    }
    if (existing == files){
	cout << "### AtRIS: all files exist, proceeding..." << endl;
	FILE * out;
        out = fopen(output.c_str(), "wb");
        if (out == 0) {
		printf("Could not open file for writing.\n");
		return 1;
	}
	// Now open each file after the other
	long int particles = 0;
	long int remaining = 0;
	FILE * fp;
	for (int k=0; k<files; k++) {
		cout << "### AtRIS: processing file " << filenames[k] << "... ";
		// temp variables to hold the data before saving
		int    one;
		ushort two;
		uchar three;
		float four;
		float five;
		int fr = 15;
		int frt = 0;
                fp = fopen(filenames[k].c_str(), "rb");
		cout << "starting file processing" << endl;
		if(fp == NULL) {cout << "fp is a NULL pointer!"; return 1;}
		while((fp != NULL) and (fr==15)) {
		   
		   // We  read the five values comprising a single data line
                   frt += fread(&one, intsize, 1, fp);
                   frt += fread(&two, ushortsize, 1, fp);
                   frt += fread(&three, ushortsize, 1, fp);
                   frt += fread(&four, floatsize, 1, fp);
                   frt += fread(&five, floatsize, 1, fp);
		   if (frt !=5) fr = frt;
		   particles++;
		   // we filter the particles
		   bool cond1 = ((pdgc==666) and (iidc==666));
		   bool cond2 = ((pdgc==one) and (iidc==666));
		   bool cond3 = ((pdgc==666) and (iidc==two));
		   bool cond4 = ((pdgc==one) and (iidc==two));
		   bool cond = cond1 or cond2 or cond3 or cond4;
		   //cout << cond1 << " " << cond2 << " " << cond3 << " " << cond4 << " " << cond << " fr: " << fr << " frt: " << frt << endl;
		   if (cond){
			fwrite(&one,intsize,1,out);
			fwrite(&two,ushortsize,1,out);
			fwrite(&three,ushortsize,1,out);
			fwrite(&four,floatsize,1,out);
			fwrite(&five,floatsize,1,out);
			remaining++;
		   }
		   frt = 0;
	    	}
		fclose(fp);
	}
    fclose(out);		
    cout << "done. remaining: " << remaining << " out of " << particles << endl;
    }
    return 0;
}
