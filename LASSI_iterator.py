# Welcome to LASSI! This is a script written in Python 2.7 that identifies and classifies signatures of selective sweeps in a study population, from haplotype, OR multilocus genotype data (though it works a little better for haplotypes)
# Please consult the user manual included with this software for operation instructions including command line and input file formatting, as well as example data

# The iterator is the central script in the LASSI pipeline [_spectrum, average_, _scan] with which users will interact, but it must be used in conjunction with the other Python (.py) scripts provided in the download package

# This script was written by Alexandre M. Harris (contact: amh522@psu.edu OR xistheway@gmail.com), please cite appropriately if you use this script
# LASSI is free to use and modify (once again, with the proper attribution)
###################################################################################################################

import sys
import subprocess

option = sys.argv[1] # Choose an option from the "pipeline" function below ["initial", "rescan", "neutavg", or "MLcalc"]

chromnum = sys.argv[2] # What chromosome are we scanning? Give its designation, such as 1, 2R, X etc.

kval = sys.argv[3] # What is the truncated frequency spectrum we're using for inference?

scanMAXID = sys.argv[4] # Provide a descriptive ID for the dataset (for the user's reference)

chromfile_fl = sys.argv[5] # Comma-separated strings: what is the name of the file *preceding*,*after* the chromosome number? For example, if your file is drosophila_chr3L_103119_hap.txt, this argument is drosophila_chr,_103119_hap.txt [[MAKE SURE _hap OR _mlg is in the filename!!!]] [[NOTE THAT INFILE NEEDS TO BE .txt FOR THIS TO WORK]]

chromfile_fl = chromfile_fl.split(",")
inchrom = chromfile_fl[0] + chromnum + chromfile_fl[1] # This is the full name of your file being scanned

if "_hap" in inchrom:
    theploidy = "hap"
elif "_mlg" in inchrom:
    theploidy = "mlg"
else:
    print("Infile name malformatted!")
    exit()

popselect = sys.argv[6] # What population are we analyzing? Use the exact formatting from the header file

if ((option == "initial") or (option == "rescan")): 
    headfile = sys.argv[7] # Name the file that has the population header...the header is needed in order to filter results properly and only keep in what is desired

    winsize = sys.argv[8] # What is the number of SNPs desired in each window? **Needs to be an odd number so that the core SNP can be unambiguously defined!**

    winshift = sys.argv[9] # By how many SNPs does the window advance with each increment? Usually good to have it be 1/10 the winsize

    output_all = sys.argv[10] # Should the frequency spectrum for each window also go an allreps_ file so it can be averaged to get the neutral background spectrum? [THIS SHOULD DEFINITELY BE "yes" for "initial" AND PROBABLY "yes" FOR "rescan", unless you have a good reason not to do it]

elif option == "neutavg":
    allreps_lines = sys.argv[7] # How many lines are in the allreps_ file to be averaged? [User needs to use wc -l from the terminal to determine this]

elif option == "MLcalc":
    out_option = sys.argv[7] # What likelihood function are we using? (1-5, see manual)

else:
    print("Pick an actual option!")
###################################################################################################################

def get_spect(inf,hdf,spectopt,ID,wsz,wsh,psl,nyn,ploidy,fcs): # Infile, spectrum option (make a new one or re-truncate an existing one), scanMAXID, yes/no to indicate whether the Kspectrum is additionally sent to an "allreps" file for further averaging into a neutral background spectrum, hap/mlg to decide what we're computing
    if spectopt == "UK": # Start everything from scratch
        spectcom = "python LASSI_spectrum_and_Kspectrum.py " + inf + " " + hdf + " " + wsz + " " + wsh + " " + psl + " " + kval + " " + nyn + " " + ID + " " + ploidy + " " + chromnum

        subprocess.call(spectcom, shell = True)
    elif spectopt == "Konly": #read in an existing U spectrum and make a new K spectrum
        if ploidy == "hap":
            countsname = fcs + "_count_Uhap_" + psl + ".txt"
            hspectname = fcs + "_spectrum_Uhap_" + psl + ".txt"
        elif ploidy == "mlg":
            countsname = fcs + "_count_Umlg_" + psl + ".txt"
            hspectname = fcs + "_spectrum_Umlg_" + psl + ".txt"

        spectcom = "python LASSI_Kspectrum_ONLY.py " + countsname + " " + hspectname + " " + psl + " " + kval + " " + nyn + " " + ID

        subprocess.call(spectcom, shell = True)

def avg_neut(inall,outavg,tlines):
    avgcom = "python LASSI_average_Kspect.py " + inall + " " + outavg + " " + tlines + " " + kval

    subprocess.call(avgcom, shell = True)

def run_ML_scan(inmaster,inct,ploidy,oopt,emin): # Assumes k=m
    MLcom = "python LASSI_scan.py " + inmaster + " " + inct + " " + wincents + " " + chromnum + " " + emin + " " + kval + " " + scanMAXID + " " + ploidy + " " + oopt

    subprocess.call(MLcom, shell = True)
###################################################################################################################

epsilon_min = str(1 / (float(kval) * 100)) # We'll optimize epsilons no smaller than this value

repneut_out = "allreps_spectra_" + scanMAXID + "_K" + kval + theploidy + "_" + popselect + ".txt"
totalneut_out = "total_spectrum_" + scanMAXID + "_K" + kval + theploidy + "_" + popselect + ".txt"

if ((option == "initial") or (option == "rescan") or (option == "MLcalc")):
    fileconstruct = inchrom[:inchrom.index(".txt")] # Easy way to find the needed files (they just get specific suffixes)

    wincents = "window_centers_" + scanMAXID + "_" + popselect + "_chr" + chromnum + ".txt" # Name of the file containing window centers for the analysis (so we can pair up signals with coordinates)

    if option == "initial":
        get_spect(inchrom,headfile,"UK",scanMAXID,winsize,winshift,popselect,output_all,theploidy,fileconstruct)

    elif option == "rescan":
        get_spect(inchrom,headfile,"Konly",scanMAXID,winsize,winshift,popselect,output_all,theploidy,fileconstruct)

    elif option == "MLcalc":
        incounts = fileconstruct + "_K" + kval + "_count_K" + theploidy + "_" + popselect + ".txt"

        run_ML_scan(totalneut_out,incounts,theploidy,out_option,epsilon_min)

elif option == "neutavg":
    avg_neut(repneut_out,totalneut_out,allreps_lines)

print("All done with %s, chromosome %s!" % (option, chromnum))