# This is the **alternate version** of the first sub-script within the greater LASSI pipeline, though I don't expect that as the user you'll need to do anything with this, since the iterator handles everything
# "Alternate" (above) simply means that we're not generating new un-truncated frequency spectra...using existing spectra, we're generating new truncations (for example, if you want to use existing data to make a K15 spectrum after already having made a K20)

# This script was written by Alexandre M. Harris (contact: amh522@psu.edu OR xistheway@gmail.com), please cite appropriately if you use this script
# LASSI is free to use and modify (once again, with the proper attribution)
###################################################################################################################

import sys
import numpy as np

Ucounts = open(sys.argv[1], "r+") # File containing unmodified counts
Ufreqs = open(sys.argv[2], "r+") # File containing unmodified spectra
popselect = sys.argv[3] # Repeated args are the same as for LASSI_spectrum_and_Kspectrum.py
kval = int(sys.argv[4])
neutyn = sys.argv[5]
scanMAXID = sys.argv[6]

if "hap" in sys.argv[1]:
    theploidy = "hap"
elif "mlg" in sys.argv[1]:
    theploidy = "mlg"

def process_Kspectra(countU, spectU): #gets the count and frequency spectra (U and K) for an input; get sampsize from Ucounts file
    countU = map(float, countU)
    spectU = map(float, spectU)

    sampsize = sum(countU)

    Kcount = countU[:kval]
    Kspect = spectU[:kval] #even if kval > len(fullspect), it'll work

    Kcount = map(lambda x: (x / sum(Kcount)) * sampsize, Kcount)
    Kspect = map(lambda x: x / sum(Kspect), Kspect)

    while len(Kcount) < kval: #should help ensure that all Koutput lines are of the same length
        Kcount.append(0)
        Kspect.append(0)

    return(Kcount, Kspect)

def output_Kprotocol(innameC, innameS, ploidy, thepop, outdata, nyn):
    outname_Kcount = innameC[:innameC.index("_count")] + "_K" + str(kval) + "_count_K" + ploidy + "_" + thepop + ".txt"

    dest_Kcount = open(outname_Kcount, "a+")

    sys.stdout = dest_Kcount

    print(" ".join(str(entry) for entry in outdata[0]))

    dest_Kcount.close()


    outname_Kplo = innameS[:innameS.index("_spectrum")] + "_K" + str(kval) + "_spectrum_K" + ploidy + "_" + thepop + ".txt"

    dest_Kplo = open(outname_Kplo, "a+")

    sys.stdout = dest_Kplo

    print(" ".join(str(entry) for entry in outdata[1]))

    dest_Kplo.close()


    if nyn == "yes": #doing the exact same thing as above, but now it'll go into an allreps_ file as well (and join window data from the other replicates)
        outname_allreps = "allreps_spectra_" + scanMAXID + "_K" + str(kval) + ploidy + "_" + thepop + ".txt"

        dest_allreps = open(outname_allreps, "a+")

        sys.stdout = dest_allreps

        print(" ".join(str(entry) for entry in outdata[1]))

        dest_allreps.close()    

while True:
    currcount = Ucounts.readline()

    if (currcount != ''): #a new line has been read (we're not at the end of the file)
        countsplit = currcount.split()

        currfreqs = Ufreqs.readline()
        freqsplit = currfreqs.split()

        Kspectra = process_Kspectra(countsplit, freqsplit)

        output_Kprotocol(sys.argv[1], sys.argv[2], theploidy, popselect, Kspectra, neutyn)
    else:
        sys.stdout = sys.__stdout__
        message = "Truncation of data files \'" + sys.argv[1] + "\' and \'" + sys.argv[2] + "\' to " + str(kval) + " is complete!"
        print(message)
        break

Ucounts.close()
Ufreqs.close()