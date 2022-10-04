# This is the first sub-script within the greater LASSI pipeline, though I don't expect that as the user you'll need to do anything with this, since the iterator handles everything

# This script was written by Alexandre M. Harris (contact: amh522@psu.edu OR xistheway@gmail.com), please cite appropriately if you use this script
# LASSI is free to use and modify (once again, with the proper attribution)
###################################################################################################################

import sys
import numpy as np
from collections import defaultdict
import gzip

np.set_printoptions(threshold=np.nan)

if ".gz" in sys.argv[1]:
    infile = gzip.open(sys.argv[1], "rb")
else:
    infile = open(sys.argv[1], "r+") # All arguments constructed in iterator
    
headfile = open(sys.argv[2], "r+")
winsize = int(sys.argv[3])
winshift = int(sys.argv[4])
popselect = sys.argv[5]
kval = int(sys.argv[6])
neutyn = sys.argv[7]
scanMAXID = sys.argv[8]
ploidy = sys.argv[9]
chromnum = sys.argv[10]

def dictfreqs(VALall, sampsize): # Input is a nested list of sequences
    countdict = defaultdict(int)

    for string in VALall:
        countdict[string] += 1

    datacounts = []
    datafreqs = []

    for item in countdict.values():
        datacounts.append(float(item))

        freq = float(item) / sampsize
        datafreqs.append(freq)

    datacounts = sorted(datacounts, reverse=True)
    datafreqs = sorted(datafreqs, reverse=True)

    datacounts = map(str, datacounts)
    datafreqs = map(str, datafreqs)

    datacounts = " ".join(datacounts)
    datafreqs = " ".join(datafreqs)

    return(datacounts, datafreqs)

def get_empir_freqs(matr, sampsize, pldy): #input is transposed 1000G
    strings = []

    for ind in range(matr.shape[0]): #iterate through number of rows, which after transpose is number of individuals

        if pldy == "hap":
            current_string_h1 = 0
            current_string_h2 = 0

            for item in matr[ind]:
                if current_string_h1 == 0: #first loop to build string
                    if item == "1":
                        current_string_h1 = "0"
                        current_string_h2 = "0"
                    elif item == "2":
                        current_string_h1 = "0"
                        current_string_h2 = "1"
                    elif item == "3":
                        current_string_h1 = "1"
                        current_string_h2 = "0"
                    elif item == "4":
                        current_string_h1 = "1"
                        current_string_h2 = "1"

                else: #concatenate current SNP
                    if item == "1":
                        current_string_h1 = current_string_h1 + "0"
                        current_string_h2 = current_string_h2 + "0"
                    elif item == "2":
                        current_string_h1 = current_string_h1 + "0"
                        current_string_h2 = current_string_h2 + "1"
                    elif item == "3":
                        current_string_h1 = current_string_h1 + "1"
                        current_string_h2 = current_string_h2 + "0"
                    elif item == "4":
                        current_string_h1 = current_string_h1 + "1"
                        current_string_h2 = current_string_h2 + "1"
                    elif item == "N":
                        current_string_h1 = current_string_h1 + item
                        current_string_h2 = current_string_h2 + item

                    if len(current_string_h2) == len(matr[ind]): #compare the number of characters (SNPs) in string sequence to length of row in array to indicate completion of haplotype
                        strings.append(current_string_h1)
                        strings.append(current_string_h2)

                # return(dictfreqs(strings, (sampsize * 2))) #CHECK TO SEE IF THIS IS CORRECT!

        elif pldy == "mlg":
            current_string_mlg = 0

            for item in matr[ind]:
                if current_string_mlg == 0: #first loop
                    current_string_mlg = item

                else: #concatenate current SNP
                    current_string_mlg = current_string_mlg + item
        

                if len(current_string_mlg) == len(matr[ind]):
                        strings.append(current_string_mlg)

    return(dictfreqs(strings, sampsize))


def update_SNPs(matchsites, SNPcurr, winbuild, splitpop, windata): #adds the SNP info to the list as we're building the analysis window [[NEEDS TO KEEP TRACK OF THE WINDOW CENTERS!]]
    popsnps = [] #list of allelic states in the study sample (extracted from the entire data at that global SNP)

    for ind in matchsites: #which indices are valid for chosen populations?
        popsnps.append(splitpop[ind + 5]) #splitpop[ind + 5] is the index for the allelic state at that SNP for a properly-formatted file

    popsnps_set = set(popsnps) #number of site types [need to make sure that the SNP is a SNP in the population]

    if ((len(popsnps_set) > 1) or (("2" or "3" or "5") in popsnps)):
        winbuild.append(SNPcurr)
        windata.append(popsnps) #for this SNP, report the allelic states; ultimately nested list of haplo/mlg

def process_spectra(spectrum, sampsize): # Gets the count and frequency spectra (U and K) for an input; sampsize should be length of match
    fullcount = spectrum[0].split()
    fullspect = spectrum[1].split()

    fullcount = map(lambda x: int(float(x)), fullcount)

    Kcount = map(float, fullcount)
    Kspect = map(float, fullspect)

    Kcount = Kcount[:kval]
    Kspect = Kspect[:kval]

    Kcount = map(lambda x: (x / sum(Kcount)) * sampsize, Kcount)
    Kspect = map(lambda x: x / sum(Kspect), Kspect)

    while len(Kcount) < kval: # Fills in zeroes if the window doesn't have "enough" SNPs
        Kcount.append(0)
        Kspect.append(0)

    return(fullcount, Kcount, fullspect, Kspect)

def output_protocol(inname, pldy, thepop, outdata, nyn):
    outname_Ucount = inname[:inname.index(".t")] + "_count_U" + pldy + "_" + thepop + ".txt"

    dest_Ucount = open(outname_Ucount, "a+")

    sys.stdout = dest_Ucount

    print(" ".join(str(entry) for entry in outdata[0]))

    dest_Ucount.close()


    outname_Kcount = inname[:inname.index(".t")] + "_K" + str(kval) + "_count_K" + pldy + "_" + thepop + ".txt"

    dest_Kcount = open(outname_Kcount, "a+")

    sys.stdout = dest_Kcount

    print(" ".join(str(entry) for entry in outdata[1]))

    dest_Kcount.close()


    outname_Uplo = inname[:inname.index(".t")] + "_spectrum_U" + pldy + "_" + thepop + ".txt"

    dest_Uplo = open(outname_Uplo, "a+")

    sys.stdout = dest_Uplo

    print(" ".join(str(entry) for entry in outdata[2]))

    dest_Uplo.close()


    outname_Kplo = inname[:inname.index(".t")] + "_K" + str(kval) + "_spectrum_K" + pldy + "_" + thepop + ".txt"

    dest_Kplo = open(outname_Kplo, "a+")

    sys.stdout = dest_Kplo

    print(" ".join(str(entry) for entry in outdata[3]))

    dest_Kplo.close()

    if nyn == "yes": #doing the exact same thing as above, but now it'll go into an allreps_ file as well (and join window data from the other replicates)
        outname_allreps = "allreps_spectra_" + scanMAXID + "_K" + str(kval) + pldy + "_" + thepop + ".txt"

        dest_allreps = open(outname_allreps, "a+")

        sys.stdout = dest_allreps

        print(" ".join(str(entry) for entry in outdata[3]))

        dest_allreps.close()

populations = (headfile.readline()).split()
headfile.close()

poparray = np.array(populations)

match = np.flatnonzero(poparray == popselect)
match = match.tolist()

if ploidy == "hap":
    thessz = len(match) * 2
elif ploidy == "mlg":
    thessz = len(match)

window_build = [] #keeps track of the window size...once it's the correct size, output results and make new window

window_data = []

outname_wincents = "window_centers_" + scanMAXID + "_" + popselect + "_chr" + chromnum + ".txt"

while True:
    position = infile.readline()

    if (position != ''): # A new line has been read (we're not at the end of the file)
        posplit = position.split()
        current_SNP = int(posplit[2]) # This indexing works for properly-formatted infiles

        if len(window_build) < winsize:
            update_SNPs(match, current_SNP, window_build, posplit, window_data)
        else:
            t_window_data = np.transpose(np.array(window_data))

            returnspects = get_empir_freqs(t_window_data, thessz, ploidy)

            thespectra = process_spectra(returnspects, thessz) #returns a tuple of length 4

            output_protocol(sys.argv[1], ploidy, popselect, thespectra, neutyn)

            winmedian = np.median(window_build) #what is the window center?

            dest_wincents = open(outname_wincents, "a+")
            sys.stdout = dest_wincents
            print(winmedian)
            dest_wincents.close()

            window_build = window_build[winshift:] #get rid of the SNPs that aren't in the new window
            window_data = window_data[winshift:]

            update_SNPs(match, current_SNP, window_build, posplit, window_data)
    else:
        t_window_data = np.transpose(np.array(window_data))

        returnspects = get_empir_freqs(t_window_data, thessz, ploidy)

        thespectra = process_spectra(returnspects, thessz)

        output_protocol(sys.argv[1], ploidy, popselect, thespectra, neutyn)

        winmedian = np.median(window_build)

        dest_wincents = open(outname_wincents, "a+")
        sys.stdout = dest_wincents
        print(winmedian)
        dest_wincents.close()

        sys.stdout = sys.__stdout__
        message = "Analysis of data file \'" + sys.argv[1] + "\' is complete! Final window contains " + str(len(window_build)) + " SNPs."
        print(message)
        break