# Final script in the LASSI pipeline, whose purpose is to take the neutral background spectrum (total_neutral_) and the counts spectra for each window (_count_K) and compute the T-statistic for each window
# Output includes T-statistic, inferred m, epsilon, p_K, likelihood option, and genomic coordinate (see manual)
# Using the output from this scanning script, you can then intersect the coordinates of genes/regulatory elements/etc and see if any peaks occur at interesting locations

# This script was written by Alexandre M. Harris (contact: amh522@psu.edu OR xistheway@gmail.com), please cite appropriately if you use this script
# LASSI is free to use and modify (once again, with the proper attribution)

#(071019) Changed output protocol of scan results--new results output with each line instead of all at once at the end; prevents the thread from timing out, but it's probably not as pythonic
###################################################################################################################

import sys
import numpy as np
import math

allspectra = sys.argv[1] # Neutral background spectrum for the dataset (total_neutral_)

samplecounts = open(sys.argv[2], "r+") # Counts spectra across all windows of the specified chromosome

centersfile = open(sys.argv[3], "r+") # Window centers for each window, so that the T-statistic output can be matched to the genomic coordinate

chromnum = sys.argv[4] # Set by iterator

epsilon_min = float(sys.argv[5]) # Should be set by the iterator to be 1/(100K); there's no hard rule but that's usually a good starting point

m_vals = int(sys.argv[6]) + 1 # Upper bound of m for the analysis; adds a +1 because of Python's indexing

scanMAXID = sys.argv[7] # Set by iterator

ploidy = sys.argv[8] # Set by iterator

out_option = sys.argv[9] # Set by iterator; which likelihood function are we using?
###################################################################################################################

def easy_lik(thespect, counts, kval): # Basic computation of the likelihood function; runs as-is for neutrality, but called as part of a larger process for sweep model
    liklist = []

    for i in range(kval):
        lik_part = counts[i] * math.log(thespect[i])
        liklist.append(lik_part)

    liklist = sum(liklist)

    return(liklist)

def alt_lik_ALL(nullspect, counts, kval, mval, epsil, oopt): # Computes the likelihood of a sweep under optimized parameters
    epsmax = nullspect[len(nullspect) - 1]
    
    if mval != kval:
        altspect = []

        tailclasses = []
        neutdiff = []
        tailinds = range(mval+1, kval+1)

        for ti in tailinds:
            try:
                the_ns = epsmax - (float(ti - mval - 1)/float(kval - mval - 1))*(epsmax - epsil)
                tailclasses.append(the_ns)
                neutdiff.append(nullspect[ti - 1] - the_ns)
            except:
                if float(ti - mval - 1) == float(kval - mval - 1):
                    tailclasses.append(epsil)
                    neutdiff.append(nullspect[ti - 1] - epsil)

        headinds = range(1, mval+1)

        for hd in headinds:
            altspect.append(nullspect[hd - 1])

        neutdiff_all = sum(neutdiff)

        for ival in headinds:
            if oopt == 1:
                theadd = ((1 / float(ival)) / sum(map(lambda x: 1 / float(x), headinds))) * neutdiff_all
            elif oopt == 2:
                theadd = ((1 / float(ival ** 2)) / sum(map(lambda x: 1 / float(x ** 2), headinds))) * neutdiff_all
            elif oopt == 3:
                theadd = (math.exp(float(-ival)) / sum(map(lambda x: math.exp(float(-x)), headinds))) * neutdiff_all
            elif oopt == 4:
                theadd = (math.exp(float(-ival**2)) / sum(map(lambda x: math.exp(float(-x**2)), headinds))) * neutdiff_all
            elif oopt == 5:
                theadd = (1 / float(mval)) * neutdiff_all

            altspect[ival - 1] += theadd

        for tc in tailclasses:
            altspect.append(tc)

    else:
        altspect = nullspect

    return(easy_lik(altspect, counts, kval))

def indata_to_array(indata): #make sure that it can be coerced to a float!
    with open(indata, "r+") as indd:
        idt = indd.readlines()

    idt = map(lambda x: x.split(), idt)

    idt = np.array(idt).astype("float")

    return(idt)

asp = indata_to_array(allspectra)

epsvals = map(lambda x: x * epsilon_min, range(1, 101))

thenull = asp.tolist()[0]

the_kval = len(thenull)

outname = scanMAXID + "_K" + str(the_kval) + "_chr" + chromnum + ploidy + ".txt"

while True:
    wincent = centersfile.readline()

    try:
        wincent = float(wincent)

        scts = samplecounts.readline().split()
        scts = np.array(scts).astype(float)
        scts = scts.tolist()

        null_lik = easy_lik(thenull, scts, the_kval)

        emax = thenull[len(thenull) - 1]
        epsvals = [ev for ev in epsvals if ev <= emax]

        alt_likelihoods_byE = []

        for vale in epsvals:
            alt_likelihoods_byM = []

            for valm in range(1, m_vals):
                altlik = alt_lik_ALL(thenull, scts, the_kval, valm, vale, int(out_option))
                alt_likelihoods_byM.append(altlik)

            likelihood_bestM = 2 * (max(alt_likelihoods_byM) - null_lik)

            if likelihood_bestM > 0:
                MLmaxM = (alt_likelihoods_byM.index(max(alt_likelihoods_byM))) + 1
            else:
                MLmaxM = 0

            alt_likelihoods_byE.append([likelihood_bestM, MLmaxM, vale])

        alt_likelihoods_byE = np.array(alt_likelihoods_byE)

        likelihood_real = max(alt_likelihoods_byE[:,0])

        out_index = np.flatnonzero(alt_likelihoods_byE[:,0] == likelihood_real)

        out_intermediate = alt_likelihoods_byE[out_index]

        if out_intermediate.shape[0] > 1:
            constarg = min(out_intermediate[:,1])

            outcons = np.flatnonzero(out_intermediate[:,1] == constarg)

            out_cons_intermediate = out_intermediate[outcons]

            if out_cons_intermediate.shape[0] > 1:
                out_cons_intermediate = out_cons_intermediate[0]

            out_intermediate = out_cons_intermediate

        outshape = out_intermediate.shape

        if len(outshape) != 1:
            out_intermediate = out_intermediate.tolist()[0]
        else:
            out_intermediate = out_intermediate.tolist()

        out_intermediate.append(thenull[len(thenull) - 1]) # Same for each window, but may be desired for reference
        out_intermediate.append(out_option)
        out_intermediate.append(wincent)

        dest = open(outname, "a+")
        sys.stdout = dest
        print(" ".join(str(entry) for entry in out_intermediate))
        dest.close()
    except:
        sys.stdout = sys.__stdout__
        message = "Scan of of data file \'" + sys.argv[2] + "\' using neutral spectrum \'" + allspectra + "\' and window centers \'" + sys.argv[3] + "\' is complete!"
        print(message)
        break

samplecounts.close()
centersfile.close()