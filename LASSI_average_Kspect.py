# This script is called after the first scan of the dataset; it is used to generate a background neutral frequency spectrum from which the likelihood function is computed (by the next/final script in the pipeline)
# The output will consist of a single line, a frequency spectrum

# This script was written by Alexandre M. Harris (contact: amh522@psu.edu OR xistheway@gmail.com), please cite appropriately if you use this script
# LASSI is free to use and modify (once again, with the proper attribution)
###################################################################################################################

import sys
import numpy as np

intable = sys.argv[1] # The allreps_ file produced by the first run of the pipeline
outfile = sys.argv[2] # Just say what you want it to be called
totallines = int(sys.argv[3]) # Use wc -l on your terminal to extract this value [this is something you need to do on your own]
kval = int(sys.argv[4]) # Recycled from previous scripts/commands

intb = open(intable, "r+")

dest = open(outfile, "a+")
sys.stdout = dest

splitlist = []
fulllist = []

theweights = []

count = 0

for line in intb:
    count += 1
    line = line.split()
    if len(line) == kval:
        splitlist.append(line)

    if ((count % 5e4 == 0) or (count == totallines)):
        splitlist = np.array(splitlist)

        splitlist = splitlist.astype(float)

        gwideK = np.mean(splitlist, axis=0)

        gwideK = gwideK.tolist()

        fulllist.append(gwideK)

        if count % 5e4 == 0:
	        theweights.append(5e4)
        elif count == totallines:
        	smallweight = count % 5e4

        	theweights.append(smallweight)

    	splitlist = []

finalout = np.average(fulllist, axis=0, weights=theweights)
finalout = finalout.tolist()
finalout = map(str, finalout)

print(" ".join(finalout))

dest.close()

intb.close()