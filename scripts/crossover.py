import sys
import csv
import collections
import math
import numpy as np
import pandas as pd

def delete_false_co(codict,mdict):
    '''Remove false positive crossover events
    based on nTypedBetween and return filtered
    dict.
    '''
    print("Crossovers before filtering: " + str(len(codict)))
    # Use set as there can be quadruple crossovers
    # This leads to duplicate bad keys
    bad_keys = set()
    for key, value in codict.items():
        #nTypedBetween 
        ntb = value[9]
        ind = value[1]
        genpos = value[2]
        # Skip last CO for individual
        # These can't be used for filtering
        if ntb != "NA":
            chr_pos1 = value[0] + "_" + value[3]
            co_index1 = codict.keys().index(key)
            # Get following crossover entry
            co_index2 = co_index1 + 1
            co_key2 = list(codict)[co_index2]
            chr_pos2 = codict[co_key2][0] + "_" + codict[co_key2][3]

            #Physical marker positions
            phypos1 = mdict[chr_pos1].split("_")[1]
            phypos2 = mdict[chr_pos2].split("_")[1]
            phychr1 = mdict[chr_pos1].split("_")[0]
            phychr2 = mdict[chr_pos2].split("_")[0]
            if phychr1 == phychr2:
                phydist = int(phypos2) - int(phypos1)
                # Use absolute numbers
                phydist = abs(phydist)
            else:
                # If adjacent marker is from different chromosome
                # Set arbitrary high physical distance
                phydist = 1000000
            # Genetic distance between breakpoints           
            gendist = float(codict[co_key2][2]) - float(genpos)
           # Recombinations separated by less than 5 markers are false
           # Recombinations <10kb are considered gene conversions
            if int(ntb) < 5 or phydist < 10000 or gendist < 0.02:
                co_index1 = codict.keys().index(key)
                # Get following crossover entry
                co_index2 = co_index1 + 1
                co_key2 = list(codict)[co_index2]
                
                bad_keys.add(key)
                if ind == codict[co_key2][1]:
                    bad_keys.add(co_key2)
                else:
                    print("Aborted. Expected paired false positives.")
                    break
    # delete all the bad entries
    for bkey in bad_keys:
        codict.pop(bkey, None) 
    print("Crossovers after filtering: " + str(len(codict)))
    return codict

# Read marker positions and crossover into dicts
# These tables are generated with r/qtl locateXO and a custom script

# chr,IND,location,left,right,ileft,iright,gleft,gright,nTypedBetween
crossovers = sys.argv[1] 
# chr,phy_chr_pos,genpos
markers = sys.argv[2]

with open(crossovers, mode='r') as co:
    reader1 = csv.reader(co)
    next(reader1, None)
    # Use chr_ind_pos as unique key and all row as values
    codict = collections.OrderedDict((row[0] + "_" + row[1] + "_" + row[2], row[0:]) for row in reader1)
with open(markers, mode='r') as mark:
    reader2 = csv.reader(mark)
    next(reader2, None)
    mdict = collections.OrderedDict((row[0] + "_" + row[2], row[1]) for row in reader2)

# Filter false positive CO
filtdict = delete_false_co(codict,mdict)
# Physical distances between CO
ind_list = []
co = {}
for key, value in filtdict.items():
    ind = value[1]
    chrom = value[0]
    lchrpos = value[0] + "_" + value[3]
    lphypos = mdict[lchrpos].split("_")[1]
    lphychr = mdict[lchrpos].split("_")[0]
    rchrpos = value[0] + "_" + value[4]
    rphypos = mdict[rchrpos].split("_")[1]
    rphychr = mdict[rchrpos].split("_")[0]
    
    if rphychr == chrom and lphychr == chrom:
        breakpoint = int((int(lphypos) + int(rphypos)) / 2)
        if chrom in co:
            # append CO to list
            breaklist = co[chrom]
            breaklist.append(breakpoint)
            co[chrom] = breaklist
        else:
            co[chrom] = [breakpoint]

    ind_list.append(ind)

# CO counts per individual
counter=collections.Counter(ind_list)
ind_list = []
frq = []
for ind in counter.most_common():
    ind_list.append(ind[0])
    frq.append(ind[1])
inddf = pd.DataFrame({'Ind' : ind_list, 'CO' : frq})
inddf.sort_values(['Ind'], ascending=True, axis=0,inplace=True)
inddf.to_csv('Individual_CO_Frequency.csv', index=False, encoding='utf-8', columns=['Ind','CO'])

# Histogram of CO across chromsomes
chromlist = []
histlist = []
binlist = []

for key, value in co.items():
    chrom = key
    breakpoints = np.array(value)
    # Round highest position up to nearest 500kb
    toprange = math.ceil(breakpoints.max() / 500000.0) * 500000.0
    # Number of bins
    binnum = int(toprange / 500000.0)
    # Calculate histogram
    hist,bins = np.histogram(breakpoints,bins=binnum,range=(0, toprange))
    histlist.extend(np.ndarray.tolist(hist))
    bins_tidy = np.ndarray.tolist(bins)[:-1]
    binlist.extend(bins_tidy)
    for bin in bins_tidy:
        chromlist.append(chrom)
codf = pd.DataFrame({'Chr' : chromlist, 'Bin' : binlist, 'Histogram' : histlist})
codf.sort_values(['Chr', 'Bin'], ascending=True, axis=0,inplace=True)
codf.to_csv('CO_hist.csv', index=False, encoding='utf-8', columns=['Chr','Bin','Histogram'])
