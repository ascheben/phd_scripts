import sys
import os
import gzip
import io
import math
import statistics
import numpy as np
from pathlib import Path
from pybedtools import BedTool

# This script was written to parse mapping statistics
# output by samtools and mosdepth for multiple files
# in a specific directory structure.

def dict2csv(mydict):
    '''Convert a nested dictionary to
    a csv table. Input is a dictionary with
    a specific expected structure. Output is
    written a csv file out.csv.
    '''
    header = [
        '#Pop',
        'GBS',
        'Sample',
        'Mapper',
        'Reference',
        'Region',
        'Breadth',
        'Depth',
        'Bases',
        'CoveredBases',
        'Loci',
        'MeanLociPerRegion',
        'SDLociPerRegion',
        'Mapped',
        'ProperlyPaired',
        'MeanReadLength',
        'MeanInsertSize',
        'InsertSizeSD',
        'ReadsOnDiffChr'
    ]
    with open('out.csv','a') as outcsv:
        headerrow = ','.join(header)
        outcsv.write(headerrow + '\n')
        for key, value in mydict.items():
            print(key)
            sid = key.split('_')
            pop = sid[0]
            gbs = sid[1]
            sample = sid[2]
            mapper = sid[3]
            ref = sid[4]

            for innerkey,innervalue in value.items():
                print(innerkey,innervalue)
                if innerkey == 'samstat':
                    pass
                else:
                    # When a chrom has zero mapped reads
                    # the chrom will not appear in the mosdepth
                    # dist file and the value is thus set to 0.
                    try:
                        breadth = innervalue['Breadth']
                    except:
                        breadth = 0
                    bases = innervalue['Bases']
                    covbases = innervalue['Covered bases']
                    depth= innervalue['Mean depth covered bases']
                    loci = innervalue['Loci']

                    if innerkey == 'total':
                        totseq = value['samstat']['Total seq']
                        propair = value['samstat']['Properly paired']
                        meanlen = value['samstat']['Mean length']
                        meanins = value['samstat']['Mean insert size']
                        insd = value['samstat']['Insert size StDev']
                        diffchr = value['samstat']['On different chr']
                        region = 'total'
                        meanloci = innervalue['Mean loci per chr']
                        sdloci = innervalue['StDev loci per chr']
                    else:
                        region = innerkey
                        meanloci = innervalue['Mean loci per 100kb']
                        sdloci = innervalue['StDev loci per 100kb']
                        totseq = propair = meanlen = meanins = insd = diffchr = 'NA'

                    outstring = [
                        pop,
                        gbs,
                        sample,
                        mapper,
                        ref,
                        region,
                        breadth,
                        depth,
                        bases,
                        covbases,
                        loci,
                        meanloci,
                        sdloci,
                        totseq,
                        propair,
                        meanlen,
                        meanins,
                        insd,
                        diffchr
                    ]
                    outstring = ','.join(map(str, outstring))
                    #print("Writing to file: ", region)
                    outcsv.write(outstring + '\n')

def mergedict(a, b, path=None):
    '''Merges nested dictionary b into
    nested dictionary a using keys.'''
    if path is None: path = []
    for key in b:
        if key in a:
            if isinstance(a[key], dict) and isinstance(b[key], dict):
                mergedict(a[key], b[key], path + [str(key)])
            elif a[key] == b[key]:
                pass # same leaf value
            else:
                raise Exception('Conflict at %s' % '.'.join(path + [str(key)]))
        else:
            a[key] = b[key]
    return a

def file2id(filepath):
    '''Parses filepath and returns a list of ID
    components to generate a unique ID and determine
    which type of analysis to run on the file.
    Input path example:
    Population1/RAD/bam/bamstats/58_bwa_v81.q20.rg.stats
    Output:
    ['Pop1','RAD','58','bwa','v81','stats']
    '''
    my_file = Path(filepath)
    if my_file.is_file():
    # file exists
        head, tail = os.path.split(my_file)
        head = head.split('/')
        # Fixed file paths required
        # Population name
        pop = head[0]
        pop = pop.replace('Population','Pop')
        # GBS method
        gbs = head[1]
        sample = tail.split('_')[0]
        suffix = tail.split('.')[-1]
        # Mapping method
        if 'bwa' in tail:
            mapper = 'bwa'
        elif 'soap' in tail:
            mapper = 'soap'
        else:
            print("Error! Mapping method not found for " + filepath)
            mapper = 'NA'
        # Reference
        if 'v81' in tail:
            ref = 'v81'
        elif 'NRG' in tail:
            ref = 'NRG'
        else:
            print("Error! Reference not found for " + filepath)
            ref = 'NA'
        # File type
        if suffix == 'stats':
            filetype = 'stats'
        elif suffix == 'txt':
            filetype = 'mosdepth'
        elif suffix == 'gz':
            filetype = 'perbase'
        else:
            print("Error! File type not found for " + filepath)
            filetype = 'NA'

        idlist = [pop,gbs,sample,mapper,ref,filetype]
        return idlist

def string2num(s):
    '''
    Convert string to int or float.
    '''
    try:
        return int(s)
    except ValueError:
        return float(s)

def get_value(myline):
    '''
    Takes a string, strips it, splits it by colon,
    gets the second element, strips leading
    whitespaces from it and converts it to a float.
    '''
    myline = myline.strip()
    myline = myline.split(':')
    myline = myline[1]
    # Remove in-line comments
    if '#' in myline:
        myline = myline.split('#')[0]
    myline = myline.lstrip()
    mynum = string2num(myline)
    return mynum

def getmstats(instats):
    '''Retrieve seven summary mappings statistics
    of interest from samtools stats output. Input
    is a file path and output is a list of integers
    and floats.
    '''
    samdic = {}
    internal_dic = {}
    # Query prefix for Summary Numbers
    qp = "SN\t"
    q_totseq = qp + "raw total sequences"
    q_prpair = qp + "reads properly paired"
    q_meanlen = qp + "average length"
    q_maxlen = qp + "maximum length"
    q_insize = qp + "insert size average"
    q_insd = qp + "insert size standard deviation"
    q_dchr = qp + "pairs on different chromosomes"

    # Query only appears once in file
    # If query appeared twice objects would be overwritten
    for line in open(instats):
        if q_totseq in line:
            totseq = get_value(line)
        if q_prpair in line:
            prpair = get_value(line)
        if q_meanlen in line:
            meanlen = get_value(line)
        if q_maxlen in line:
            maxlen = get_value(line)
        if q_insize in line:
            insize = get_value(line)
        if q_insd in line:
            insd = get_value(line)
        if q_dchr in line:
            dchr = get_value(line)
    internal_dic['Total seq'] = totseq
    internal_dic['Properly paired'] = prpair
    internal_dic['Mean length'] = meanlen
    internal_dic['Mean insert size'] = insize
    internal_dic['Insert size StDev'] = insd
    internal_dic['On different chr'] = dchr
    samdic['samstat'] = internal_dic
    return samdic

def perbasestats(bedgz):
# takes a single .per-base.bed mosdepth output file
# calculates mean coverage of non-zero bases
    scaffdict = {}
    sum_covbases = 0
    sum_bases = 0
    chrset = set()
    # Use wrapper to get text and not binary
    f = io.TextIOWrapper(gzip.open(bedgz, 'r'))
    for l in f:
        l_arr=l.rstrip().split("\t")
        scaff = l_arr[0]
        chrset.add(scaff)
        start = int(l_arr[1])
        end = int(l_arr[2])
        cov = int(l_arr[3])
        if cov > 0:
            # Add chr to dict if not in dict
            bases = end - start
            covbases = bases * cov
            if scaff in scaffdict:
                sum_covbases = scaffdict[scaff][1] + covbases
                sum_bases = scaffdict[scaff][0] + bases
                mean_cov = sum_covbases / sum_bases
                scaffdict[scaff] = [sum_bases,sum_covbases,mean_cov]
            else:
                mean_cov = covbases / bases
                scaffdict[scaff] = [bases,covbases,mean_cov]
    totbases = 0
    totcovbases = 0
    for key, value in scaffdict.items():
        # Dict to nest in scaffdict
        statsdict = {}
        statsdict['Bases'] = value[0]
        statsdict['Covered bases'] = value[1]
        statsdict['Mean depth covered bases'] = round(value[2],4)
        scaffdict[key] = statsdict
        totbases = totbases + value[0]
        totcovbases = totcovbases + value[1]
    totmeancov = round((totcovbases / totbases),4)
    statsdict = {}
    statsdict['Bases'] = totbases
    statsdict['Covered bases'] = totcovbases
    statsdict['Mean depth covered bases'] = totmeancov
    scaffdict['total'] = statsdict
    f.close()
    # Add missing chr to dict
    # This ensures all chr are always there
    for chrom in chrset:
        if chrom not in scaffdict:
            emptydic = {}
            emptydic['Bases'] = 0
            emptydic['Covered bases'] = 0
            emptydic['Mean depth covered bases'] = 0
            scaffdict[chrom] = emptydic
    return scaffdict

def radloci(bedgz):
    '''
    Use pybedtools to get estimated
    callable RAD loci positions. Return
    number of loci.
    '''
    bedcov = BedTool(bedgz)
    filtcov = []
    locidic = {}
    chrsize = {}
    chrset = set()
    for line in bedcov:
        mychr = line[0]
        chrset.add(mychr)
        # Locus must have minimum coverage of 4
        if int(line[3]) >3:
            bedrow = (line[0],int(line[1]),int(line[2]))
            filtcov.append(bedrow)
        endpoint = int(line[2])
        if mychr in chrsize:
            endpoints = chrsize[mychr]
            endpoints.append(endpoint)
            chrsize[mychr] = endpoints
        else:
            chrsize[mychr] = [endpoint]
    # Replace list of endpoints with max chr position
    for key, value in chrsize.items():
        maxpos = max(value)
        chrsize[key] = maxpos
    # The chrsize dict now contains chr lengths

    filtbed = BedTool(filtcov)
    # Merge regions to get loci
    # Distance of 100 should merge
    # properly paired reads with insert size < 300
    loci = filtbed.merge(d=100)
    # Counter list for all loci
    locicnt = 0
    # Counter list for chr loci
    chrcnts = []
    # Dic to hold chr loci stats
    for l in loci:
        locicnt = locicnt + 1
        scaff = l[0]
        start = int(l[1])
        stop = int(l[2])
        locmid = int((start + stop) / 2)
        if scaff in locidic:
            scaffcnt = locidic[scaff][0] + 1
            locpos = locidic[scaff][1]
            locpos.append(locmid)
            locidic[scaff] = [scaffcnt, locpos]
        else:
            locidic[scaff] = [1,[locmid]]

    for key, value in locidic.items():
             # Append loci per chr to list
             # Ignore unplaced and scaff
             if 'npla' not in key and 'scaff' not in key:
                chrcnts.append(value[0])
             # Get chr len from dict
             chrlen = chrsize[key]
             # Calculate max bin size with chr len
             toprange = math.ceil(chrlen / 100000.0) * 100000.0
             # Number of 100kb bins
             binnum = int(toprange / 100000.0)
             # Calculate histogram
             hist,bins = np.histogram(value[1],bins=binnum,range=(0, toprange))
             hist = np.ndarray.tolist(hist)
             mhist = round(statistics.mean(hist),4)
             sdhist = round(statistics.stdev(hist),4)
             lchrdic = {}
             lchrdic['Loci'] = value[0]
             lchrdic['Mean loci per 100kb'] = mhist
             lchrdic['StDev loci per 100kb'] =   sdhist
             locidic[key] = lchrdic
    try:
        lchrdic = {}
        lchrdic['Loci'] = locicnt
        lchrdic['Mean loci per chr'] = round(statistics.mean(chrcnts),4)
        lchrdic['StDev loci per chr'] = round(statistics.stdev(chrcnts),4)
        locidic['total'] = lchrdic
    except:
        lchrdic['Loci'] = 'NA'
        lchrdic['Mean loci per chr'] = 'NA'
        lchrdic['StDev loci per chr'] = 'NA'
        locidic['total'] = lchrdic

    for chrom in chrset:
        if chrom not in locidic:
            emptydic = {}
            emptydic['Loci'] = 0
            emptydic['Mean loci per 100kb'] = 0
            emptydic['StDev loci per 100kb'] = 0
            locidic[chrom] = emptydic
    print(locidic)
    return locidic
def covbreadth(intxt):
    '''Takes a mosdepth.dist.txt file as input
    and calculates the breadth of coverage per
    chromosome. Returns dict of chromosome breadths.
    '''
    bdic = {}
    chrset = set()
    with open(intxt,'r') as f:
        for line in f:
            line = line.strip()
            line = line.split('\t')
            chrom = line[0]
            chrset.add(chrom)
            if line[1] == '1':
                internal_dic = {}
                internal_dic['Breadth'] = string2num(line[2])
                bdic[line[0]] = internal_dic
    # Add zero-breadth chromosomes
    for chromosome in chrset:
        if chromosome not in bdic:
            internal_dic = {}
            internal_dic['Breadth'] = 0
            bdic[chromosome] = internal_dic
    return bdic


maindic = {}

inlist = sys.argv[1]
with open(inlist,'r') as f:
    for line in f:
        line = line.strip()
        infile = line
        idlist = file2id(infile)
        #try:
            # Sample ID
        sid = (idlist[2] + "_" + idlist[0] + "_" + idlist[1] +
               "_" + idlist[3] + "_" + idlist[4])
        # Get type of statistics file
        filetype = idlist[5]
        # Extract stats from file as dict
        if filetype == 'stats':
            stats = getmstats(infile)
            #print(sid, stats)
        elif filetype == 'perbase':
            meancov = perbasestats(infile)
            locicnt = radloci(infile)
            # Merge the nested dicts
            stats = mergedict(meancov,locicnt)
            #print(sid,allstats)
        elif filetype == 'mosdepth':
            stats = covbreadth(infile)

        if sid in maindic:
            new_entry = mergedict(maindic[sid],stats)
            maindic[sid] = new_entry
        else:
            maindic[sid] = stats
        #except:
        #    print("Error!")

dict2csv(maindic)
