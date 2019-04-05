#!/usr/bin/env python
#this will merge identical values in the first column and sum the corresponding numbers in the second column
import sys
import pandas as pd
faidx = sys.argv[1]
        
#df = pandas.DataFrame.from_csv(faidx, sep='\t')
df = pd.read_csv(faidx, sep='\t', names = ["Scaff", "Len"])
sumgroup = df.groupby("Scaff").sum()
print(sumgroup)
