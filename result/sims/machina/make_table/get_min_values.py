import pandas as pd
import sys 

filename = sys.argv[1]
table = pd.read_csv(filename, delimiter='\t', names = ['seed', 'patterns', 'mig', 'comig', 'seedsites', 'pattern', 'ub', 'lb', 'time'])
table = table[table['mig'].apply(lambda x: unicode(x).isnumeric())]
table = table[table['comig'].apply(lambda x: unicode(x).isnumeric())]
table['mig'] = table['mig'].apply(lambda x: int(x))
table['comig'] = table['comig'].apply(lambda x: int(x))

table['pattern'] = table['pattern'].apply(lambda x: x.replace('PS', 'S'))

for i, p in enumerate(table['patterns'].unique()):
    if i != 0: print " & "
    df = table[table['patterns'] == p]
    
    sortedDF = df.sort_values(['mig', 'comig'])
    minm = sortedDF.iloc[0]['mig']
    mincm = sortedDF.iloc[0]['comig']
    minVS = sortedDF[(sortedDF['mig'] == minm) & (sortedDF['comig'] == mincm)]['pattern'].unique()
    if len(minVS) == 1: print minVS[0]
    else:
        print "{"
        for i,v in enumerate(minVS):
            if i != 0: print ", "
            print v 
        print "}"
    print " & "
    
    #index = df['mig'].argmin()
    #print "(" + str(df.loc[index]['mig'])+", " + str(df.loc[index]['comig'])+")"
    print "(" + str(minm)+", " + str(mincm)+")"
    

