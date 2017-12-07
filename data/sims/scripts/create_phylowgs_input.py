import sys
import pandas as pd


###
#   Formats the simulated data into the input format for phylowgs
###

if len(sys.argv) < 3:
    print "USAGE: python create_phylowgs_input.py [data_file] [output_directory]" 

filename = sys.argv[1]
outdir = sys.argv[2]
output_file = outdir+'/'+filename.split('/')[-1].split('.')[0]+'.phylowgs'

data = pd.read_table(filename, skiprows=3)
with open(output_file, 'w') as out:
    out.write("\t".join(['id', 'gene', 'a', 'd', 'mu_r', 'mu_v'])+'\n')
    for i,v in enumerate(data['character_label'].unique()):
        sub = data[data['character_label']==v]
        idx = 's'+str(i)
        gene = str(v)
        a = ','.join(map(str,sub['ref']))
        d = ','.join(map(str,sub['ref']+sub['var']))
        out.write("\t".join(map(str, [idx, gene, a, d, 0.999, 0.5]))+'\n')

