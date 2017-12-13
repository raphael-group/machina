import sys
import pandas as pd


###
#   Formats the simulated data into the input format for phylowgs
###

if len(sys.argv) < 3:
    print "USAGE: python create_pastri_input.py [data_file] [output_directory]" 

filename = sys.argv[1]
outdir = sys.argv[2]
output_file = outdir+'/'+filename.split('/')[-1].split('.')[0]+'.pastri'

cluster_file = sys.argv[3]
output_proposal=outdir+'/'+filename.split('/')[-1].split('.')[0]+'.pastri_prop'

data = pd.read_table(filename, skiprows=3)
A = "> A\n"
D = "> D\n"
n_mut = len(data['character_index'].unique())
n_sample = len(data['#sample_index'].unique())
A += str((n_mut, n_sample)) + '\n'
D += str((n_mut, n_sample)) + '\n'
with open(output_file, 'w') as out:
    for i,v in enumerate(data['character_index'].unique()):
        sub = data[data['character_index']==v]
        a = ' '.join(map(str,sub['var'])) + '\n'
        d = ' '.join(map(str,sub['ref']+sub['var'])) + '\n'
        A += a
        D += d
        #out.write("\t".join(map(str, [idx, gene, a, d, 0.999, 0.5]))+'\n')

    out.write(A+'\n')
    out.write(D)


v_dict = {}
with open(cluster_file) as f:
    for i,line in enumerate(f):
        line = map(int, line.strip().split(';'))
        for l in line: v_dict[l] = i

n_cluster = i+1
data['cluster'] = data['character_label'].apply(lambda x: v_dict[x])
cluster_data = data[['cluster', '#sample_index', 'ref', 'var']].groupby(['cluster','#sample_index']).sum().reset_index()


A = "> Alpha\n"
D = "> Beta\n"

A += str((n_cluster, n_sample)) + '\n'
D += str((n_cluster, n_sample)) + '\n'
for i,v in enumerate(cluster_data['cluster'].unique()):
            sub = cluster_data[cluster_data['cluster']==v]
            a = ' '.join(map(str,sub['var'])) + '\n'
            d = ' '.join(map(str,sub['ref']+sub['var'])) + '\n'
            #d = ' '.join(map(str,sub['ref'])) + '\n'
            A += a
            D += d

with open(output_proposal, 'w') as out:
    out.write(A+'\n')
    out.write(D)
