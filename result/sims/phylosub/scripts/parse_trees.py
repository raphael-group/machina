import sys
import json
from pandas import read_table

#filename='example_data.summ.json'
#edgelist_file='edgelist.txt'
#input_file='reads_seed0.tsv'
#leaflabel_file="leaf_labels.txt"

if len(sys.argv) < 5: 
    print "Usage: phylowgs_results machina_input_file output_edgelist output_leaflabel" 
    exit(1)

filename = sys.argv[1]
input_file = sys.argv[2]
edgelist_file = sys.argv[3]
leaflabel_file = sys.argv[4]

with open(filename) as f:
    doc = "".join(f.readlines())

parsed_json = json.loads(doc)

###
#   Find the best tree
###

best = float('-inf')
best_v = None
for v in parsed_json['trees']:
    ll = float(parsed_json['trees'][v]['llh'])
    if ll > best:
        best_v = v
        best = ll

struct = parsed_json['trees'][best_v]['structure']
      
###
#   Read in the tree file for the best tree to get mut assgmnts
###
result_dir = "/".join(filename.split('/')[:-1]) 
with open(result_dir+"/"+best_v+'.json') as f:
    tree = "".join(f.readlines())

tree_json = json.loads(tree)
mut_assignments = tree_json['mut_assignments']

####
#    Get usages for the tree
####
freqs = {}
for w in parsed_json['trees'][best_v]['populations']:
        cp = parsed_json['trees'][best_v]['populations'][w]['cellular_prevalence'] 
        freqs[w] = cp

bin_usage = {}
usage = {}
for w in parsed_json['trees'][best_v]['populations']:
    freq = freqs[w]
    if w in struct:
        children = struct[w]
        sum_children = [sum([freqs[str(c)][i] for c in children]) for i in range(len(freq))]
    else:
        sum_children = [0]*len(freq)
    bin_usage[w] = [int((f - s) > 0.05) for f, s in zip(freq, sum_children)]
    usage[w] = [(f - s) for f, s in zip(freq, sum_children)]
      

###
#   Get the set of vertices
###
vs = []
for v in struct:
    vs.append(int(v))
    for w in struct[v]:
        vs.append(int(w))
        
vs = set(vs)

###
#   Read in input file to get character labels and sites
###

character_map={}
site_map={}
data = read_table(input_file, skiprows=3)
for i,v in enumerate(data['character_index'].unique()):
    label = data[data['character_index']==v]['character_label'].iloc[0]
    character_map[v] = label

for i,v in enumerate(data['#sample_index'].unique()):
    label = data[data['#sample_index']==v]['anatomical_site_label'].iloc[0]
    site_map[v] =label


###
#   Name vertices according to what mutations occured in vertices
###
name_map = {}
for i in range(max(vs)+1):
    try:
        ssms = mut_assignments[str(i)]['ssms']
    except KeyError:
        continue

    ssm_names = [character_map[int(l[1:])] for l in ssms]
    name=";".join(map(lambda x: str(x), ssm_names))
    name_map[str(i)]=name

name_map['0'] = 'GL'


###
#   Write out the edge list for the tree
###

with open(edgelist_file, 'w') as edgelist, open(leaflabel_file, 'w') as leaf:
    for v in struct:
        for w in struct[v]:
            name_v = name_map[str(v)]
            name_w = name_map[str(w)]
            edgelist.write("\t".join((name_v, name_w))+'\n')
            
            # Leaf attachements
            samples = bin_usage[unicode(w)]
            sites = set([site_map[i] for i,s in enumerate(samples) if s])
            for s in sites:
                edgelist.write("\t".join((name_w, name_w+"_"+s))+"\n")
                leaf.write("\t".join((name_w+"_"+s, s))+"\n")

print "Writing:", edgelist_file
print "Writing:", leaflabel_file

