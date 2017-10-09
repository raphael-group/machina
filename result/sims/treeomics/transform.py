import sys
import re

def parse(filename):
	res = set([])
	with open(filename) as f:
		for line in f:
			if line.find("\\item") != -1:
				s = line.rstrip("\n").split("\\item")
				for ss in s:
					g = re.match(r'.*{(.*) to (.*)}:.*', ss)
					if g:
						source = g.group(1).replace(" ", "_")
						target = g.group(2).replace(" ", "_")
						res.add((source, target))
	return res
				#print line, s

if __name__ == "__main__":
	if len(sys.argv) != 2:
		sys.stderr.write("Usage %s <FILE>\n" % sys.argv[0])

	filename = sys.argv[1]
	edge_list = parse(filename)
	sources = set()
	targets = set()
	for edge in edge_list:
		sources.add(edge[0])
		targets.add(edge[1])
		print edge[0], edge[1]

	leaves = targets.difference(sources)
	for leaf in leaves:
		sys.stderr.write(leaf + " " + leaf.split("_")[0] + "\n")
