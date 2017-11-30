#!/usr/bin/python
import sys

def parse_loci(filename):
    result = []
    with open(filename) as f:
        for line in f:
            s = line.rstrip("\n").split("\t")
            chromosome = int(s[0])
            position = int(s[1])
            change = s[2]
            name = s[3]

            result.append((chromosome, position, change, name))
    return result

def parse_reads(filename):
    SAMPLE_LABEL_IDX = 1
    MUTATION_LABEL_IDX = 5
    REF_IDX = 6
    VAR_IDX = 7

    result_tot = []
    result_var = []
    samples = set()
    with open(filename) as f:
        nrAnatomicalSites = int(f.readline().rstrip("\n").split()[0])
        nrSamples = int(f.readline().rstrip("\n").split()[0])
        nrMutations = int(f.readline().rstrip("\n").split()[0])
        f.readline()
        for line in f:
            s = line.rstrip("\n").split("\t")
            if s[SAMPLE_LABEL_IDX] not in samples:
                samples.add(s[SAMPLE_LABEL_IDX])
                result_tot.append([s[SAMPLE_LABEL_IDX], []])
                result_var.append([s[SAMPLE_LABEL_IDX], []])
            ref = int(s[REF_IDX])
            var = int(s[VAR_IDX])
            idx = int(s[MUTATION_LABEL_IDX])
            result_tot[-1][-1].append(ref + var)
            result_var[-1][-1].append(var)
    return result_tot, result_var, nrMutations

def write_reads(f, reads, nrMutations, loci):
    nrSamples = len(reads)
    samples = [r[0] for r in reads]
    f.write("\t".join(["Chromosome", "Position", "Change", "Gene"] + samples) + "\n")
    f.write("\t".join(["#Purity", ".", ".", "."] + ["1.0"] * nrSamples) + "\n")
    for i in range(nrMutations):
        s = map(str, list(loci[i])) + [str(r[1][i]) for r in reads]
        s[3] = 'mut' + str(i)
        f.write("\t".join(s) + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s <reads> <loci>\n" % sys.argv[0])
        sys.exit(1)

    reads_filename = sys.argv[1]
    loci_filename = sys.argv[2]

    loci = parse_loci(loci_filename)

    result_tot, result_var, nrMutations = parse_reads(reads_filename)

    write_reads(sys.stdout, result_var, nrMutations, loci)
    write_reads(sys.stderr, result_tot, nrMutations, loci)

