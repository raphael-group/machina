#!/usr/bin/python
import sys

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s <reads> <output_dir>\n" % sys.argv[0])
        sys.exit(1)

    reads_filename = sys.argv[1]
    output_dir = sys.argv[2]

    SAMPLE_LABEL_IDX = 1
    MUTATION_LABEL_IDX = 5
    REF_IDX = 6
    VAR_IDX = 7
    MULT = 5

    sample = ""
    ff = None
    with open(reads_filename) as f:
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        for line in f:
            s = line.rstrip("\n").split("\t")
            if sample != s[SAMPLE_LABEL_IDX]:
                if ff != None:
                    ff.close()
                sample = s[SAMPLE_LABEL_IDX]
                ff = open(output_dir + "/" + sample + ".tsv", "w")
                ff.write("\t".join(["mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn", "genotype"]) + "\n")
            pos = int(s[MUTATION_LABEL_IDX])
            for j,i in enumerate(range(pos*1000, pos*1000 + MULT)):
                if int(s[VAR_IDX]) >= MULT:
                    var = int(s[VAR_IDX]) + j - MULT / 2
                    ff.write("\t".join([str(i), s[REF_IDX], str(var), "2", "1", "1", "AB"]) + "\n")
                else:
                    ff.write("\t".join([str(i), s[REF_IDX], s[VAR_IDX], "2", "1", "1", "AB"]) + "\n")

        ff.close()
