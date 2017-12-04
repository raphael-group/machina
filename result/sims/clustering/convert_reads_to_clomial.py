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

    samples = set()
    positions = set()
    alt = {}
    tot = {}
    with open(reads_filename) as f:
        f.readline()
        f.readline()
        f.readline()
        f.readline()
        for line in f:
            s = line.rstrip("\n").split("\t")
            sample = s[SAMPLE_LABEL_IDX]
            samples.add(sample)
            if sample not in alt:
                alt[sample] = {}
                tot[sample] = {}

            pos = s[MUTATION_LABEL_IDX]
            positions.add(pos)
            alt[sample][pos] = s[VAR_IDX]
            tot[sample][pos] = str(int(s[VAR_IDX]) + int(s[REF_IDX]))

    with open(output_dir + "/total.tsv", "w") as f:
        f.write("\t".join(["position"] + list(samples)) + "\n")

        for pos in positions:
            f.write("chr1:" + pos)
            for sample in samples:
                f.write("\t" + tot[sample][pos])
            f.write("\n")

    with open(output_dir + "/alt.tsv", "w") as f:
        f.write("\t".join(["position"] + list(samples)) + "\n")

        for pos in positions:
            f.write("chr1:" + pos)
            for sample in samples:
                f.write("\t" + alt[sample][pos])
            f.write("\n")
