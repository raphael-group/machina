#!/usr/bin/python
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: %s <reads>\n" % sys.argv[0])
        sys.exit(1)

    reads_filename = sys.argv[1]

    SAMPLE_LABEL_IDX = 1
    MUTATION_LABEL_IDX = 5
    REF_IDX = 6
    VAR_IDX = 7

    with open(reads_filename) as f:
        m = int(f.readline().split()[0])
        k = int(f.readline().split()[0])
        n = int(f.readline().split()[0])
        f.readline()
        sample = ""
        ok = False
        string = ""
        print " ".join(["sample"] + map(str, range(n)))
        for line in f:
            s = line.rstrip("\n").split("\t")
            if sample != s[SAMPLE_LABEL_IDX]:
                sample = s[SAMPLE_LABEL_IDX]
                if ok:
                    print string
                string = sample
                ok = False

            vaf = float(s[VAR_IDX]) / (float(s[VAR_IDX]) + float(s[REF_IDX]))
            if vaf >= 0.01:
                string += " " + str(1)
                ok = True
            else:
                string += " " + str(0)
        if ok:
            print string
