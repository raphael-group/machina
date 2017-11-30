#!/usr/bin/python
import sys

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: %s <MACHINA_PMH_CTI_RESULT>\n" % sys.argv[0])
        sys.exit(1)

    trees = set()
    min_score = (100,100,100)
    with open(sys.argv[1]) as f:
        for line in f:
            s = line.rstrip("\n").split("\t")
            idx = int(s[0].rstrip("-"))
            if s[2] != '-':
                score = (int(s[2]), int(s[3]), int(s[4]))
                pattern = s[1].lstrip("(").rstrip(")").split(", ")[-1]
                if pattern != "R": continue
                if score < min_score:
                    min_score = score
                    trees = set()
                if score == min_score:
                    trees.add((idx, pattern, s[5]))

    for t in trees:
        print ",".join(map(str, [t[0], t[1], t[2]] + list(min_score)))
