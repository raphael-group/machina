#!/usr/bin/python
import sys
import os

def parseSamples(filename):
    PATIENT_IDX = 0
    SAMPLE_IDX = 1
    ANATOMICAL_SITE_IDX = 7

    sampleToAnatomicalSite = {}
    with open(filename) as f:
        f.readline()
        for line in f:
            s = line.rstrip("\n").split("\t")
            s[ANATOMICAL_SITE_IDX] = s[ANATOMICAL_SITE_IDX].replace(" ", "_")
            sampleToAnatomicalSite[(s[PATIENT_IDX], s[SAMPLE_IDX])] = s[ANATOMICAL_SITE_IDX]

    return sampleToAnatomicalSite

def parsePatient(filename):
    with open(filename) as f:
        indexToSample = f.readline().rstrip("\n").split("\t")[2:]
        indexToCharacter = []
        characterPresences = []
        matrix = []
        for line in f:
            s = line.rstrip("\n").split("\t")
            characterLabel = s[0].replace(" ", "_")
            characterPresences.append(s[1].split(" "))
            indexToCharacter.append(characterLabel)
            values = map(float, s[2:])
            matrix.append(values)
        return indexToSample, indexToCharacter, characterPresences, matrix

if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.stderr.write("Usage %s <SAMPLES.txt> <PATIENT.txt> <PRIMARY> <width>\n" % sys.argv[0])
        sys.exit(1)

    samplesFilename = sys.argv[1]
    sampleToAnatomicalSite = parseSamples(samplesFilename)

    patientFilename = sys.argv[2]
    indexToSample, indexToCharacter, characterPresences, F = parsePatient(patientFilename)

    patient = os.path.basename(patientFilename).rstrip(".txt")

    primary = sys.argv[3]

    width = float(sys.argv[4])

    #print indexToSample
    #print indexToCharacter
    #print F
    #print sampleToAnatomicalSite

    indexToAnatomicalSite = []
    anatomicalSiteToIndex = {}
    for p, sample in enumerate(indexToSample):
        assert (patient, sample) in sampleToAnatomicalSite
        anatomicalSite = sampleToAnatomicalSite[(patient, sample)]
        if anatomicalSite not in indexToAnatomicalSite:
            if primary == anatomicalSite:
                for s in anatomicalSiteToIndex:
                    anatomicalSiteToIndex[s] += 1
                anatomicalSiteToIndex[primary] = 0
                indexToAnatomicalSite = [primary] + indexToAnatomicalSite
            else:
                anatomicalSiteToIndex[anatomicalSite] = len(indexToAnatomicalSite)
                indexToAnatomicalSite.append(anatomicalSite)

    m = len(indexToAnatomicalSite)
    k = len(indexToSample)
    n = len(indexToCharacter)

    print m, "#anatomical sites"
    print k, "#samples"
    print n, "#mutations" 
    print "\t".join(["#sample_index","sample_label","#anatomical_site_index","anatomical_site_label","character_index","character_label","f_lb","f_ub"])
    for p in range(k):
        sample = indexToSample[p]
        assert (patient, sample) in sampleToAnatomicalSite
        anatomicalSite = sampleToAnatomicalSite[(patient, sample)]
        assert anatomicalSite in anatomicalSiteToIndex
        for i in range(n):
            val = F[i][p]
            lb = min(1, max(val - width, 0))
            if sample not in characterPresences[i]:
                lb = 0
            ub = max(0, min(val + width, 1))
            print "\t".join([str(p), indexToSample[p], str(anatomicalSiteToIndex[anatomicalSite]), anatomicalSite, str(i), indexToCharacter[i], str(lb), str(ub)])
