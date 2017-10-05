#!/usr/bin/python
import sys

def parse_CNA(filename):
    intervals = []
    with open(filename) as f:
        for line in f:
            if line.startswith("#"): continue
            line = line.strip().split(",")
            sample = line[0]
            if sample.split("-")[0] != "A22": continue
            chrm = line[1]
            if chrm == "X": continue
            startpos = line[2]
            endpos = line[3]
            maj1 = line[8]
            min1 = line[9]
            frac1 = line[18]
            maj2 = line[11]
            min2 = line[12]
            frac2 = line[19]

            interval = [sample, chrm, int(startpos), int(endpos), (int(maj1), int(min1), float(frac1))]

            #if maj2 != "NA":
            #    interval.append((int(maj2), int(min2), float(frac2)))

            # only consider heterozygous diploid regions
            if maj1 == "1" and min1 == "1" and maj2 == "NA":
                intervals.append(interval)
    return intervals

def parse_SNV(filename):
    SNVs = {}
    with open(filename) as f:
        for line in f:
            if line.startswith("#"): continue
            line = line.strip().split(',')

            try:
                patient, sample, chrm, pos = line[:4]
                status = line[14]
                cluster = line[13]

                if patient != "A22": continue
                if sample.upper() == "A22-L": continue
                if chrm == "X": continue
                if status == "not-validated": 
                    vaf = float(line[6]) # WGS
                    depth = int(line[7]) # WGS
                else:
                    vaf = float(line[11]) # targeted
                    depth = int(line[12]) # targeted

                pos = int(pos)
                annotation = line[-1].split("|")[0]

            except:
                print line
                raise
            if (chrm, pos) not in SNVs:
                SNVs[(chrm,pos)] = []

            SNVs[(chrm,pos)].append([sample, vaf, depth, annotation, cluster])
    return SNVs

def match(CNAs, SNVs):
    interval_to_SNVs = {}
    SNVs_matched = {}
    for locus in SNVs:
        chrm, pos = locus
        SNVs_matched[(chrm, pos)] = []
        for snv_sample in SNVs[locus]:
            sample, vaf, depth, annotation, cluster = snv_sample
            found = False
            for cna in CNAs:
                int_sample, int_chrm, startpos, endpos, C1 = cna[:5]
                if len(cna) > 5: 
                    C2 = cna[5]
                    Cs = [C1,C2]
                else:
                    Cs = [C1]

                if sample.upper() == int_sample.upper() and chrm == int_chrm and pos > startpos and pos <= endpos:
                    SNVs_matched[locus].append(snv_sample + Cs)
                    found = True

                    if sample.upper() not in interval_to_SNVs:
                        interval_to_SNVs[sample.upper()] = {}
                    t = tuple(cna[1:4])
                    if t not in interval_to_SNVs[sample.upper()]:
                        interval_to_SNVs[sample.upper()][t] = []
                    interval_to_SNVs[sample.upper()][t].append(snv_sample[3])
                    break
        if SNVs_matched[locus] == []:
            del SNVs_matched[locus]


    return SNVs_matched

def writeSciClone(SNVs_matched):
    # get samples
    samples = {}
    for locus in SNVs_matched:
        for call in SNVs_matched[locus]:
            sample = call[0]
            samples[sample] = open("sciclone/" + sample + ".tsv", "w")
            #chr, pos, ref_reads, var_reads, vaf
            samples[sample].write("\t".join(["chr", "pos", "ref_reads", "var_reads", "vaf"]) + "\n")
        break

    for locus in SNVs_matched:
        chrm, pos = locus
        for call in SNVs_matched[locus]:
            sample = call[0]
            vaf = call[1]
            total_reads = call[2]
            var_reads = int(vaf * total_reads)
            ref_reads = int((1 - vaf) * total_reads)

            f = samples[sample]
            #chr, pos, ref_reads, var_reads, vaf
            f.write("\t".join(map(str, [chrm, pos, ref_reads, var_reads, vaf * 100])) + "\n")

    for sample in samples:
        f = samples[sample]
        f.close()

def writePyClone(SNVs_matched):
    # get samples
    samples = {}
    for locus in SNVs_matched:
        for call in SNVs_matched[locus]:
            sample = call[0]
            samples[sample] = open("pyclone/" + sample + ".tsv", "w")
            #chr, pos, ref_reads, var_reads, vaf
            samples[sample].write("\t".join(["mutation_id", "ref_counts", "var_counts", "normal_cn", "minor_cn", "major_cn"]) + "\n")
        break

    for locus in SNVs_matched:
        chrm, pos = locus
        for call in SNVs_matched[locus]:
            sample = call[0]
            vaf = call[1]
            total_reads = call[2]
            symbol = call[3]
            var_reads = int(vaf * total_reads)
            ref_reads = int((1 - vaf) * total_reads)
            identifier = "%s_chr%s:%d" % (symbol, chrm, pos)

            f = samples[sample]
            #chr, pos, ref_reads, var_reads, vaf
            f.write("\t".join(map(str, [identifier, ref_reads, var_reads, 2, 1, 1])) + "\n")

    for sample in samples:
        f = samples[sample]
        f.close()

def writeMACHINA(SNVs_matched):
    # get samples
    samples = {}
    for locus in SNVs_matched:
        for call in SNVs_matched[locus]:
            sample = call[0]
            samples[sample] = len(samples)
        break

    mutation = {}
    print "\t".join(["#sample_index", "sample_label", "anatomical_site_index", "anatomical_site_label", "character_index", "character_label", "ref", "var"])
    for locus in SNVs_matched:
        chrm, pos = locus
        if len(SNVs_matched[locus]) < 10: continue
        for call in SNVs_matched[locus]:
            sample = call[0]
            vaf = call[1]
            total_reads = call[2]
            symbol = call[3]
            var_reads = int(vaf * total_reads)
            ref_reads = int((1 - vaf) * total_reads)
            identifier = "%s:chr%s:%d" % (symbol, chrm, pos)
            if identifier not in mutation:
                mutation[identifier] = len(mutation)

            print "\t".join(map(str, [samples[sample], sample, 
                samples[sample], sample, mutation[identifier], identifier, ref_reads, var_reads]))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.stderr.write("Usage: %s <COPY_NUMBER.csv> <VARIANTS.csv>\n" % sys.argv[0])
        sys.exit(1)

    CNAs = parse_CNA(sys.argv[1])
    #for cna in CNAs:
    #    print cna

    SNVs = parse_SNV(sys.argv[2])
    #for locus in SNVs:
    #    print locus, SNVs[locus]

    SNVs_matched = match(CNAs, SNVs)
    #writeSciClone(SNVs_matched)
    #writePyClone(SNVs_matched)

    writeMACHINA(SNVs_matched)
