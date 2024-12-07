import argparse
import numpy as np
import matplotlib.pyplot as plt

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--ids", type=str, help="Input tsv file with variant ID and respective AF value", required=True)
    parser.add_argument("-t", "--table", type=str, help="Intersection tsv file", required=True)
    parser.add_argument("-o", "--output", type=str, required=True, help="Output name for SAGA-only plot")
    parser.add_argument("-s", "--sample", type=str, required=True, help="Sample name (intersection name)")
    def list_of_ints(arg):
        return list(map(int, arg.split(',')))
    parser.add_argument("-c", "--callsets", type=list_of_ints, required=True, help="Values for callset indexes (for example: -c 2,4 represents an intersection of 4 callsets where the first 2 are from a patient)")

    args = parser.parse_args()

    return args

def plotting(id_file, table, output, sample, callsets):

    afs = dict()
    with open(id_file, 'r') as vcf:
        for line in vcf:
            columns = line.strip().split()
            if columns[1] != ".":
                afs[columns[0]] = float(columns[1])

    afs_saga = []
    afs_both = []

    idxs_comp = [9 + i for i in range(callsets[0] + 1)]
    idxs_add = [9 + i for i in range(callsets[0] + 1, callsets[1])]

    with open(table, 'r') as t:
        for line in t:
            columns = line.strip().split()
            if all([columns[i] == "True" for i in idxs_comp]):
                ids = columns[9 + callsets[0] + callsets[1]].split(";")
                if all([columns[i] == "False" for i in idxs_add]):
                    for j in ids:
                        if j in afs:
                            if afs[j] > 0:
                                afs_saga.append(afs[j])
                else:
                    for j in ids:
                        if j in afs:
                            if afs[j] > 0:
                                afs_both.append(afs[j])

    _, bins = np.histogram(afs_saga, bins=20)
    logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
    plt.hist(afs_saga, bins=logbins, edgecolor='black')
    plt.xscale('log')
    plt.title("SAGA-only AF distribution for sample {}".format(sample))
    plt.xlabel("Allele Frequency")
    plt.ylabel("Number of variants")
    plt.tight_layout()
    plt.savefig(output)
    plt.close()

    if callsets[1] - callsets[0] > 0:
        _, bins = np.histogram(afs_both, bins=20)
        logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        plt.hist(afs_both, bins=logbins, edgecolor='black')
        plt.xscale('log')
        plt.title("SAGA and HGSVC AF distribution for sample {}".format(sample))
        plt.xlabel("Allele Frequency")
        plt.ylabel("Number of variants")
        plt.tight_layout()
        out_both = output.replace("saga","both")
        plt.savefig(out_both)
        plt.close()


def main():
    args = parse_arguments()
    plotting(args.ids, args.table, args.output, args.sample, args.callsets)

if __name__ == "__main__":
    main()