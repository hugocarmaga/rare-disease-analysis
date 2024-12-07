import argparse
import matplotlib.pyplot as plt
import numpy as np

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--gaf", type=str, help="GAF file with alignments to plot (Sample name should be in the file name as in 'sample_something.gaf')", required=True)
    parser.add_argument("-o", "--output", required=True, help="Output file name.")
    parser.add_argument("-s", "--sample", type=str, required=True, help="Sample name (intersection name)")

    args = parser.parse_args()

    return args

def plotting(gaf, output, sample):
    idents = []

    with open(gaf,'r') as g:
        for line in g:
            columns = line.strip().split()
            idents.append(float(columns[15].split(":")[2]))

    _, _, bars = plt.hist(idents, bins=np.arange(0.85,1,0.01), edgecolor='black')
    plt.xlabel("Sequence identity")
    plt.ylabel("Number of variants")
    plt.title("Sequence identity distribution for sample {} alignments".format(sample))
    plt.text(plt.xlim()[0] * 1.01, plt.ylim()[1] * 0.9, 'Total of {} variants'.format(len(idents)), fontsize = 12, bbox = dict(facecolor = 'white', alpha = 0.5))
    plt.bar_label(bars)
    plt.savefig(output)

def main():
    args = parse_arguments()
    plotting(args.gaf, args.output, args.sample)

if __name__ == "__main__":
    main()