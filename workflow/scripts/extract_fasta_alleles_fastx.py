import argparse
import gzip
from pyfastx import Fasta

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", type=str, help="FASTA file with the reference sequence", required=True)
    parser.add_argument("-v", "--vcf", type=str, help="VCF file with the variants to extract", required=True)
    parser.add_argument("-s", "--shoulder", type=int, help="Number of bp to add to the ends of the variant allele", required=False, default=1000)
    parser.add_argument("-o", "--out", type=str, help="Output prefix for the FASTA file with the variant allele sequences ", required=True)

    args = parser.parse_args()

    return args

def extract_vars(fasta, vcf, shoulder, out):
    o = open(out, 'w')
    with gzip.open(vcf,'rt') as v:
        for line in v:
            i=1
            if not line.startswith("#"):
                columns = line.split()

                header = columns[2]
                pos = int(columns[1])
                start = max(0, pos - shoulder - 1)

                svtype = columns[7].split(";")[1].split("=")[1]
                if svtype == "INS":
                    end = min(len(fasta[columns[0]]), pos + shoulder)

                    sequence = fasta.fetch(columns[0],(start, pos)) + columns[4][1:] + fasta.fetch(columns[0], (pos, end))

                elif svtype == "DEL":
                    size = len(columns[3])
                    end = min(len(fasta[columns[0]]), pos + shoulder + size)

                    sequence = fasta.fetch(columns[0],(start, pos)) + fasta.fetch(columns[0], (pos + size, end))

                o.write(">" + header + "\n")
                o.write(sequence + "\n")

            if i % 100 == 0:
                print("Read {} lines".format(i))
            i += 1

    o.close()


def main():
    args = parse_arguments()
    fasta = Fasta(args.fasta)
    extract_vars(fasta, args.vcf, args.shoulder, args.out)

if __name__ == "__main__":
    main()


