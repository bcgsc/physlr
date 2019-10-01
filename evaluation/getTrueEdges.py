import sys
import argparse
import timeit
import os

def areOverlapping(r, s):
    overlap = max(0, min(r[1], s[1]) - max(r[0], s[0]))
    if overlap == 0:
        return (False, overlap)
    return (True, overlap)

def main():

    argparser = argparse.ArgumentParser()
    argparser.add_argument(
        "-t", "--thread", action="store", dest="threads", type=int,
        default=min(16, os.cpu_count()),
        help="number of threads [16 or number of CPU]")
    argparser.add_argument(
        "-i", "--input", action="store", dest="input", type=str,
        help="input tigmint molecule bed file")
    argparser.add_argument(
        "-o", "--output", action="store", dest="output", type=str,
        help="ouput edge path")
    args = argparser.parse_args()
    t0 = timeit.default_timer()
    molecule = {}
    f = open(args.input, "r")
    h = open(args.output, "w")

    print(
        int(timeit.default_timer() - t0),
        "Collecting barcodes",\
            file=sys.stderr)

    for i in f:
        columns = i.split("\t")
        new_mol = columns[3]
        chromosome = columns[0]
        if chromosome not in molecule:
            molecule[chromosome] = {}
        while new_mol in molecule[chromosome]:
            new_mol = new_mol + "_"
        molecule[chromosome][new_mol] = (int(columns[1]), int(columns[2]))

    print(
        int(timeit.default_timer() - t0),
        "Intersecting barcodes",\
            file=sys.stderr)

    for chromosome in molecule:
        print(
            int(timeit.default_timer() - t0),
            "chromosome :", chromosome,\
                file=sys.stderr)
        molecule_keys = list(molecule[chromosome].keys())
        num_molecules = len(molecule_keys)
        print(
            int(timeit.default_timer() - t0),
            "# molecule in chromosome :", num_molecules,\
                file=sys.stderr)
        checkpoint = int(num_molecules / 10)
        checkpoint_num = 1
        for i in range(num_molecules):
            if (i + 1) % checkpoint == 0:
                print(
                    int(timeit.default_timer() - t0),
                    "Done", checkpoint_num * 10, "%",\
                        file=sys.stderr)
                checkpoint_num = checkpoint_num + 1
            molecule_n = molecule_keys[i]
            n = molecule[chromosome][molecule_n]
            for j in range(i+1, len(molecule_keys)):
                molecule_m = molecule_keys[j]
                m = molecule[chromosome][molecule_m]
                overlapResults = areOverlapping(n, m)
                if overlapResults[0]:
                    h.write(str(molecule_n).rstrip("_")+"\t"+str(molecule_m).rstrip("_")+"\t"+str(overlapResults[1])+"\n")
    f.close()
    h.close()

if __name__ == "__main__":
    main()
