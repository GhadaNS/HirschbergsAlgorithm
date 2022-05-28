import argparse
import nmw
import hb


parser = argparse.ArgumentParser()
# Creating optional arguments
parser.add_argument("-t", action="store_true", help="Print i & j iterations")
parser.add_argument("-f", action="store_true", help="The input is a file")
parser.add_argument("-l", action="store_true", help="The file input sequences are lines & not characters")
# Creating positional arguments
parser.add_argument("gap", type=int, help="Parameter gap penalty \"g\" for the construction of F Score Matrix")
parser.add_argument("match", type=int, help="Parameter match \"m\" for the construction of F Score Matrix")
parser.add_argument("diff", type=int, help="Parameter diff \"d\" for the construction of F Score Matrix")
parser.add_argument("A", help="The 1st sequence or file to be aligned")
parser.add_argument("B", help="The 2nd sequence or file to be aligned")
# Getting Arguments
args = parser.parse_args()
g = args.gap
m = args.match
d = args.diff
t = args.t
l = args.l
f = args.f

align_nmw = nmw.nmw(g, m, d)
align_hb = hb.hb(g, m, d, t)

if args.f:
    with open(args.A, 'r') as A, open(args.B, 'r') as B:
        if l:  # a and b are lists of lines
            a = A.read().splitlines()
            b = B.read().splitlines()

            print("\nAlignments using the NEEDLEMAN-WUNSCH Alignment -----------\n")
            f = align_nmw.F(a, b)  # Creation of the score matrix "f"
            w = []
            z = []
            ww, zz = align_nmw.EnumerateAlignments_lines(a, b, f, w, z)  # Getting the list of all possible alignments
            for w, z in zip(ww, zz):  # By alignments as lists
                for i, j in zip(w, z):  # By lines in the alignments lists
                    if i == j:
                        print("=", i, "\n=", j)
                    else:
                        print("<", i, "\n>", j)
            
            print("\nAlignments using the HIRSCHBERG Alignment -----------------\n")
            ww, zz = align_hb.Hirschberg_lines(a, b)
            for w, z in zip(ww, zz):
                for i, j in zip(w, z):
                    if i == j:
                        print("=", i, "\n=", j)
                    else:
                        print("<", i, "\n>", j)
        else:  # a and b are strings
            a = A.read().splitlines()
            b = B.read().splitlines()
            a = ''.join(a)  # Joining the lines so the whole file can be aligned as a single sequence
            b = ''.join(b)

            print("\nAlignments using the NEEDLEMAN-WUNSCH Alignment -----------\n")
            f = align_nmw.F(a, b)
            w = ""
            z = ""
            ww, zz = align_nmw.EnumerateAlignments(a, b, f, w, z)
            for w, z in zip(ww, zz):
                print(w)
                print(z)

            print("\nAlignments using the HIRSCHBERG Alignment -----------------\n")
            ww, zz = align_hb.Hirschberg(a, b)
            alignments = zip(ww, zz)
            alignments = list(set(alignments))  # Eliminates duplicate Alignments
            for w, z in alignments:
                print(w)
                print(z)
else:  # A & B arguments are already strings
    a = args.A
    b = args.B

    print("\nAlignments using the NEEDLEMAN-WUNSCH Alignment -----------\n")
    f = align_nmw.F(a, b)
    w = ""
    z = ""
    ww, zz = align_nmw.EnumerateAlignments(a, b, f, w, z)
    for w, z in zip(ww, zz):
        print(w)
        print(z)

    print("\nAlignments using the HIRSCHBERG Alignment -----------------\n")
    ww, zz = align_hb.Hirschberg(a, b)
    alignments = zip(ww, zz)
    alignments = list(set(alignments))  # Eliminates duplicate Alignments
    for w, z in alignments:
        print(w)
        print(z)
