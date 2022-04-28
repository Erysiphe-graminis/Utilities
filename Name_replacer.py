

import argparse

parser = argparse.ArgumentParser(description='Takes an input file -i, and an index -x, and replaces all instances of string a with string b, then outputs -o')
parser.add_argument("-i", "--infile", dest="infile", help="path to a file to substitute the text within.",action="store")
parser.add_argument("-x", "--index", dest="indexfile", help="path to a .tsv file where the replacement strings can be found",action="store")
parser.add_argument("-a", "--column_a", dest="acol", help="column in the inxex file of original strings (0-indexed)",action="store")
parser.add_argument("-b", "--column_b", dest="bcol", help="column in the inxex file of original strings (0-indexed)",action="store")
parser.add_argument("-o", "--outfile", dest="outfile", help="path to the output file",action="store")
(options) = parser.parse_args()

columna = []
columnb = []

readfile = open(options.indexfile,"r")
contents = readfile.readlines()
for line in contents:
    sline = line.strip().split("\t")
    print("Now on line beginning " + sline[0])
    columna.append(sline[int(options.acol)])
    columnb.append(sline[int(options.bcol)])

readfile.close()

readfile = open(options.infile,"r")
writefile = open(options.outfile,"w")
contents = readfile.readlines()
for line in contents:
    rline = line
    for entry in range(len(columna)):
        rline = rline.replace(columna[entry],columnb[entry])
    writefile.write(rline)
readfile.close()
writefile.close()
