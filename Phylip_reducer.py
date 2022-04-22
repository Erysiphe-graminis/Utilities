#Takes a .phy (relaxed input) and processes it to remove uninformative positions based on either or both of depth and diversity. Currently takes nucleotides only


import optparse

parser = optparse.OptionParser(description='Takes an input file in phylip format (relaxed) and outputs a file with positions reduced to only those which are informative')
parser.add_option("-i", "--infile", dest="infile", help="path to a file to substitute the text within.",action="store")
parser.add_option("-o", "--outfile", dest="outfile", help="path to the output file",action="store")
parser.add_option("-v", "--variants", dest="variants", help="Minimum different nucleotides necessary at a position. Default of 2. Raising this will exponentially reduce the size of the file",action="store", default = "2")
parser.add_option("-b", "--upperbound", dest="bound", help="Maximum different nucleotides necessary at a position. Defaults to off (-1). Giving 2, for example, will remove a position where the options are A, T and W",action="store", default = "-1")
parser.add_option("-d", "--depth", dest="depth", help="Minimum depth of a position to include it (ie, number of non-missing nucleotides)",action="store", default = "0")
parser.add_option("-y", "--variability", dest="SNPs", help="Minimum number of accessions with a SNP for the position to be considered. Operates by comparing each SNP to the total number of accessions and removing positions where no output is over the given threshold (eg, for 10 accessions, 8/2 both score 2, but 7/2/1 score 3/2/1)",action="store", default = "1")
parser.add_option("-x", "--exclude", dest="exclude", help="Comma seperated list of accessions to ignore when deciding on positions to keep, IE if every accession but one has the same \"SNP\" then that SNP would be passed over if that accession is listed here. The accessions will still be read in and printed out",action="store", default = "")
parser.add_option("-s", "--strict", dest="strict", help="Output in strict Phlylip mode (all sequences are length 10",action="store_true", default = False)
(options, args) = parser.parse_args()

#hardcoded for now, fix later

min_vars = int(options.variants)
min_present = int(options.depth)
upperbound = int(options.bound)
lowerbound = int(options.SNPs)
excludelist = options.exclude.split(",")

nucleotide_pairs = {"A":"T","T":"A","G":"C","C":"G","R":"Y","Y":"R","S":"S","W":"W","K":"M","M":"K","B":"V","D":"H","H":"D","V":"B","N":"N",".":"-","-":"-","a":"t","t":"a","g":"c","c":"g","r":"y","y":"r","s":"s","w":"w","k":"m","m":"k","b":"v","d":"h","h":"d","v":"b","n":"n"}

nucleotides = nucleotide_pairs.keys()
nonsynonymousnucleotides = ["A","T","G","C","R","Y","S","W","K","M","B","V","D","H"]

readfile = open(options.infile, "r")
contents = readfile.readlines()
readfile.close()
line1= contents[0].strip().split(" ")
datasets = int(line1[0])
characters = int(line1[1])
#reads in file
print("Read in a file with " + str(characters) + " nucleotides and " + str(datasets) + " accessions") 
print("Evaluating positions")
positionlist = []
nonsynonymouscount = []
poscovcount = []
for value in range(characters+1):
    positionlist.append(set())
    nonsynonymouscount.append(0)
    poscovcount.append({})
    for nucl in nonsynonymousnucleotides:
        poscovcount[value][nucl] = 0

#indexes all positions
for line in contents[1:]:
    if not line[:-characters+1].strip() in excludelist:
        for value, nucleotide in enumerate(line.strip()[-(characters+1):]):
            nucl = nucleotide.upper()
            positionlist[value].add(nucl)
            if nucl in nonsynonymousnucleotides:
                nonsynonymouscount[value] += 1
                poscovcount[value][nucl] += 1
if not len(contents)-1 == datasets:
    print("Error! You have " + str(len(contents)-1) + " datasets, but the file header says there are " + str(datasets) + " You will need to fix this")

keptpositions = []
#compares to each option and keeps positions which fulfill all requirements

for position, nucleotide in enumerate(positionlist):
    nucleotide.discard("N")
    nucleotide.discard("-")
    if len(nucleotide) >= min_vars:
        if len(nucleotide) <= upperbound or upperbound == -1:
            if nonsynonymouscount[position] >= min_present:
                for nucleo in nonsynonymousnucleotides:
                    #print((datasets/2)-(abs((datasets/2)-(poscovcount[position][nucleo]))))
                    if (datasets/2)-(abs((datasets/2)-(poscovcount[position][nucleo]))) >= lowerbound:
                        keptpositions.append(position)
                        break

#write outputs
print("Kept " + str(len(keptpositions)) + " positions out of " + str(characters))
writefile = open(options.outfile, "w")
writefile.write(str(datasets) + "\t" + str(len(keptpositions)))
for line in contents[1:]:
    line = line.strip()
    print("writing output for accession " + line[:-(characters)])
    if options.strict == False:
        writefile.write("\n" + line[:-(characters)])
    else:
        writefile.write("\n")
        newname = line[:-characters]
        if len(newname) > 10:
            newname = newname[:10]
        while len(newname) < 10:
            newname = newname + " "
        writefile.write(newname)
    for position in keptpositions:
        writefile.write(line[-(characters-position+1)])
writefile.close()
