

import os
import argparse

parser = argparse.ArgumentParser(description='Process the output of earlier scripts to combine muscle outputs')
parser.add_argument("-i", "--infile", dest="infile", help="path to a .tsv file where the names of accessions are in the first column, and the path to where their CDS .fa file is stored (inclusive) is in the fourth.",action="store")
parser.add_argument("-c", "--genefile", dest="cdpath", help="path to input the list of reference genes from", action="store")
parser.add_argument("-m", "--musclepath", dest="inpath", help="path to input the muscle alignments from \"./gene_alignments\"",action="store", default = "./gene_alignments")
parser.add_argument("-o", "--outfile", dest="outfile", help="file to output the concatenated alignment to, defaults to \"./Concatenated_Muscle_alignments.phy\"",action="store", default = "./Concatenated_Muscle_alignments.phy")
parser.add_argument("-a", "--namecol", dest="namecol", help="Override the column where names are stored in the inventory file",action="store", default = "1")
(options) = parser.parse_args()

if not os.path.isfile(options.infile):
    print("Can't find a file at " + options.infile + " exiting")
    exit()



if not os.path.isfile(options.cdpath):
    print("Can't find a file at " + options.cdpath + " exiting")
    exit()




genelist = []
geneset = set()
genedict = {}
accessions = []


#identify accessions
readfile = open(options.infile,"r")
contents = readfile.readlines()
readfile.close()
for line in contents[1:]:
    sline = line.strip().split("\t")
    accessions.append(sline[int(options.namecol)-1])

accessions.append("Reference")



#identify total gene list

readfile = open(options.cdpath,"r")
contents = readfile.readlines()
readfile.close()
for line in contents:
    sline = line.strip()
    if sline[0] == ">":
        geneID = sline.split(">")[1].split("...")[0].split(" ")[0]
        genelist.append(geneID)
        geneset.add(geneID)

if not len(genelist) == len(geneset):
    print("Error: Some genes are duplicated")



#read in gene sequences
for gene in genelist:
    genedict[gene] = {}
    if not os.path.isfile(options.inpath + "/" + gene + ".aln"):
        print("No alignment detected for gene " + gene)
        for accession in accessions:
            genedict[gene][accession] = ""
        continue
    readfile = open(options.inpath + "/" + gene + ".aln","r")
    contents = readfile.readlines()
    contents.append(">End")
    readfile.close()
    geneseq = ""
    accession = "noaccession"
    lengthset = set()
    for line in contents:
        sline = line.strip()
        if sline[0] == ">":
            if not accession == "noaccession":
                genedict[gene][accession] = geneseq
                lengthset.add(len(geneseq))
            for query in accessions:
                if sline[1:len(query)+1] == query:
                    accession = query
                    continue
            geneseq = ""
        else:
            geneseq = geneseq + sline
    if not len(lengthset) == 1:
        print("Something is wrong with the alignments for gene " + gene + " they are not the same length")
        print(lengthset)
#print them out again

writefile = open(options.outfile,"w")

for accession in accessions:
    finalseq = ""
    for gene in genelist:
        finalseq = finalseq + genedict[gene][accession]
    if accession == accessions[0]:
        writefile.write(str(len(accessions)) + " " + str(len(finalseq)))
    writefile.write("\n" + accession + " ")
    writefile.write(finalseq)

writefile.close()
