import os
import argparse

parser = argparse.ArgumentParser(description='Process the output of earlier scripts to add pfam data, where it exists to genes - returns a .tsv with extra columns (also adds the chromosome).')
parser.add_argument("-g", "--genes", dest="infile", help="path to a .tsv file where the first column is gene IDs.",action="store")
parser.add_argument("-a", "--annotations", dest="gfffile", help="path to a .gff file where these genes can be found",action="store")
parser.add_argument("-d", "--database", dest="pfamfile", help="path to a .txt file downloaded from the pfam website, listing all domains and their annotations",action="store")
parser.add_argument("-o", "--outfile", dest="outfile", help="path to the output file",action="store")
(options) = parser.parse_args()


geneset = set()
genelist = []
genedict= {}

readfile = open(options.infile, "r")
print("Reading in original file")
originalfile = readfile.readlines()
returnfile = []
for line in originalfile:
    sline = line.strip().split("\t")
    returnfile.append(sline)
    geneset.add(sline[0])

genelist = list(geneset)
for gene in genelist:
    genedict[gene] = {"pfam":[],"Annotation":[],"Designation":"","Chromosome":""}


readfile.close()
print("reading in pfam file")
pfamdatabase = {}
pfamset = set()
readfile = open(options.pfamfile, "r")
contents = readfile.readlines()
for line in contents[1:]:
    sline = line.strip().split("\t")
    if not sline[0] in pfamset:
        pfamdatabase[sline[0]] = sline[2]
        pfamset.add(sline[0])
print("added " + str(len(pfamset)) + " pfam entries")
readfile.close()
print("Reading in gff3 file")
readfile = open(options.gfffile, "r")
contents = readfile.readlines()
for line in contents[1:]:
    sline = line.strip().split("\t")
    if sline[2] == "mRNA":
        desc = sline[8]
        desc = desc.split(";")
        ID = desc[0][3:]
        if ID in geneset:
            genedict[ID]["Chromosome"] = sline[0]
            if len(desc) > 2:
                pfamlist = desc[3].replace(",","").split("PFAM:")[1:]
                genedict[ID]["pfam"] = pfamlist
                for entry in pfamlist:
                    if entry in pfamset:
                        genedict[ID]["Annotation"].append(pfamdatabase[entry])
                    else:
                        genedict[ID]["Annotation"].append("None found")
                        print("No pfam annotation found for entry " + entry + " in gene " + ID)
            genedict[ID]["Designation"] = desc[2]
readfile.close()
print("Processing")
for entry in returnfile[1:]:
    ID = entry[0]
    entry.append(genedict[ID]["Chromosome"])
    entry.append(genedict[ID]["Designation"])
    entry.append(",".join(genedict[ID]["pfam"]))
    entry.append(",".join(genedict[ID]["Annotation"]))
print("Writing")
writefile = open(options.outfile, "w")
writefile.write(str("\t".join((returnfile[0])) + "\t" + "Chromosome" + "\t" + "Designation" + "\t" + "pfam_domains" + "\t" + "pfam_annotations"))
for line in returnfile[1:]:
    writefile.write("\n" + str("\t".join(line)))
writefile.close()

