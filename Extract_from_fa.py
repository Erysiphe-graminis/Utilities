#Role of script is to provide a file with identifiers and extract those from a fasta
#TO DO - update to also permit .gff extraction, and .gff mediated extraction from genomes

import optparse

parser = optparse.OptionParser()
parser.add_option("-g", "--genes", dest="invfile", help="path to text file listing the gene IDs to be extracted, will take the first element of a .tsv file per line. Does not accept headers",action="store", default="na")
parser.add_option("-i", "--inventory", dest="genomefile", help="path to the genome file to be extracted from OR to an inventory file giving the genomes by name in column 1 and the output folder in column 4, toggle with the -m flag",action="store", default="na")
parser.add_option("-m", "--multiple-genomes", dest="multimode", help="Flag to tell the software that the sequences are to be extracted from more than one genome, genome inventory file should be written as a .tsv where the first column is the accession name and the fourth is the path to the accession, with a header",action="store_true", default="False")
parser.add_option("-p", "--phylip-out", dest="phylipmode", help="Flag to tell the software that the sequences are to be output in phylip (relaxed) format, otherwise a set of files, one for each genome, will be created",action="store_true", default="False")
parser.add_option("-o", "--outfile", dest="outfile", help="Name of the output file",action="store", default="Genes.fa")
parser.add_option("-n","--names",dest = "names", help="Prefix, if any, between the output files locations and accession names",default = "SNP_process_v10")
(options, args) = parser.parse_args()

invfile = options.invfile
accessions = []
accessions.append([])
accessions.append([])

if options.multimode == True:
    readfile = open(options.genomefile, "r")
    contents = readfile.readlines()
    readfile.close()
    print("Opening genome inventory file " + options.genomefile + "\nSkipping header line")
    for line in contents[1:]:
        sline = line.strip().split("\t")
        accessions[0].append(sline[0])
        accessions[1].append(sline[3] + "/" + options.names + sline[0] + "_full_CDS.fa")
else:
    accessions[0] = ["Accession"]
    accessions[1] = [options.genomefile]

print(accessions)

genedict = {}
genelist = []

readfile = open(invfile, "r")
contents = readfile.readlines()
readfile.close()
print("Reading in list of genes")
for line in contents:
    sline = line.strip().split("\t")[0]
    if len(sline)>0:
        if sline[0] == ">":
            sline = sline[1:]
        genelist.append(sline)
        genedict[sline] = {"Maxlen":0}
testset = set(genelist)

print("Done, found " + str(len(genelist)) + " genes")

for accession in range(len(accessions[0])):
    print(accession)
    readfile = open(accessions[1][accession], "r")
    accname = accessions[0][accession]
    print("Reading in from accession " + accname)
    contents = readfile.readlines()
    readfile.close()
    contents.append(">EOF")
    currentgene = ""
    currentseq = ""
    genesfound = set()
    for line in contents:
        sline = line.strip()
        if len(sline)>0:
            if sline[0] == ">":
                if not currentgene == "":
                    genedict[currentgene][accname] = currentseq
                    genesfound.add(currentgene)
                    if len(currentseq) > genedict[currentgene]["Maxlen"]:
                        genedict[currentgene]["Maxlen"] = len(currentseq)
                    currentseq = ""
                if sline[1:] in genelist:
                    currentgene = sline[1:]
                else:
                    currentgene = ""
            elif not currentgene == "":
                currentseq += sline
    tempset = testset - genesfound
    if len(tempset) > 0:
        print("Could not find " + str(tempset))

if options.phylipmode == True:
    writefile = open(options.outfile, "w")
    firsttime = True
    for accession in accessions[0]:
        tempseq = ""
        accname = accession
        print(accname)
        for gene in genelist:
            print(gene)
            tempseq += genedict[gene][accname]
        if firsttime == True:
            firsttime = False
            writefile.write(str(len(accessions[0])) + " " + str(len(tempseq)))
        writefile.write("\n" + accname + " " + tempseq)
    writefile.close()
else:
    for accession in range(len(accessions[0])):
        writename = options.outfile + "_" + accessions[0][accession]
        writefile = open(writename, "w")
        accname = accessions[0][accession]
        for gene in genelist:
            writefile.write("\n>" + gene + "\n" + genedict[gene][accname])
        writefile.close()





