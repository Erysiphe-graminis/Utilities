
readfile = open("/home/student/Jonathan/RsIIB4/Table_S2_Core_genes_list.tsv","r")

contents = readfile.readlines()

genelist = []

for line in contents[1:]:
    sline = line.strip().split("\t")
    genelist.append(sline[11][3:])
print("Identified " + str(len(genelist)) + " genes")
geneset = set(genelist)
if not len(geneset) == len(genelist):
    print("Warning, " + str(len(geneset)) + " unique gene names identified")
    genelist = list(geneset)
readfile.close()

readfile = open("/home/student/Jonathan/RsIIB4/GMI1000/GCF_000009125.1_ASM912v1_cds_from_genomic.fna", "r")
contents = readfile.readlines()
currentname = ""
currentsequence = ""
sequencetoextract = False
genedict = {}
for gene in genelist:
    genedict[gene] = ""
contents.append("\n>END_OF_FILE")
for line in contents:
    sline = line.strip()
    if line[0] == ">":
        sline = sline.split(" ")
        if sequencetoextract == True:
            genedict[currentname] = currentsequence
            sequencetoextract = False
            currentsequence = ""
        for entry in sline:
            if entry[:6] == "[locus":
                currentname = entry[14:-1]
        if currentname in geneset:
            sequencetoextract = True
    elif sequencetoextract == True:
        currentsequence = currentsequence + sline
readfile.close()

writefile = open("/home/student/Jonathan/RsIIB4/Core_Reference_genes.fa","w")
for gene in genelist:
    if len(genedict[gene]) < 1:
        print("Warning, gene ID " + gene + " has no associated sequence")
    else:
        writefile.write("\n>" + gene + "\n" + genedict[gene])
writefile.close()


