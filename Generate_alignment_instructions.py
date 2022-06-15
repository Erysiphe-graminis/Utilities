from optparse import OptionParser
import os
from datetime import datetime
now = datetime.now()
#Options for instructions - need to provide a master file detailing the input and output locations, as well as add flags for common choices (eg, limit outputs to 5) - also add functionality for copying over appropriate genome files?

parser = OptionParser()
parser.add_option("-i", "--inventory", dest="infile", help="path to inventory file detailing data to be fed into hisat2, and outputs",action="store",type="string")
parser.add_option("-s", "--short", dest="short", help="activate to only make a \"short\" file which does the first run", action="store_true",default=False)
parser.add_option("-o", "--outfile", dest="outfile", help="path to and name of bash file (without extension) to be created with the instructions given", action="store",type="string", default="./Alignment_instructions")
parser.add_option("-g", "--genomefile", dest="genomefile", help="Path to and name of the genome file (including extension) to be referenced by SNP calling", action="store",type="string", default="Nonegiven")
parser.add_option("-r", "--hisat2", dest="hisatindex", help="path to and name of the hisat2-index files to be referenced (without the.N.ht2 extension)", action="store",type="string", default="Nonegiven")
parser.add_option("-t", "--threads", dest="threads", help="threads to give hisat2 and pigz", action="store",type="string", default="1")
parser.add_option("-x", "--additional", dest="additional", help="additional options to pass to hisat, give them here \"enclosed in quotes\"", action="store",type="string", default="NA")
parser.add_option("-d", "--SNPdepth", dest="SNPdepth", help="a minimum depth requirement for SNPs, default = 3", action="store",type="string", default="3")
parser.add_option("-c", "--compress", dest="gzip", help="pigz output files, default = off", action="store_true", default=False)
parser.add_option("--sam-only", dest="samonly", help="only produce .sam files, no .bam or other processing, default = off", action="store_true", default=False)
parser.add_option("--keep-bam", dest="keepbam", help="keep the .bam file, as well as the pileups, default = off", action="store_true", default=False)
parser.add_option("--keep-sam", dest="keepsam", help="keep the .sam file, as well as the pileups, default = off", action="store_true", default=False)

(options, args) = parser.parse_args()

additionaloptions = "NA"
if not options.additional == "NA":
    additionaloptions = " " +options.additional

if not os.path.exists(options.infile):
    print("No such input file at "+options.infile)
    print("Script will continue, but proceed at your own risk")
if not os.path.exists(options.genomefile):
    print("No such genome file at "+options.genomefile)
    print("Script will continue, but proceed at your own risk")
if not os.path.exists(options.hisatindex + ".1.ht2"):
    print("No such index at "+options.hisatindex)
    print("Script will continue, but proceed at your own risk")

SNPdepth = options.SNPdepth

print("Reading inventory")
readfile = open(options.infile, "r")
contents = readfile.readlines()
readfile.close()

inventorylist = []
inventory = {}


#Go through each member of the inventory, harvest relevant information, store it in a nested library (can update model with more info later if necessary) - also creates a list of the files in the order they went into the inventory - can be sorted if desired
linelength = len(contents[0].split("\t"))
print("Skipping header line")
for count, line in enumerate(contents[1:]):
    if not len(line.split("\t")) == linelength:
        print("Error, line " + str(count) + " has a different number of components to the header, skipping")
        continue
    sline = line.strip().split("\t")
    startdir = sline[1]
    enddir = sline[2]
    name = sline[0]
    inventorylist.append(name)
    inventory[name] = {"startdir":startdir,"enddir":enddir,"name":name,}
print("Identified " + str(len(inventorylist)) + " accessions in file")
timestamp = now.strftime("_%Y_%m_%d_%H_%M")
outfile = options.outfile + timestamp + ".sh"
print("Writing to outfile " + outfile) 
writefile = open(outfile,"w")



maxaccessions = len(inventorylist)
if options.short == True:
    maxaccessions = 1
print("Generating instructions for " + str(maxaccessions) + " accessions")
writefile.write("#!/bin/bash\n")


valid_nuc_extensions = ["fa","fq","fasta","fastq"]

for accession in inventorylist[0:maxaccessions]:
    print(accession)
    filelist = os.listdir(inventory[accession]["startdir"])
    forwardlist = []
    reverselist = []
    unassignedlist = []
    finalname = inventory[accession]["enddir"] + "/" + inventory[accession]["name"]
    finaldir = inventory[accession]["enddir"] + "/"
    copylist = []
    for entry in filelist:
        splentry = entry.split(".")
        if len(splentry) > 1:
            compression = False
            if splentry[-1] == "gz":
                compression = True
            if splentry[-(1+ int(compression))] in valid_nuc_extensions:
                copylist.append(entry)
                if splentry[-(2+int(compression))][-1] in ("1","F"):
                    forwardlist.append(finaldir  + ".".join(splentry[:-int(compression)]))
                elif splentry[-(2+int(compression))][-1] in ("2","R"):
                    reverselist.append(finaldir  + ".".join(splentry[:-int(compression)]))
                else:
                    unassignedlist.append(finaldir  + ".".join(splentry[:-int(compression)]))
            else:
                print("Excluding file " + entry + " as file extensions are not in the recognized list (.fa, .fq, .fasta, .fastq + optionally .gz)")
    forwardlist.sort()
    reverselist.sort()
    unassignedlist.sort()
    totallist = forwardlist + reverselist + unassignedlist
    writefile.write("\nmkdir "+ finaldir)
    for entry in copylist:
        writefile.write("\ncp " + inventory[accession]["startdir"] + "/" + entry + " " + finaldir)
    writefile.write("\npigz -d -p 10 "+finaldir+"*")
    forwardlist = ",".join(forwardlist)
    reverselist = ",".join(reverselist)
    unassignedlist = ",".join(unassignedlist)
    writefile.write("\nhisat2 -q --no-unal --summary-file " + finaldir + inventory[accession]["name"]+"_hisat2_summary.txt --new-summary -p " + options.threads + " -x " + options.hisatindex + " -1 " + forwardlist + " -2 " + reverselist + " -S "+ finalname + "_aligned.sam")
    if len(unassignedlist) > 0:
        writefile.write(" -U " + unassignedlist)
    if not additionaloptions == "NA":
        writefile.write(additionaloptions)
    if not options.samonly == True:
        writefile.write("\nsamtools view -@ " +options.threads + " -S -b "+ finalname + "_aligned.sam > " +finalname + "_aligned.bam\nsamtools sort -@ " + options.threads + " " + finalname + "_aligned.bam -o " + finalname + "_aligned.sorted.bam\nsamtools mpileup -o " + finalname + "_aligned.sorted.bam.pileup -f " + options.genomefile + " " + finalname + "_aligned.sorted.bam")
        writefile.write("\nvarscan pileup2snp " + finalname + "_aligned.sorted.bam.pileup --min-coverage " + SNPdepth + " > " + finalname + "_aligned.sorted.bam.pileup2snp")
        writefile.write("\nvarscan pileup2indel " + finalname + "_aligned.sorted.bam.pileup --min-coverage 5 > " + finalname + "_aligned.sorted.bam.pileup2indel")
        if options.gzip == True:
            writefile.write("\npigz --best -p "+ options.threads + " " + finalname + "_aligned.sorted.bam.pileup")
        if not options.keepsam == True:
            writefile.write("\nrm "+finalname+"_aligned.sam")
        elif options.gzip == True:
            writefile.write("\npigz --best -p " + options.threads + " "+finalname+"_aligned.sam")
        writefile.write("\nrm "+finalname+"_aligned.bam")
        if not options.keepbam == True:
            writefile.write("\nrm "+finalname+"_aligned.sorted.bam")
        elif options.gzip == True:
            writefile.write("\npigz --best -p " + options.threads + " "+finalname+"_aligned.sorted.bam")
    for entry in totallist:
        writefile.write("\nrm "  + entry)
writefile.close()
