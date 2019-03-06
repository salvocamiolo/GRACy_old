import os
import sys
import time
#Note!
#Remove -test in the webin command in order to upload in the official website
#The script should be launched within the folder containing the fw files



inputFolder = "/home3/scc20x/Desktop/NicoCollectionReads/UGompelsReads/UrsulaResubmission" #Change from Gui
filesToSubmit = "./toSubmit.txt" #Change from Gui
createProject = "no" #Change from Gui it may be no and be present in the readsInfo file
createSample = "no" #Change from Gui. It may be no and be present in the readsInfo file
projectInfo = "projectInfo.txt"
sampleInfo = "samplesInfo.txt"
readsInfo = "readsInfo.txt"

sampleAccession = {}
runAccession = {}
projectAccessionNumber = ""
experimentAccession = {}
sampleList = []



def createNewProject(projectFile):
    if os.path.isfile(projectFile)==True:
        projectInfoFile = open(projectFile)
        projectInfoFile.readline()
        line = projectInfoFile.readline().rstrip()
        projectInfoFields = line.split("\t")
        projectAlias = projectInfoFields[0]
        projectName = projectInfoFields[1]
        projectTitle = projectInfoFields[2]
        projectDescription = projectInfoFields[3]
        project_xml = open("project.xml","w")
        project_xml.write("<PROJECT_SET>\n")
        project_xml.write("<PROJECT alias=\""+projectAlias+"\">\n")
        project_xml.write("<NAME>"+projectName+"</NAME>\n")
        project_xml.write("<TITLE>"+projectTitle+"</TITLE>\n")
        project_xml.write("<DESCRIPTION>"+projectDescription+"</DESCRIPTION>\n")
        project_xml.write("<SUBMISSION_PROJECT>\n<SEQUENCING_PROJECT></SEQUENCING_PROJECT>\n</SUBMISSION_PROJECT>\n<PROJECT_LINKS>\n<PROJECT_LINK>\n<XREF_LINK>\n<DB></DB>\n<ID></ID>\n</XREF_LINK>\n</PROJECT_LINK>\n</PROJECT_LINKS>\n</PROJECT>\n</PROJECT_SET>")
        #IMPORTANT! To change username and password in GRACy
        project_xml.close()
        
        print "curl -u Webin-50760:Elisaegizia14 -F \"SUBMISSION=@submission.xml\" -F \"PROJECT=@project.xml\" \"https://www.ebi.ac.uk/ena/submit/drop-box/submit/\" >projectReceipt"
        os.system("curl -u Webin-50760:Elisaegizia14 -F \"SUBMISSION=@submission.xml\" -F \"PROJECT=@project.xml\" \"https://www.ebi.ac.uk/ena/submit/drop-box/submit/\" >projectReceipt")
        receiptFile = open("projectReceipt")

        while True:
            line = receiptFile.readline().rstrip()
            if not line:
                break
            fields = line.split("\"")
            for item in fields:
                if "PRJE" in item:
                    receiptFile.close()
                    return item
                    break

        


def createNewSample(sampleName):
    if not sampleName in sampleAccession:
        sampleAccession[sampleName] = ""
    infile = open(sampleInfo)
    infile.readline() #Read header

    sampleFound = 0
    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        alias = fields[0]
        center = fields[1]
        title = fields[2]
        taxonID = fields[3]
        scientificName = fields[4]
        geographicLocation = fields[5]
        hostCommonName = fields[6]
        hostHealthState = fields[7]
        isolationSource = fields[8]
        sex = fields[9]
        hostScientificName = fields[10]
        collectorName = fields[11]
        collectingInstitution = fields[12]
        isolate = fields[13]
        
        if alias == sampleName:
            sampleFound = 1
            xmlfile = open(alias+"_sample.xml","w")
            xmlfile.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<SAMPLE_SET>\n<SAMPLE alias=\""+alias+"\" center_name=\""+center+"\">\n")
            xmlfile.write("<TITLE>"+title+"</TITLE>\n<SAMPLE_NAME>\n<TAXON_ID>"+taxonID+"</TAXON_ID>\n<SCIENTIFIC_NAME>"+scientificName+"</SCIENTIFIC_NAME>\n")
            xmlfile.write(" <COMMON_NAME></COMMON_NAME>\n</SAMPLE_NAME>\n<SAMPLE_ATTRIBUTES>\n<SAMPLE_ATTRIBUTE>\n")
            xmlfile.write(" <TAG>geographic location (country and/or sea)</TAG>\n<VALUE>"+geographicLocation+"</VALUE>\n")
            xmlfile.write(" </SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n<TAG>host common name</TAG>\n<VALUE>"+hostCommonName+"</VALUE>\n")
            xmlfile.write(" </SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n<TAG>host health state</TAG>\n<VALUE>"+hostHealthState+"</VALUE>\n")
            xmlfile.write(" </SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n<TAG>isolation source host associated</TAG>\n<VALUE>"+isolationSource+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
            xmlfile.write("<TAG>host sex</TAG>\n<VALUE>"+sex+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
            xmlfile.write("<TAG>host scientific name</TAG>\n<VALUE>"+hostScientificName+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
            xmlfile.write("<TAG>collector name</TAG>\n<VALUE>"+collectorName+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
            xmlfile.write("<TAG>collecting institution</TAG>\n<VALUE>"+collectingInstitution+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n<SAMPLE_ATTRIBUTE>\n")
            xmlfile.write("<TAG>isolate</TAG>\n<VALUE>"+isolate+"</VALUE>\n</SAMPLE_ATTRIBUTE>\n</SAMPLE_ATTRIBUTES>\n</SAMPLE>\n</SAMPLE_SET>\n")
            xmlfile.close()

            os.system("curl -u Webin-50760:Elisaegizia14 -F \"SUBMISSION=@submission.xml\" -F \"SAMPLE=@"+alias+"_sample.xml\" \"https://www.ebi.ac.uk/ena/submit/drop-box/submit/\" >sampleReceipt.txt")
            receiptFile = open("sampleReceipt.txt")
            while True:
                line = receiptFile.readline().rstrip()
                if not line:
                    break
                fields = line.split("\"")
                for item in fields:
                    if "ERS" in item or "SRS" in item or "DRS" in item:
                        receiptFile.close()
                        return item
                        break
            
            #sampleAccession[sampleName] = sampleAccessionNumber
            


    if sampleFound == 0:
        print "Sample",sampleName,"was not found in the provided sample file."
        print "Now exit...."
        exit()

    infile.close()




#*********************************************************
#************** Start main algorithm *********************
#*********************************************************



#Create project if needed
if createProject=="yes":
    projectAccessionNumber = createNewProject(projectInfo)


#Submit reads
toSubmitFile = open(filesToSubmit)
while True:
    sampleName = toSubmitFile.readline().rstrip()
    if not sampleName:
        break

    if not sampleName in experimentAccession:
        experimentAccession[sampleName] = ""

    if not sampleName in runAccession:
        runAccession[sampleName] = ""

    sampleList.append(sampleName)


    #Create samples if needed
    if createSample == "yes":
        sampleAccessionNumber = createNewSample(sampleName)
        time.sleep(3)
    
    #Grab reads information from readsInfo file
    infile = open(readsInfo)
    infile.readline()
    foundRecord = 0
    while True:
        line = infile.readline().rstrip()
        if not line:
            break
        fields = line.split("\t")
        projectAccessionField = fields[0]
        sampleAccessionField = fields[1]

        sampleAlias = fields[2]
        instrument = fields[3]
        insertSize = fields[4]
        librarySource = fields[5]
        librarySelection = fields[6]
        libraryStrategy = fields[7]
        fq1 = fields[8]
        fq2 = fields[9]

        if createProject == "yes":
            projectAccessionField = projectAccessionNumber
        if createSample == "yes":
            sampleAccessionField = sampleAccessionNumber

        sampleAccession[sampleName] = sampleAccessionField

        if sampleAlias == sampleName:
            foundRecord = 1
            os.system("ln -s "+inputFolder+"/"+fq1)
            os.system("ln -s "+inputFolder+"/"+fq2)

            manifest = open(sampleName+"_manifestFile.txt","w",buffering=0)
            manifest.write("INSTRUMENT\t"+instrument+"\nINSERT_SIZE\t"+insertSize+"\nLIBRARY_SOURCE\t"+librarySource+"\nLIBRARY_SELECTION\t"+librarySelection+"\n")
            manifest.write("STUDY\t"+projectAccessionField+"\nSAMPLE\t"+sampleAccessionField+"\nNAME\t"+sampleName+"_hcmv\n")
            manifest.write("LIBRARY_STRATEGY\t"+libraryStrategy+"\nFASTQ\t"+fq1+"\nFASTQ\t"+fq2+"\n")
            manifest.close()
            os.system("~/Software/jre1.8.0_191/bin/java -jar ./webin-cli-1.6.0.jar -context reads -userName Webin-50760 -password Elisaegizia14  -manifest "+sampleName+"_manifestFile.txt -submit >fastqReceipt")
            os.system("rm -f "+fq1+" "+fq2)

            receiptFile = open("fastqReceipt")
            while True:
                line = receiptFile.readline().rstrip()
                if not line:
                    break
                fields = line.split(" ")
                print fields
                for item in fields:
                    if "ERR" in item or "SRR" in item or "DRR" in item: 
                        if not "ERROR" in item:
                            runAccession[sampleName] = item

                    if "ERX" in item or "SRX" in item or "DRX" in item: 
                        experimentAccession[sampleName]=item

                    
            receiptFile.close()
            



    print "Finished"
    print sampleName,sampleAccession[sampleName],runAccession[sampleName],experimentAccession[sampleName]
    sys.stdin.read(1)
    os.system("rm -f *manifestFile.txt *Receipt* *.report")




outfile = open("accessionNumbersTable.txt","w")
outfile.write("SampleAlias\tSampleAccession\tRunAccession\tExperimentAccession\n")

for item in sampleList:
    print item,sampleAccession[item],runAccession[item],experimentAccession[item]
    outfile.write(item+"\t"+sampleAccession[item]+"\t"+runAccession[item]+"\t"+experimentAccession[item]+"\n")

outfile.close()




    
    