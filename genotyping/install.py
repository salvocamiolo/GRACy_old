#!/usr/bin/python
import os
os.system("pwd >installationFolder")
infile = open("installationFolder")
outfile = open("header","w")

instFolder = infile.readline().rstrip()

outfile.write("#!/usr/bin/python\n")
outfile.write("installationFolder = \""+instFolder+"\"\n")
outfile.close()


os.system("sed '1,2d' hcmvGenotypingGUI.py > temp.py")
os.system("cat header  temp.py >installedVersion.py")
os.system("mv installedVersion.py hcmvGenotypingGUI.py")
os.system("chmod +x hcmvGenotypingGUI.py")
os.system("rm header installationFolder temp.py -f")