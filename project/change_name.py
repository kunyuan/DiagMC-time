#!/usr/bin/python
import os
old="0.50"
new="0.70"
for file in os.listdir("."):
    if file.startswith(old):
        newfile=new+file[4:] 
        os.system("mv "+file+" "+newfile)
