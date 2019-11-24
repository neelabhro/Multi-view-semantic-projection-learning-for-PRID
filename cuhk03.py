import os
import os.path
import shutil
from os import listdir
from time import sleep


a = listdir('.')
a.remove('cuhk03.py')

names = []
count = {}

for j, i in enumerate(a):
    name = i[:i.find('_', i.find('_')+1)+2]
    try:
        count[name]+=1
    
    except:
        names.append(i)
        count[name] = 1 

# print(names)
# print(count)

for j, i in enumerate(count):
    while(count[i] < 5):
        if(int(i[-1]) == 1):
            shutil.copyfile("./"+names[j], "./"+names[j][:-6]+"0"+str(1+count[i])+".png")
            print("created "+names[j][:-6]+"0"+str(1+count[i])+".png")
        else:
            if(count[i] == 4):
                shutil.copyfile("./"+names[j], "./"+names[j][:-6]+str(6+count[i])+".png")
                print("created "+names[j][:-6]+str(6+count[i])+".png")
            else:
                shutil.copyfile("./"+names[j], "./"+names[j][:-6]+"0"+str(6+count[i])+".png")
                print("created "+names[j][:-6]+"0"+str(6+count[i])+".png")
        count[i]+=1
