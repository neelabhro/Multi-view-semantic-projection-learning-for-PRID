import os
import os.path
import shutil
from os import listdir
from time import sleep


a = listdir('.')
for i in a:
    if(not os.path.exists('./'+i[-6:-4])):
        os.mkdir('./'+i[-6:-4])
    shutil.copyfile("./"+i, './'+i[-6:-4]+'/'+i)
