import sys
import os

args = sys.argv
print(args)
num = int((len(args)-2)/3)

os.system("cd ../")

for i in range(2,num+2):
    dir = "/".join(args[2*num+i].split("/")[:-1])
    command = "salmon quant -i {0} -l A -1 {1} -2 {2} -p 8 -o {3}".format(args[1], args[i], args[num+i], dir)
    os.system(command)
