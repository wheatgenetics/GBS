fname='/Users/mlucas/Desktop/GBS_Links.txt'
outFile='/Users/mlucas/Desktop/GBS_wget_commands.txt'
wgetList=[]
gbsPart=''
httpPart=''

with open(fname) as f:
    content = f.readlines()
print('')
linecount=0
for l in content:
    if l[0:3]=='GBS':
        gbsPart=l
    elif l[0:4]=='http':
        httpPart=l
        wgetCommand ='wget -c ' + httpPart[0:-1] + ' -O ' + gbsPart[0:-2]+'\n'
        wgetList+=wgetCommand
        linecount += 1
        print(linecount,wgetCommand)
    else:
        wgetCommand=''

with open(outFile,'w') as o:
    for l in wgetList:
        o.write(l)

