import sys

if len(sys.argv)<2:
    print 'error'
    sys.exit()

file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]
file4 = sys.argv[4]
outfile = sys.argv[5]

f1 = open(file1)
f2 = open(file2)
f3 = open(file3)
f4 = open(file4)
ofh = open(outfile,'w')
count1 = 0
count2 = 0
count3 = 0
count4 = 0
for line in f1:
    count1 += 1
#print count1
    
for line in f2:
    count2 += 1

for line in f3:
    count3 += 1

for line in f4:
    count4 += 1

count5 = count1 + count2
print >> ofh , count5

count6 = count3 + count4
print >> ofh , count6
