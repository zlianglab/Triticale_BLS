import sys
import os

mdict = {}

with open(sys.argv[1],'r') as fh:
    fh.readline()
    for line in fh:
        new = line.strip().split('\t')
        mdict[new[0]] = new[1]

mfold = sys.argv[2]

mfiles = os.listdir(mfold)

for mfile in mfiles:
    a = '{}/{}/tomtom.tsv'.format(mfold, mfile)
    with open(a,'r') as fh:
        fh.readline()
        for line in fh:
            if len(line.strip()) > 1:
                if line.startswith('#'):continue
                new = line.strip().split('\t')
                target = new[1]
                if target in mdict:
                    name = mdict[target]
                else:
                    name = target
                print (mfile, new[0], target, name)
