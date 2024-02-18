#!/usr/bin/env python3
#Author: CY Gao

import argparse
import re
from collections import defaultdict


class CutSignalp():

    def __init__(self,fasta,posFile,column,output,module):
        self.fasta=fasta
        self.posFile=posFile
        self.column=column
        self.output=output
        self.module=module

    def makeidex(self):
        with open(self.fasta) as f:
            key_index=defaultdict(list)
            line_num=1
            for line in f:
                line=line.rstrip()
                if not line:
                    line_num+=1
                    continue
                if line.startswith(">"):
                    header=line[1:].strip()
                    header=re.split(r"[ ,;\t]",header)[0]
                    #print(header)
                    line_num+=1
                else:
                    key_index[header].append(line_num)
                    line_num+=1
        return key_index


    def readSeq(self):
        key_index=self.makeidex()
        with open(self.fasta) as f:
            lines=f.readlines()
            for key in key_index.keys():
                seq=[lines[i-1].strip() for i in key_index[key]]
                seq="".join(seq)
                yield key,seq
                

    def readPos(self):
        posdict={}
        with open(self.posFile,'r') as f:
            for line in f:
                if not line.startswith("#"):
                    line=line.strip().split("\t")
                    geneid=line[0]
                    geneid=re.split(r"[ ,;\t]",geneid)[0]
                    #print(len(geneid))
                    try:
                        cs=line[(int(self.column)-1)]
                        cs=re.search(r'(\d+)-\d+',cs)
                        if cs:
                            cs=cs.group(1)
                            cspos=int(cs)
                            posdict[geneid]=cspos
                    except IndexError:
                        pass

        return posdict

    
    def cutSP(self):
        posdict=self.readPos()
        with open(f'{self.output}','w') as w1:
            for geneid,seq in self.readSeq():
                #print(len(geneid))
                if posdict.get(geneid):
                    if self.module=="Protein":
                        w1.write(">"+geneid+"\n"+"M"+seq[posdict[geneid]:]+"\n")
                    else:
                        cspos=posdict[geneid]*3
                        w1.write(">"+geneid+"\n"+"ATG"+seq[cspos:]+"\n")


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--fasta',type=str,help='输入序列的fasta文件')
    parser.add_argument('--prediction',type=str,help='输入信号肽预测结果文件')
    parser.add_argument('--output',type=str,help='结果输出文件')
    parser.add_argument('--column',type=str,default="5",help='信号肽预测剪切位点所在列号')
    parser.add_argument('--module',type=str,default="Protein",help='选择输入的fasta文件是蛋白序列/CDS序列。默认为蛋白序列')
    args=parser.parse_args()
    cutsspp=CutSignalp(args.fasta,args.prediction,args.column,args.output,args.module)
    cutsspp.cutSP()

if __name__=="__main__":
    main()




