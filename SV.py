from pybedtools import BedTool
import re,vcf
from collections import OrderedDict


class SV():
    def __init__(self,file,sample_name=None,gt='hg19',gf=None,sample_flag=False,expansion=500):
        self._gt=gt
        self._sample_flag=sample_flag
        if sample_flag:

            svs,flank1,flank2=self.makeTargets(file,sample_name)
            self._svs=BedTool(svs)
            self._svdict=self.filter(self._svs)
            self._gf=gf
            self._flank1=self.expand(flank1,expansion)
            self._flank2=self.expand(flank2,expansion)
            #print flank1[10069]
            #self._flank1=self.filter(flank1)
            #self._flank2=self.filter(flank2)
            #self._GC=self.gcContent(self._svdict,5)
            self._RLCR=self.rlcrPercent()
        else:
            #svs=BedTool().random(l=100000,n=2400,g='config/'+self._gt+'.genome')
            self._svs=BedTool(file)
            #self._svdict=self.filter(svs)

    def makeTargets(self,file,sample_name):
        if file.endswith('bed'):
            return self.makeTargetsBed(file,sample_name)
        elif file.endswith('vcf.gz'):
            return self.makeTargetsVcf(file,sample_name)
        else:
            print ('incorrect format!')
            exit(1)

    def makeTargetsBed(self,bed_file,sample_name):
        with open(bed_file,'r') as fp:
            nskip=0
            lines=fp.readlines()
            lines=[line for line in lines if line.strip()!='']
            info=re.split('\t| ',lines[0].strip())
            if re.match('\d+',info[1])==None:
                lines=lines[1:]
            nlines=len(lines)
            targets=[]
            targets_flank1=[]
            targets_flank2=[]
            for i in range(nlines):
                info=re.split('\t| ',lines[i].strip())
                #chr20   51908751   51917850    DUP    NA12878   index
                ncol=len(info)
                if ncol<3:
                    print ('Invalid bed format!')
                    exit(1)
                elif ncol<5:
                    info[0]=info[0].replace('chr','').replace('23','X').replace('24','Y')
                    targets.append(tuple(info[:4]+[i-nskip]))
                    if self._sample_flag:
                        targets_flank1.append((info[0],info[1],info[1],info[3],i-nskip))
                        targets_flank2.append((info[0],info[2],info[2],info[3],i-nskip))
                elif info[4]==sample_name:
                    info[0]=info[0].replace('chr','').replace('23','X').replace('24','Y')
                    targets.append(tuple(info+[i-nskip]))
                    if self._sample_flag:
                        targets_flank1.append((info[0],info[1],info[1],info[3],i-nskip))
                        targets_flank2.append((info[0],info[2],info[2],info[3],i-nskip))
                else:
                    nskip+=1
        return targets,targets_flank1,targets_flank2

    def makeTargetsVcf(self,vcf_file,sample_name):
        i=0
        targets=[]
        targets_flank1=[]
        targets_flank2=[]
        for sv in vcf.Reader(filename=vcf_file):
            sample_genotype=sv.genotype(sample_name)
            if sample_genotype.gt_type==None or sample_genotype.gt_type==0:
                continue
            if sv.INFO['SVTYPE']=='DEL' or sv.INFO['SVTYPE']=='DUP':
                targets.append((sv.CHROM,sv.POS,sv.INFO['END'],sv.INFO['SVTYPE'],i))
                targets_flank1.append((sv.CHROM,sv.POS,sv.POS,sv.INFO['SVTYPE'],i))
                targets_flank2.append((sv.CHROM,sv.INFO['END'],sv.INFO['END'],sv.INFO['SVTYPE'],i))
                i+=1
            #chr20   51908751   51917850    DUP    index    
        return targets,targets_flank1,targets_flank2

    def expand(self,svs,nbp):
        return BedTool(svs).slop(b=nbp,g='config/'+self._gt+'.genome')

    def filter(self,svs):
        black_beds=BedTool('config/'+self._gt+'_excluded.bed.gz')
        svs=svs.subtract(black_beds)
        svs_dict=OrderedDict()
        index=-1 if self._sample_flag else 0 
        for sv in svs:
            if svs_dict.get(sv[index])==None:
                svs_dict[sv[index]]=[sv]
            else:
                svs_dict[sv[index]].append(sv)
        return svs_dict

    def gcContent(self,svdict,rnd=1):
        GC={}
        for i in svdict:
            #print svdict[i]
            gc=0
            length=0
            #print BedTool(svdict[i]).nucleotide_content(fi=self._gf)
            for j in BedTool(svdict[i]).nucleotide_content(fi=self._gf):
                gc+=int(j[-5])+int(j[-4])
                length+=int(j[-1])
            #print gc,length
            GC[i]=rnd*round(gc*100/length/rnd)
        return GC

    def rlcrPercent(self):
        rlcr={} #repetitive and low complexity region
        for target in self._svs:
            i=target[-1]
            if self._svdict.get(i)==None:
                rlcr[i]=1
            elif len(self._svdict[i])==1:
                rlcr[i]=0
            else:
                nbps=0
                for target in self._svdict[i]:
                    nbps+=int(target[2])-int(target[1])
                rlcr[i]=1-nbps*1.0/(int(self._svs[int(i)][2])-int(self._svs[int(i)][1]))
        return rlcr





