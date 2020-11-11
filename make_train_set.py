from pybedtools import BedTool
import vcf

def make_trains(sv,tsv,tag,f):
	count=0
	for tar in tsv:
		count+=1
		#print count
		flag=False
		for tar1 in sv:
			if overlap(tar,tar1):
				flag=True
		if flag:
			f.write(tar[0]+'\t'+tar[1]+'\t'+tar[2]+'\t'+tar[3]+'\t'+tag+'\t'+'1\n')
		else:
			f.write(tar[0]+'\t'+tar[1]+'\t'+tar[2]+'\t'+tar[3]+'\t'+tag+'\t'+'0\n')


def overlap(t1,t2):
	s1=int(t1[1])
	e1=int(t1[2])
	s2=int(t2[1])
	e2=int(t2[2])
	if t1[0]!=t2[0] or e1<s2 or e2<s1:
		return False
	else:
		olap=min(e1,e2)-max(s1,s2)
		ratio=olap*1.0/(e1-s1)
		if ratio>0.1:
			return True
		else:
			return False

def checkChrom(chr):
	chrchrom=['chr'+str(i) for i in range(1,23)]+['chrX','chrY']
	chrom=[str(i) for i in range(1,23)]+['X','Y']
	if chr in chrchrom or chr in chrom:
		return True
	else:
		return False



def makeTargetsVcf(vcf_file,sample_name,flag=True):
    i=0
    targets=[]
    #targets_flank1=[]
    #targets_flank2=[]
    for sv in vcf.Reader(filename=vcf_file):
        sample_genotype=sv.genotype(sample_name)
        if flag:
            if sample_genotype.gt_type==None or sample_genotype.gt_type==0:
                continue
        if sv.INFO['SVTYPE']=='DEL' or sv.INFO['SVTYPE']=='DUP':
        	if checkChrom(sv.CHROM):
				targets.append((sv.CHROM,sv.POS,sv.INFO['END'],sv.INFO['SVTYPE'],i))
            	#targets_flank1.append((sv.CHROM,sv.POS,sv.POS,sv.INFO['SVTYPE'],i))
            	#targets_flank2.append((sv.CHROM,sv.INFO['END'],sv.INFO['END'],sv.INFO['SVTYPE'],i))
            	i+=1

    return BedTool(targets)




gsd='/data/xwu/gsd/ALL.wgs.mergedSV.v8.20130502.svs.genotypes.vcf.gz'
lump='/home/xwu/opsv/lumpy/NA12878high_lumpy_sort.vcf.gz'
dell='/home/xwu/opsv/delly/NA12878high_delly.vcf.gz'

print ('sv')
sv=makeTargetsVcf(gsd,'NA12878')
print (len(sv))
print ('lsv')
lsv=makeTargetsVcf(lump,'NA12878',False)
print (len(lsv))
print ('dsv')
dsv=makeTargetsVcf(dell,'NA12878')
print (len(dsv))
f=open('train_set.bed','w')
f.write('chromsome\tstart\tend\tgenotype\tsample\tcategory\n')
make_trains(sv,lsv,'NA12878',f)
make_trains(sv,dsv,'NA12878',f)
f.close()


