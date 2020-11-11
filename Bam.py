import pysam
import numpy as np
from Snv import median

class Bam():
	def __init__(self, bam_file):
		self._fn=bam_file
		self._fp=pysam.AlignmentFile(bam_file,'rb')
		if not self._fp.check_index():
			print (bam_file,'has no index !')
			exit(1)
		self._total_reads=self._fp.mapped
		self._chr_prefix='chr' if 'chr' in ''.join(self._fp.references) else ''
		self._chrs=[self._chr_prefix+str(i) for i in range(1,23)]+[self._chr_prefix+'X',self._chr_prefix+'Y']
		lengths=[self._fp.get_reference_length(i) for i in self._chrs]
		self._lengths=dict(zip(self._chrs,lengths))
		self._gt=self.getVersion()

	# def chr_length(chr_length_file,bed_flag=False):
	# 	with open(chr_length_file,'w') as f:
	# 		for c in self._chrs:
	# 			if bed_flag:
	# 				s=c+'\t0\t'+str(self._lengths[c])+'\n'
	# 			else:
	# 				s=c+'\t'+str(self._lengths[c])+'\n'
	# 			f.write(s)


	def rpkm(self,sv):
		rpkm={}
		for i in sv._svdict:
			nreads=0
			npos=0
			for target in sv._svdict[i]:
				reads=self._fp.fetch(self._chr_prefix+target[0],int(target[1]),int(target[2]))
				npos+=int(target[2])-int(target[1])
				for read in reads:
					if self.checkRead(read)==2:
						nreads+=1
			if npos==0:
				rpkm[i]=0
			else:
				rpkm[i]=10**9*nreads*1.0/npos/self._total_reads
		return rpkm

	def rpkm2(self,sv):
		rpkm={}
		for target in sv._svs:
			nreads=0
			npos=0
			reads=self._fp.fetch(self._chr_prefix+target[0],int(target[1]),int(target[2]))
			npos+=int(target[2])-int(target[1])
			for read in reads:
				if self.checkRead(read)==2:
					nreads+=1
			if npos==0:
				rpkm[target[-1]]=0
			else:
				rpkm[target[-1]]=10**9*nreads*1.0/npos/self._total_reads
		return rpkm


	def coverage(self,sv,length=0):
		coverage={}
		for i in sv._svdict:
			cov=[]
			for target in sv._svdict[i]:
				deps=pysam.depth("-a", "-Q" "10", "-r", '{}:{}-{}'.format(self._chr_prefix+target[0],target[1],target[2]), "-l", str(length), self._fn).strip().split('\n')
				for dep in deps:
					dep=dep.strip().split('\t')
					if len(dep)<3:
						continue
					cov.append(int(dep[2]))
				coverage[i]=median(cov,1)
		return coverage


	def SL(self,sv,count=10000):
		isize={}
		rlength={}
		for i in sv._svdict:
			ct=0
			size=[]
			length=[]
			for target in sv._svdict[i]:
				reads=self._fp.fetch(self._chr_prefix+target[0],int(target[1]),int(target[2]))
				for read in reads:
					if self.checkRead(read)==2:
						ct+=1
						size.append(abs(read.template_length))
						length.append(read.query_alignment_length)
						if ct>count:
							break
				if ct>count:
					break
			isize[i]=median(size)+3*std(size)
			rlength[i]=median(length)
		return (isize,rlength)

	def SL2(self,sv,count=10000):
		isize={}
		rlength={}
		for target in sv._svs:
			ct=0
			size=[]
			length=[]
			reads=self._fp.fetch(self._chr_prefix+target[0],int(target[1]),int(target[2]))
			for read in reads:
				if self.checkRead(read)==2:
					ct+=1
					size.append(abs(read.template_length))
					length.append(read.query_alignment_length)
					if ct>count:
						break
			isize[target[-1]]=median(size)+3*std(size)
			rlength[target[-1]]=median(length)
		return (isize,rlength)

	def checkRead(self,read,mq=10):
		if read.is_duplicate or read.reference_id!=read.next_reference_id or read.is_qcfail or read.is_reverse==read.mate_is_reverse or read.is_unmapped or read.mate_is_unmapped:
			return 0
		elif read.is_secondary or read.is_supplementary or read.mapping_quality<mq or not read.is_paired or not read.is_proper_pair:
			return 1
		else:
			return 2


	def discordant(self,sv,thr,expan=500):
		discordant={}
		for i in sv._svdict:
			ndis=0
			ncon=0
			nsp=0
			flank=sv._flank1[i]+sv._flank2[i] if sv._flank1.get(i)!=None and sv._flank2.get(i)!=None else []
			for target in flank:
				reads=self._fp.fetch(self._chr_prefix+target[0],int(target[1]),int(target[2]))
				for read in reads:
					if self.checkRead(read)==0:
						continue
					elif self.checkRead(read)==2:
						ncon+=1
					else:
						if read.is_paired:
							if (read.template_length>thr and (abs(read.reference_start-int(sv._svs[int(i)][1]))<expan and abs(read.next_reference_start-int(sv._svs[int(i)][2]))<expan 
							or abs(read.reference_start-int(sv._svs[int(i)][2]))<expan and abs(read.next_reference_start-int(sv._svs[int(i)][1]))<expan)):
								ndis+=1
						if read.is_secondary or read.is_supplementary:
							align = read.get_tag("SA").split(',')
							if (align[0]==read.reference_name and (abs(read.reference_start-int(sv._svs[int(i)][1]))<expan and abs(int(align[1])-int(sv._svs[int(i)][2]))<expan 
							or abs(read.reference_start-int(sv._svs[int(i)][2]))<expan and abs(int(align[1])-int(sv._svs[int(i)][1]))<expan)):
								nsp+=1
			discordant[i]=(ndis,nsp,ncon)
		return discordant

	def discordant2(self,sv,thr,expan=500):
		discordant={}
		for target in sv._svs:
			i=target[-1]
			ndis=0
			ncon=0
			nsp=0
			flank=[sv._flank1[int(i)],sv._flank2[int(i)]] 
			for tar in flank:
				reads=self._fp.fetch(self._chr_prefix+tar[0],int(tar[1]),int(tar[2]))
				for read in reads:
					if self.checkRead(read)==0:
						continue
					elif self.checkRead(read)==2:
						ncon+=1
					else:
						if read.is_paired:
							if (read.template_length>thr and (abs(read.reference_start-int(sv._svs[int(i)][1]))<expan and abs(read.next_reference_start-int(sv._svs[int(i)][2]))<expan 
							or abs(read.reference_start-int(sv._svs[int(i)][2]))<expan and abs(read.next_reference_start-int(sv._svs[int(i)][1]))<expan)):
								ndis+=1
						if read.is_secondary or read.is_supplementary:
							align = read.get_tag("SA").split(',')
							if (align[0]==read.reference_name and (abs(read.reference_start-int(sv._svs[int(i)][1]))<expan and abs(int(align[1])-int(sv._svs[int(i)][2]))<expan 
							or abs(read.reference_start-int(sv._svs[int(i)][2]))<expan and abs(int(align[1])-int(sv._svs[int(i)][1]))<expan)):
								nsp+=1
			discordant[i]=(ndis,nsp,ncon)
		return discordant

	def getVersion(self):
		chr1_length=self._fp.lengths[0]
		if chr1_length==249250621:
			return 'hg19'
		elif chr1_length==248956422:
			return 'hg38'
		elif chr1_length==247249719:
			return 'hg18'
		else:
			return 'Unknown genome version !'
			exit(1)

def std(l):
	if len(l)==0:
		return float('nan')
	else:
		return np.std(l)


