
f=open('train_set.bed','r')
new_f=open('new_train_set.bed','w')
count={'DEL':0,'DUP':0}
l=f.readline()
while l!='':
	ls=l.strip().split('\t')
	genotype=ls[3]
	category=ls[5]
	if genotype=='genotype':
		new_f.write(l)
	else:
		if category=='1':
			count[genotype]+=1
			new_f.write(l)
		else:
			if count[genotype]>0:
				new_f.write(l)
				count[genotype]-=1
	l=f.readline()
f.close()
new_f.close()

