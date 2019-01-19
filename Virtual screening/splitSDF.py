# Useful script for splitting a concatenated SDF (for eg: those from PubChem)
# Courtesy of someone on StackOverflow

f= "compound" # name of SDF file from database
split_number= 1000 # required #molecules per chopped file
number_of_sdfs=0
i=0
j=0
f2=open(f+'_'+str(j)+'.sdf','w')
for line in open(f+'.sdf'):
	f2.write(line)
	if line[:4] == "$$$$":
		i+=1
	if i > number_of_sdfs:
		number_of_sdfs += split_number
		f2.close()
		j+=1
		f2=open(f+'_'+str(j)+'.sdf','w')
print(i)
