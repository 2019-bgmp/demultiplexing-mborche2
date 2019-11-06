#Read in files to be demultiplexed
import argparse
def get_args():
	parser = argparse.ArgumentParser(description="DemultiplexingProgram")
	parser.add_argument("-R1", "--ForwardRead", required=True)
	parser.add_argument("-R2", "--ForwardBarcode", required=True)
	parser.add_argument("-R3", "--ReverseBarcode", required=True)
	parser.add_argument("-R4", "--ReverseSequence", required=True)
	return parser.parse_args()

args = get_args()
R1_File = args.ForwardRead
R2_File = args.ForwardBarcode
R3_File = args.ReverseBarcode
R4_File = args.ReverseSequence

import gzip

#Store all of the barcodes
barcodelist = ["GTAGCGTA", "CGATCGAT", "GATCAAGG", "AACAGCGA", "TAGCCATG", "CGGTAATC", "CTCTGGAT", "TACCGGAT", "CTAGCTCA", "CACTTCAC", "GCTACTCT", "ACGATCAG", "TATGGCAC", "TGTTCCGT", "GTCCTAAG", "TCGACAAG", "TCTTCGAC", "ATCATGCG", "ATCGTGGT", "TCGAGAGT", "TCGGATTC", "GATCTTGC", "AGAGTCCA", "AGGATAGC"]
#Set count variables and initialize dictionaries
PairedCount=0
IndexHoppedCount=0
UnknownLowQCount=0
IndexPairsDictionary={}
IndexHoppedPairsDictionary={}

#R1_File = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R1_001.fastq.gz"
#R2_File = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R2_001.fastq.gz"
#R3_File = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R3_001.fastq.gz"
#R4_File = "/projects/bgmp/shared/2017_sequencing/1294_S1_L008_R4_001.fastq.gz"
#R1_File = "testR1.fastq.gz"
#R2_File = "testR2.fastq.gz"
#R3_File = "testR3.fastq.gz"
#R4_File = "testR4.fastq.gz"

with gzip.open(R1_File, "rt") as f1, gzip.open(R2_File, "rt") as f2, gzip.open(R3_File, "rt") as f3, gzip.open(R4_File, "rt") as f4:
	#Forward Sequence
	while True:
		R1L1 = f1.readline().strip()
		if R1L1 == "":
			break
		R1L2 = f1.readline().strip()
		R1L3 = f1.readline().strip()
		R1L4 = f1.readline().strip()
		#Forward Barcode
		R2L1 = f2.readline().strip()
		R2L2 = f2.readline().strip()
		R2L3 = f2.readline().strip()
		R2L4 = f2.readline().strip()
		#Reverse Barcode
		R3L1 = f3.readline().strip()
		R3L2 = f3.readline().strip()
		R3L3 = f3.readline().strip()
		R3L4 = f3.readline().strip()
		#Reverse Sequence
		R4L1 = f4.readline().strip()
		R4L2 = f4.readline().strip()
		R4L3 = f4.readline().strip()
		R4L4 = f4.readline().strip()
		#Change reverse barcode to reverse compliment
		Barcode_Reverse_Complement = R3L2[::-1]
		New_Barcode=[]
		for index in range(len(Barcode_Reverse_Complement)):
			New_Barcode.append(0)
		for index in range(len(Barcode_Reverse_Complement)):
			if Barcode_Reverse_Complement[index] == "A":
				New_Barcode[index] = "T"
			elif Barcode_Reverse_Complement[index] == "T":
				New_Barcode[index] = "A"
			elif Barcode_Reverse_Complement[index] == "C":
				New_Barcode[index] = "G"
			elif Barcode_Reverse_Complement[index] == "G":
				New_Barcode[index] = "C"
		Actual_Reverse_Complement = "".join(str(i) for i in New_Barcode)
		scores = []
		for letter in R2L4:
			scores.append(0)
		if R2L2 in barcodelist:
			for index in range(len(R2L4)):
				scores[index] = ord(str(R2L4[index]))-33
				if min(scores) >= 30:
					for index in range(len(R2L4)):
						scores[index] = ord(str(R3L4[index]))-33
						if min(scores) >= 30:
							if R2L2 == Actual_Reverse_Complement:
								filenme = R2L2+"-"+R3L2
								filenme2 = R2L2+"-"+R3L2+"_R"+".fastq"
								filenme3 = R2L2+"-"+R3L2+".fastq"
								with open(filenme3, "a+") as fil:
									fil.write(R1L1)
									fil.write("\n")
									fil.write(R1L2)
									fil.write("\n")
									fil.write(R1L3)
									fil.write("\n")
									fil.write(R1L4)
									fil.write("\n")
									fil.close()
								with open(filenme2, "a+") as fil:
									fil.write(R4L1)
									fil.write("\n")
									fil.write(R4L2)
									fil.write("\n")
									fil.write(R4L3)
									fil.write("\n")
									fil.write(R4L4)
									fil.write("\n")
									fil.close()
								if filenme not in IndexPairsDictionary:
									IndexPairsDictionary[filenme] = 1
								else:
									IndexPairsDictionary[filenme] += 1
								PairedCount += 1
							elif Actual_Reverse_Complement in barcodelist:
								filenme = R2L2+"-"+R3L2
								with open("Index_Hopped.fastq", "a+") as fil:
									fil.write(R1L1)
									fil.write("\n")
									fil.write(R1L2)
									fil.write("\n")
									fil.write(R1L3)
									fil.write("\n")
									fil.write(R1L4)
									fil.write("\n")
									fil.close()
								with open("Index_Hopped_R.fastq", "a+") as filh:
									filh.write(R1L1)
									filh.write("\n")
									filh.write(R1L2)
									filh.write("\n")
									filh.write(R1L3)
									filh.write("\n")
									filh.write(R1L4)
									filh.write("\n")
									filh.close()
								if filenme not in IndexHoppedPairsDictionary:
									IndexHoppedPairsDictionary[filenme] = 1
								else:
									IndexHoppedPairsDictionary[filenme] += 1
								IndexHoppedCount += 1
							else:
								with open("Unknown_And_Low_Quality.fastq", "a+") as fil:
									fil.write(R1L1)
									fil.write("\n")
									fil.write(R1L2)
									fil.write("\n")
									fil.write(R1L3)
									fil.write("\n")
									fil.write(R1L4)
									fil.write("\n")
									fil.close()
								with open("Unknown_And_Low_Quality_R.fastq", "a+") as fil:
									fil.write(R1L1)
									fil.write("\n")
									fil.write(R1L2)
									fil.write("\n")
									fil.write(R1L3)
									fil.write("\n")
									fil.write(R1L4)
									fil.write("\n")
									fil.close()
								UnknownLowQCount += 1
						else:
							with open("Unknown_And_Low_Quality.fastq", "a+") as fil:
								fil.write(R1L1)
								fil.write("\n")
								fil.write(R1L2)
								fil.write("\n")
								fil.write(R1L3)
								fil.write("\n")
								fil.write(R1L4)
								fil.write("\n")
								fil.close()
							with open("Unknown_And_Low_Quality_R.fastq", "a+") as fil:
								fil.write(R1L1)
								fil.write("\n")
								fil.write(R1L2)
								fil.write("\n")
								fil.write(R1L3)
								fil.write("\n")
								fil.write(R1L4)
								fil.write("\n")
								fil.close()
							UnknownLowQCount += 1
				else:
					with open("Unknown_And_Low_Quality.fastq", "a+") as fi:
						fi.write(R1L1)
						fi.write("\n")
						fi.write(R1L2)
						fi.write("\n")
						fi.write(R1L3)
						fi.write("\n")
						fi.write(R1L4)
						fi.write("\n")
						fi.close()
					with open("Unknown_And_Low_Quality_R.fastq", "a+") as fih:
						fih.write(R1L1)
						fih.write("\n")
						fih.write(R1L2)
						fih.write("\n")
						fih.write(R1L3)
						fih.write("\n")
						fih.write(R1L4)
						fih.write("\n")
						fih.close()
					UnknownLowQCount += 1

#Print Count variables and summary information
print(IndexPairsDictionary)
with open("Summary_Demultiplex.txt", "w+") as summ:
	summ.write("Index Matched Counts:")
	matched=str(PairedCount)
	summ.write(matched)
	for key,value in IndexPairsDictionary.items():
		info = str(key)+":"+str(value)
		summ.write("\n")
		summ.write(info)
	summ.write("\n")
	summ.write("Index Hopped Counts")
	hopped=str(IndexHoppedCount)
	summ.write(hopped)
	for key,value in IndexHoppedPairsDictionary.items():
		info = str(key)+":"+str(value)
		summ.write("\n")
		summ.write(info)
	summ.write("\n")
	summ.write("Unknown and Low Quality Count:")
	unknwn = str(UnknownLowQCount)
	summ.write(unknwn)
	countdictionary ={}
	for item in barcodelist:
		countdictionary[item] = 0
	for key in IndexHoppedPairsDictionary:
		for item in barcodelist:
			if str(item) in str(key):
				value = IndexHoppedPairsDictionary[key]
				countdictionary[item] += value
	for key in IndexPairsDictionary:
		for item in barcodelist:
			if str(item) in str(key):
				value = IndexPairsDictionary[key]
				countdictionary[item] += value
	summ.write("\n")
	summ.write("Percent of total reads per barcode:")
	for key,value in countdictionary.items():
		info = str(key)+":"+str(value)
		summ.write("\n")
		summ.write(info)
	summ.close()
