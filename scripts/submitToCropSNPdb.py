import time
import os
import sys
from mysqlqueries import *
import argparse
from argparse import RawTextHelpFormatter
import os.path

try:
    import configparser 
except ImportError:
    import ConfigParser  # Python 2.7

__version__ = "0.1"
__author__ = "Armin Scheben"

def parse_user_args():
	"""
	This function parses the arguments provided by the user
	:return: a dictionary having a key for each arguments
	:rtype: dict
	"""
	parser = argparse.ArgumentParser(
	description='Welcome to CropSNPdb Submitter %s.\n',
	usage='python submitToCropSNPdb.py -g [GENOTYPE_FILE] -c [CULTIVAR_FILE] -a [ARRAY_NAME]',
	formatter_class=RawTextHelpFormatter, add_help=False)

	mandatory = parser.add_argument_group('mandatory arguments')

	mandatory.add_argument(
	    '-g', '--genotype', dest='genotype', required=True, metavar='GENOTYPE TABLE', help='A CSV table of accessions and genotype calls ')

	mandatory.add_argument(
	    '-c', '--cultivar', dest='cultivar', required=True, 
	    help='A CSV table of cultivars and corresponding accessions and reference IDs')

	mandatory.add_argument(
	'-a', '--array', dest='array', required=True, metavar='PATH', 
	help='Brassica60k || Wheat90k')

	optional = parser.add_argument_group('optional arguments')

	optional.add_argument('-h', '--help', action="help", help="Show this help message and exit")

	return vars(parser.parse_args())


#Takes specie, path to file with SNP data, path to file with cultivar data and the SNParray used
#Opens the files and looks for the cultivars and puts the cultivar in a dic with the accession
#If a cultivar isn't know yet adds it to the database
#Afterwards loops through the snp file and puts SNPs into the database 
#Submits to database after 100000 SNPs
def SubmitFiles(SNPFilePath,cultivarFilePath,SNPArrayName):
	accessionWithCultivar = {}
	cultivarFile = open(cultivarFilePath,'r')
	countGenotypes = 0
	#specieId = SpeciesIdWithSpecies(Species)
	allCultivars = []
	accessionToAdd = []
	cultivarsToAdd = []
	allCultivars = getAllCultivars()
	for line in cultivarFile:
		if "Sample_ID" not in line:
			line = line.replace('\n','')
			line = line.replace('\r','')
			splitLine = line.split(',')
			accession = splitLine[0]
			cultivar = splitLine[1]
			refid = splitLine[2]
			speciesID = splitLine[3]

			if cultivar not in allCultivars and cultivar != "":
				addCultivarWithSpeciesID(cultivar,speciesID)
				allCultivars.append(cultivar)
			elif cultivar not in allCultivars and cultivar == "":
				cultivar = accession
				allCultivars.append(cultivar)
				addCultivarWithSpeciesID(cultivar,speciesID)

			cultivarId = CultivarIdWithCultivar(cultivar)
			accessionToAdd.append((accession,cultivarId,refid))
			accessionWithCultivar[accession] = cultivar
	#There is no check to kick out duplicate accessions with same cultivar and RefID
	addAccession(accessionToAdd)
	cultivarFile.close()

	SNPFile = open(SNPFilePath,'r')
	SNPRelations = []
	i = 0
	beginTime = time.time()
	SNPArrayId = SNPArrayIdWithName(SNPArrayName)
	allAccessions = getAllAccessions()
	allSNP = list(getAllSNPsWithSNPArrayId(SNPArrayId)) 
	for line in SNPFile:
		line = line.replace('\n','')
		line = line.replace('\r','')
		splitLine = line.split(',')
		if splitLine[0] == 'SNP_ID':
			accessions = splitLine[1:]
		else:
			SNPTableId = ""
			SNPId = splitLine[0]
			genotypes = splitLine[1:]
			SNPFound = False
			j = 0
			while(not SNPFound and j < len(allSNP)):
				databaseSNP = allSNP[j]
				if(SNPId == databaseSNP[1]):
					SNPTableId = databaseSNP[0]
					del allSNP[j]
					SNPFound = True
				j = j + 1
			if(SNPTableId!=""):
				for accession,genotype in zip(accessions,genotypes):
					accessionId = ""
					accessionFound = False
					j = 0
					while(not accessionFound and j < len(allAccessions)):
						databaseAccession = allAccessions[j]
						if(accession == databaseAccession[1]):
							accessionId = databaseAccession[0]
							accessionFound = True
						j = j + 1
					SNPRelations.append((accessionId,SNPTableId,genotype))	
					if(i % 100000 is 0 and i != 0):
						print(i)
						print(time.time()-beginTime)
						addAccessionSNPRelation(SNPRelations)
						SNPRelations = []
					i = i + 1
	addAccessionSNPRelation(SNPRelations)
	SNPFile.close()	
def main():
	params = parse_user_args()
	genofile = params['genotype']
	cultivarfile = params['cultivar']
	arrayname = params['array']
	#There is no check to verify array name exists
	#Duplicate accessions are not handled, just added redundantly
	if os.path.isfile(genofile) and os.path.isfile(cultivarfile): 
		SubmitFiles(genofile, cultivarfile, arrayname)
	else:
		print("Error loading input files. Check file pathes are correct.")
if __name__ == "__main__":
	main()
