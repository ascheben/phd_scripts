import time
import os
import sys
import MySQLdb

db = MySQLdb.connect("localhost",
"curator",
"SNPDatabase",
"mydb" )

def addCultivars(listOfCultivar):
	cursor = db.cursor()
	try:
		cursor.executemany("INSERT INTO Cultivar VALUES (NULL,%s,%s)",listOfCultivar)	
		db.commit()
		print("Added cultivars: " + str(listOfCultivar))
	except:
		db.rollback()
	cursor.close()

def addCultivarWithSpeciesID(Accession,speciesID):
	cursor = db.cursor()
	try:
		cursor.execute("INSERT INTO Cultivar VALUES(NULL,%s,%s)",(Accession,speciesID))
		db.commit()
	except:
		db.rollback()
	cursor.close()

def addAccession(listWithAccession):
	cursor = db.cursor()
	try:
		cursor.executemany("INSERT INTO Accession VALUES(NULL,%s,%s,%s)",listWithAccession)	
		db.commit()
		print("Added accessions: " + str(listWithAccession))
	except:
		print("Exception occured for " + str(listWithAccession))
		db.rollback()
	cursor.close()

def SpeciesIdWithSpecies(speciesName):
	cur = db.cursor(MySQLdb.cursors.DictCursor)
	cur.execute("SELECT idSpecie FROM Specie WHERE SpecieName=%s",(speciesName,))
	result_set = cur.fetchall()
	idSpecies = result_set[0]
	idSpecies = idSpecies["idSpecie"]
	cur.close()
	return idSpecies

def addAccessionSNPRelation(SNPRelations):
	cursor = db.cursor()
	try:
		cursor.executemany("INSERT INTO Accession_has_SNP VALUES (%s,%s,%s)",SNPRelations)
		db.commit()
	except:
		print("Exception occurred while inserting SNPs")
		db.rollback()
	cursor.close()

def getAllAccessions():
	cur = db.cursor(MySQLdb.cursors.DictCursor)
	accessionList = []
	cur.execute("SELECT * FROM Accession")
	result_set = cur.fetchall()
	for row in result_set:
		tupleList = row["idAccession"],row["Accession"],row["Cultivar_idCultivar"]
		accessionList.append(tupleList)
	cur.close()
	return accessionList

def SNPArrayIdWithName(SNPArrayName):
	cur = db.cursor(MySQLdb.cursors.DictCursor)
	cur.execute("SELECT idSNPArray FROM SNPArray WHERE SNPArrayName=%s",(SNPArrayName,))
	result_set = cur.fetchall()
	print(result_set)
	idArray = result_set[0]
	idArray = idArray["idSNPArray"]
	cur.close()
	return idArray

def getAllSNPsWithSNPArrayId(SNPArrayId):
	cur = db.cursor(MySQLdb.cursors.DictCursor)
	SNP = []
	cur.execute("SELECT * FROM SNP WHERE SNPArray_idSNPArray=%s", (SNPArrayId,))
	result_set = cur.fetchall()
	for row in result_set:
		tupleList = row["SNPTableid"],row["SNPid"]
		SNP.append(tupleList)
	cur.close()
	return SNP

def getAllCultivars():
	cur = db.cursor(MySQLdb.cursors.DictCursor)
	allCultivars = []
	cur.execute("SELECT CultivarName FROM Cultivar")
	result_set = cur.fetchall()
	for row in result_set:
		tupleList = row["CultivarName"]
		allCultivars.append(tupleList)
	cur.close()
	return allCultivars

def CultivarIdWithCultivar(cultivar):
	cur = db.cursor(MySQLdb.cursors.DictCursor)
	cur.execute("SELECT idCultivar FROM Cultivar WHERE CultivarName=%s",(cultivar,))
	result_set = cur.fetchall()
	idCultivar = result_set[0]
	idCultivar = idCultivar["idCultivar"]
	cur.close()
	return idCultivar
