import sys
from flaskext.mysql import MySQL
sys.path.insert(0, "../")
from app import app

mysql = MySQL(app)

db = mysql.connect()

cursor = db.cursor()

manifest = sys.argv[1]
sp = sys.argv[2]
sp_list = sp.split(",")
mylines = []
with open(manifest) as f:
    for line in f:
        row = line.strip().split(',')
        mylines.append(row) 

count = 0
for SNP in mylines:
    try:
        SNPid = str(SNP[0])
        #Chromosome = str(SNP[1])
        #Coordination = str(SNP[2])
        Reference = str(SNP[1])
        Strand = str(SNP[2])
        SourceSeq = str(SNP[3])
        SNPArray_idSNPArray = int(SNP[4])
        sp_count = 5
        for i in range(len(sp_list)):
            SpeciesID = sp_list[i]
            Chromosome = str(SNP[sp_count])
            Coordinate = str(SNP[sp_count+1])
            sp_count = sp_count + 2 
        # KEN: NEED to check whether this (Chromosome,SNPid) is already in the database, if yes, skip the insert.
            cursor.execute('''INSERT into SNP (SNPid,Reference, Strand, SourceSeq, SNPArray_idSNPArray) values (%s, %s, %s, %s, %s, %s, %s)''',(SNPid, Reference, Strand, SourceSeq, SNPArray_idSNPArray))
            cursor.execute('''INSERT into SNP_Coordinates (Chromosome, Coordinate, Specie_idSpecie) values (%s, %s, %s)''',(Chromosome, Coordinate, SpeciesID))
        db.commit()
        count = count + 1
    except:
        db.rollback()
        print("There was an error")
db.close()
print(str(count) + " SNPs added to database. Exiting!")
