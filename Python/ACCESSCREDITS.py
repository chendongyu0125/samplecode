import requests
import unicodedata
import csv
import numpy
import os
import codecs,cStringIO

#Defining for unicodewriter

class UTF8Recoder:
	def __init__(self, f, encoding):
		self.reader = codecs.getreader(encoding)(f)
	def __iter__(self):
		return self
	def next(self):
		return self.reader.next().encode("utf-8")

class UnicodeReader:
	def __init__(self, f, dialect=csv.excel, encoding="utf-8-sig", **kwds):
		f = UTF8Recoder(f, encoding)
		self.reader = csv.reader(f, dialect=dialect, **kwds)
	def next(self):
		'''next() -> unicode
		This function reads and returns the next line as a Unicode string.
		'''
		row = self.reader.next()
		return [unicode(s, "utf-8") for s in row]
	def __iter__(self):
		return self

class UnicodeWriter:
	def __init__(self, f, dialect=csv.excel, encoding="utf-8-sig", **kwds):
		self.queue = cStringIO.StringIO()
		self.writer = csv.writer(self.queue, dialect=dialect, **kwds)
		self.stream = f
		self.encoder = codecs.getincrementalencoder(encoding)()
	def writerow(self, row):
		'''writerow(unicode) -> None
		This function takes a Unicode string and encodes it to the output.
		'''
		self.writer.writerow([s.encode("utf-8") for s in row])
		data = self.queue.getvalue()
		data = data.decode("utf-8")
		data = self.encoder.encode(data)
		self.stream.write(data)
		self.queue.truncate(0)
	def writerows(self, rows):
		for row in rows:
			self.writerow(row)



os.chdir("/Users/rodrigoazuero/Dropbox/Didris/SPADIES")
url="http://spadies.mineducacion.gov.co/spadies/spadies/consultas"

#Guardamos únicamente las variables de ICFES-examen de estado para realizarla a través de 
#años-semestres. 

#63->Periodo
#64->Número de semestre cursado
#37->Periodo de ingreso

name = 'ICETEXACCES.csv'


#This part of the file takes the number of students with 
#'TIPO DE CREDITO ICETEX RECIBIDO'. We will also use a different version below. 

with open('universities.csv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
	for institut in spamreader:
		print(institut)
		if int(institut[0]) == 1101:
			initial = 0
			types = 'w'
		else:
			initial = 1
			types = 'a'
		fil='!2--,'+str(institut[0])
		modo='0'
		stri='c=4&diferenciados=73,63,&fil='+fil+'!64--,1&modo=0'
		print(stri)
		r = requests.post(url, data=stri)
		sol = requests.post(url, data=stri)
		JSONDATA=sol.json()
		if JSONDATA[0:14] != 'La consulta no':
			print('here')
			rows = len(JSONDATA[3]['valores'])
			print('here2')
			univ = int(institut[0])*numpy.ones((rows-initial,1))
			print('here3')
			A = numpy.array(JSONDATA[3]['valores'][initial:rows])
			A = numpy.hstack((univ,A))
			print(A)
			with open(name, types) as fp:
				a = UnicodeWriter(fp,quoting=csv.QUOTE_ALL)
				a.writerows(A)





#This part of the file uses a different form of getting the students who
#have a credit with icetex. We do it by asking what is the number of credits 
#from icetex they have. 


name = 'ICETEXACCES2.csv'
with open('universities.csv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=' ', quotechar='|')
	for institut in spamreader:
		print(institut)
		if int(institut[0]) == 1101:
			initial = 0
			types = 'w'
		else:
			initial = 1
			types = 'a'
		fil='!2--,'+str(institut[0])
		modo='0'
		stri='c=4&diferenciados=30,63,&fil='+fil+'!64--,1&modo=0'
		print(stri)
		r = requests.post(url, data=stri)
		sol = requests.post(url, data=stri)
		JSONDATA=sol.json()
		if JSONDATA[0:14] != 'La consulta no':
			print('here')
			rows = len(JSONDATA[3]['valores'])
			print('here2')
			univ = int(institut[0])*numpy.ones((rows-initial,1))
			print('here3')
			A = numpy.array(JSONDATA[3]['valores'][initial:rows])
			A = numpy.hstack((univ,A))
			print(A)
			with open(name, types) as fp:
				a = UnicodeWriter(fp,quoting=csv.QUOTE_ALL)
				a.writerows(A)
