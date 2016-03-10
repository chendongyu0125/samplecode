#install.packages("RCurl")
#install.packages("XML")
library(RCurl)
library(XML)
library(tictoc)

#Initializing array of info
            


NAME<- vector(mode="character", length=1)
DATE<- vector(mode="character", length=1)
CLASIFICACION<- vector(mode="character", length=1)
UNIV <- vector(mode="character", length=20)
ART <- vector(mode="character", length=1000)
LIB <- vector(mode="character", length=1000)
CAPLIB <- vector(mode="character", length=1000)
NUMBART<-vector(mode="character", length=1)
NUMLIB<-vector(mode="character", length=1)
NUMCAPLIB<-vector(mode="character", length=1)
CERTIFICADO<-vector(mode="character", length=1)
LUGAR<-vector(mode="character", length=1)
DATATOTAL<-cbind(t(NAME),t(DATE),t(CLASIFICACION),t(UNIV),t(ART),t(LIB),t(CAPLIB),t(NUMBART),t(NUMLIB),t(NUMCAPLIB),t(CERTIFICADO),t(LUGAR))
DATAITER<-DATATOTAL

for (grupo in 17200:17000){
  #Obteniendo el número del grupo

  if (grupo<10){
    numero=paste0("0000000000000",grupo)
  }
  if (grupo<100 && grupo>=10){
    numero=paste0("000000000000",grupo)
  }
  if (grupo<1000 && grupo>=100){
    numero=paste0("00000000000",grupo)
  }
  if (grupo<10000 && grupo>=1000){
    numero=paste0("0000000000",grupo)
  }
  if (grupo<100000 && grupo>=10000){
    numero=paste0("000000000",grupo)
  }

  #Getting the url
  url = paste0("http://scienti.colciencias.gov.co:8080/gruplac/jsp/visualiza/visualizagr.jsp?nro=",
             numero)
  
  webpage=htmlTreeParse(url, isHTML=TRUE)


  results = xmlChildren(xmlRoot(webpage))
  a=results$body
  #--------

  #Obtaining the information
  #First of all the name, if it is not empty, we proceed:
  Name<-a[1]$span[1]$text
  if ( (class(Name[1])!="NULL") ){
    Name<-as.character(Name)[6]
    Date<-as.character(a[4]$table[2]$tr[2]$td[1]$text)[6]
    Clasificacion<-as.character(a[4]$table[8]$tr[2]$td[1]$text)[6]
    Certificado<-as.character(a[4]$table[5]$tr[2]$td[1]$text)[6]
    Lugar<-as.character(a[4]$table[3]$tr[2]$td[1]$text)[6]
    Numart<-length(a[18]$table)-1
    Numlib<-length(a[20]$table)-1
    Numcaplib<-length(a[22]$table)-1
    if (length(Clasificacion)==0){
      Clasificacion="NULL"
    }
    #Extracting name of universities
    UNILIST<-a[6]$table
    NUMBuniv=length(UNILIST)
    #Generating empty df of university afiliations
    UNIVITER <- UNIV
    
    for (uu in 2:NUMBuniv){
      univ<-as.character(a[6]$table[uu][1]$tr[1]$td[1]$text)[6]
      UNIVITER[uu-1]=univ
    }
    #Extracting Articles. Generate dataframe of articles
    ARTITER=ART
    #Iterating over articles
    if (Numart>0){
      for (aart in 1:Numart){
        article=as.character(a[18]$table[aart+1]$tr[1]$td[5]$text)[6]
        article2=as.character(a[18]$table[aart+1]$tr[1]$td[3]$text)[6]
        article3=paste0(article,",",article2)
        ARTITER[aart]=article3
      }
    }
    
    #Extracting BOOKS Generate dataframe of articles
    LIBITER=LIB
    #Iterating over articles
    if (Numlib>0){
      for (llib in 1:Numlib){
        libro=as.character(a[20]$table[llib+1]$tr[1]$td[5]$text)[6]
        libro2=as.character(a[20]$table[llib+1]$tr[1]$td[3]$text)[6]
        libro3=paste0(libro,",",libro2)
        LIBITER[llib]=libro3
      }
    }
    
    
    #Extracting BOOKS Generate dataframe of articles
    CAPLIBITER=CAPLIB
    #Iterating over articles
    if (Numcaplib>0){
      for (ccaplib in 1:Numcaplib){
        caplibro=as.character(a[22]$table[ccaplib+1]$tr[1]$td[5]$text)[6]
        caplibro2=as.character(a[22]$table[ccaplib+1]$tr[1]$td[3]$text)[6]
        caplibro3=paste0(caplibro,",",caplibro2)
        CAPLIBITER[ccaplib]=caplibro3
      }
    }

    infogroup=cbind(Name,Date,Clasificacion,t(UNIVITER),t(ARTITER),t(LIBITER),t(CAPLIBITER),Numart,Numlib,Numcaplib,Certificado,Lugar)
    DATATOTAL=rbind(DATATOTAL,infogroup)
   write.csv(DATATOTAL, file="/Users/rodrigoazuero/Dropbox/Didris/COLCIENCIAS/Colcienciasdata12.csv")
    fin=cbind(grupo,Name)
    print(fin)

  }#End if conditionof nombre found
}#End for loop 
