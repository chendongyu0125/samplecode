from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from bs4 import BeautifulSoup
import re
import pandas as pd
import os
import urllib2
import requests
import lxml.html as lh
import csv
import time
#chrome driver should be on path. MAC: /usr/local/bin
driver = webdriver.Chrome()

#Chose poverty line in 2015: 
InputPov=2.6

url = "http://iresearch.worldbank.org/PovcalNet/povDuplicateWB.aspx"
driver.get(url)
driver.find_element_by_id("txtPovertyLine").clear()
povline = driver.find_element_by_id("txtPovertyLine")
povline.send_keys(str(InputPov))


driver.find_element_by_id("SubmitValue").click()



#Need to wait until page goes. 
print('Waiting 10 seconds until page is loaded correctly')
time.sleep(4)
soup1=BeautifulSoup(driver.page_source, 'lxml')
soup2=soup1.find(id="aggrOutputArea")



table_tag = soup2.select("table")[0]
tab_data = [[item.text for item in row_data.select("th,td")]
                for row_data in table_tag.select("tr")]


t2=tab_data[3]
t3=t2[56]
t4=float(t3)

print('Headcount poverty rate (%) in 2015 was '+t3+' for poverty line selected of '+str(InputPov))



outfile = open("table_data.csv","w")
writer = csv.writer(outfile)


for data in tab_data:
	writer.writerow(data)
	print(' '.join(data))



