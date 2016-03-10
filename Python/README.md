# Python
#ACCESCREDITS.py
This file scrapes data from the website http://spadies.mineducacion.gov.co/spadies/JSON.html containing information about higher education in Colombia. 
Note that the requests for the web scrape are done via "post" rather than "get" requests. 

#gmapsdistance.py
The function gmapsdistance uses the Google Maps Distance Matrix API to compute the distance and time between two points. An API key is not necessary to perform the query but the function supports its usage. If an API key is being used the Distance Matrix API should be enabled in the Google Developers Console. Google maps must be able to find both the origin and the destination in order for the function to run. If the origin or destination contains multiple words, they should be separated by a plus sign (+). The distance is returned in meters and the time in seconds.
Four different modes of transportation are allowed: bicycling, walking, driving, transit.
I am the author and maintainer of the gmapsdistance package in R. This version in Python allows to compute a matrix of distances between points. Say you have a point of "N" locations nodes and you want to compute all the distance and travel time between each pair of points. This version in Python allows such a computation. Note that the input should be a file containing "N" locations, the output will be a csv file containing "N*(N+1)/2" distances or travel times between each pair of locations. 
As an example, I include the file "Cluster8INPUT.csv" which includes geographic locations in Ecuador. By running this file you will get the walking distance between each pair of nodes. It can be modified to get the distance by car, bike or public transportation and the corresponding time. 

