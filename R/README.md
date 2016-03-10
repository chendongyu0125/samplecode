# R

#scrarticles

This file scrapes information from the website of the Administrative Department of Science, Technology and Innovation for Colombia (COLCIENCIAS) in order to obtain the number of articles published by each university in Colombia by year, category, country published, and more categories. 

#gmapsdistance

The function `gmapsdistance` uses the [Google Maps Distance Matrix API](https://developers.google.com/maps/documentation/distance-matrix/intro?hl=en) to compute the distance and time between two points. An [API key](https://developers.google.com/maps/documentation/distance-matrix/get-api-key#key) is not necessary to perform the query but the function supports its usage. If an API key is being used the Distance Matrix API should be enabled in the Google Developers Console. Google maps must be able to find both the origin and the destination in order for the function to run. If the origin or destination contains multiple words, they should be separated by a plus sign (+). The distance is returned in meters and the time in seconds. 
More information about gmapsdistance available in the [github website](https://github.com/rodazuero/gmapsdistance) and the [website at the CRAN-repository](https://cran.r-project.org/web/packages/gmapsdistance/index.html).