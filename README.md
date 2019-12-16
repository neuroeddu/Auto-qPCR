# Auto-qPCR

Currently the website folder contains python script using flask to create a local server

We are deploying that program using google cloud computing and the scripts for that process are in the "cloud" folder.

A user guide in in progress.







# Orginal python script instructions

The python script for to run processing from the terminal is still avialable

Step to run a qPCR analysis

1. a) save the raw qPCR output data as a .csv file. Place the file(s) in a data folder.

   b) copy the project_config.conf file into the data folder.
   
   c) Update the project_config.conf file to match your situation:
   
         1. choose your model: absolute, relative (delta CT), stability (delta delta CT)
   
         2. Specify Control Genes (house keeping/ configuration genes)
         
         3. If you are using delta-delta CT (the stability test) specify your reference sample.
         
         4. Save the file. 
2. open command window (by typing cmd in your windows bar).
3. navigate to the folder with the program files. Type: "cd /Desktop/Auto-qPCR/program/idlab/
4. type "python idlab.py"      (this will run the program)
5. A window will open asking you to select the location of your raw data.  Navigate in the window to select your data folder and pick OK.

The analysis will run and the plots will all be displayed. The plots will be save and an xls file with the output data will be saved in the "data" folder.
