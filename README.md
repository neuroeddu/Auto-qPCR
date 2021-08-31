# Auto-qPCR
A python program to process spreadsheet output directly from qPCR thermocycler.

See the web app version https://auto-q-pcr.com/

Python scripts to run the program our found in this repository. 

# To run a local GUI
1. navigate to the folder 'local-GUI'. You must be within the local-GUI folder or the program won't run.
$ cd local-GUI/
$ ls 
application  main.py  requirements.txt  venv


2. in command line type: python main.py  - this will run the program
$ python main.py
 * Serving Flask app "application" (lazy loading)
 * Environment: production
   WARNING: This is a development server. Do not use it in a production deployment.
   Use a production WSGI server instead.
 * Debug mode: off
 * Running on http://127.0.0.1:5000/ (Press CTRL+C to quit)


either a local browers will open or you will have the link in your terminal click open link and the GUI will open in your webbrowser (firefox, chrome).

All the source scripts are available. Python scripts using flask to create a local server which works to process qPCR data and perform statistics.

The file requirments.txt list packages needed to run the local server.

3. In your local web browser follow the instructions and enter all the user input boxes. 



The origninal version run in command line using a configuration file.  
# Command Line python script instructions

The initial version of this analysis program using command line can be found in the folder 'auto-q-pcr-command_line'.

The document 'Perpare your computer pythonScript.docx' details how to install the required packages.

Steps to run a qPCR analysis from the command line. You must install python and required packages.

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

Program conception: Rhalena Thomas and Gilles Maussion
Program and web app design and management: Rhalena Thomas
Command line data input and absolute model: Iveta Demirova
Relative models and genomic instability: Eddie Cai  
Web app and interface development, statistic and plotting: Gracia Gu

This repository is managed by Rhalena Thomas.  The auto-qPCR webapp is maintained by Rhalena Thomas. 

