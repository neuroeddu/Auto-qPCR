# Auto-qPCR

Step to run a qPCR analysis

1. a) save the raw qPCR output data as a .csv file. Place the file(s) in a data folder.

   b) copy the project_config.conf file into the data folder.
   
   c) If you are using different Control Genes or a different sample order than the example file changet hose to match the experiment you want to analysis.  Save the file. 
2. open command window (do they by typing cmd in your windows bar).
3. navigate to the folder with the program files. Type: "cd /Desktop/Auto-qPCR/program/idlab/
4. type "python idlab.py"      (this will run the program)
5. A window will open asking you to select the location of your raw data.  Navigate in the window to select your data folder and pick OK.

The analysis will run and the plots will all be displayed. The plots will be save and an xls file with the output data will be saved in the "data" folder.
