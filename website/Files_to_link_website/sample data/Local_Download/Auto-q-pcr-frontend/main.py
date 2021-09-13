# Auto-q-PCR is a program for analysis of qPCR data for absolute and relative quantification
# Copyright (C) 2021 Rhalena Thomas, Eddie Cai, Gracia Gu and Iva Demirova
#This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation,  version 3 of the License.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.


#########################################################################################################

# application is a directory containing scripts for each model and the AutoqPCR.py script which reads in the arguments from the GUI and runs the appropriate scripts

from application import app

if __name__ == '__main__':
	app.run(debug=False) 

