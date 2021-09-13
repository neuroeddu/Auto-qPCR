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


from flask import Flask
from flask_mail import Mail

mail = Mail()

app = Flask(__name__)

app.static_folder = 'static'
app.secret_key = b'6Leez0scAAAAABp0A4Fi85tPDhz70dgRVXcX4g3f'

app.config["MAIL_SERVER"] = 'Exchange'
app.config["MAIL_PORT"] = 465
app.config["MAIL_USERNAME"] = 'neuroeddu.mni@mcgill.ca'
app.config["MAIL_PASSWORD"] = '*********'
app.config['MAIL_USE_SSL'] = True

app.config['UPLOAD_FOLDER'] = 'static/files/'

app.config['MAX_CONTENT_LENGTH'] = 200 * 1024 * 1024 # 200MB limit

mail.init_app(app)

from application import routes
