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

# app.config["MAIL_SERVER"] = 'Exchange'
# app.config["MAIL_PORT"] = 465
# app.config["MAIL_USERNAME"] = 'neuroeddu.mni@mcgill.ca'
# app.config["MAIL_PASSWORD"] = ''
# app.config['MAIL_USE_SSL'] = True

app.config['MAIL_SERVER']='smtp.mailtrap.io'
app.config['MAIL_PORT'] = 2525
app.config['MAIL_USERNAME'] = '011bcce9b3c8a9'
app.config['MAIL_PASSWORD'] = '95d90cbf56699b'
app.config['MAIL_USE_TLS'] = True
app.config['MAIL_USE_SSL'] = False

#https://stackoverflow.com/questions/60910104/smtp-authentication-error-while-while-sending-mail-from-outlook-using-python-lan

app.config['UPLOAD_FOLDER'] = 'static/files/'

app.config['MAX_CONTENT_LENGTH'] = 200 * 1024 * 1024 # 200MB limit

mail.init_app(app)

from application import routes
