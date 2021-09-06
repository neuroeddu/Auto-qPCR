from flask import Flask
from flask_mail import Mail

mail = Mail()

app = Flask(__name__)

app.static_folder = 'static'
app.secret_key = b'6Leez0scAAAAABp0A4Fi85tPDhz70dgRVXcX4g3f'

app.config["MAIL_SERVER"] = 'smtp.gmail.com'
app.config["MAIL_PORT"] = 465
app.config["MAIL_USERNAME"] = 'autoqpcr@gmail.com'
app.config["MAIL_PASSWORD"] = 'neuroeddu'
app.config['MAIL_USE_SSL'] = True

app.config['UPLOAD_FOLDER'] = 'static/files/'

app.config['MAX_CONTENT_LENGTH'] = 200 * 1024 * 1024 # 200MB limit

mail.init_app(app)

from application import routes
