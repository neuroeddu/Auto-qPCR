from flask import Flask
from flask_mail import Mail

mail = Mail()

app = Flask(__name__)

app.static_folder = 'static'
app.secret_key = b'6Le3c7oZAAAAADF6jktQ2xuxnb1I1tlODKQwaWxU'

app.config["MAIL_SERVER"] = 'smtp.gmail.com'
app.config["MAIL_PORT"] = 465
app.config["MAIL_USERNAME"] = 'autoqpcr@gmail.com'
app.config["MAIL_PASSWORD"] = 'neuroeddu'
app.config['MAIL_USE_SSL'] = True

app.config['MAX_CONTENT_LENGTH'] = 200 * 1024 * 1024 # 200MB limit

mail.init_app(app)

from application import routes