#!/usr/bin/python
import sys
import logging
logging.basicConfig(stream=sys.stderr)

sys.path.append("/var/www/Auto-qPCR/website/Auto-q-pcr-frontend/")
activate_this = '/var/www/Auto-qPCR/website/Auto-q-pcr-frontend/venv/bin/activate_this.py'
with open(activate_this) as file_:
    exec(file_.read(), dict(__file__=activate_this))

from application import app as application
application.secret_key = 'Add your secret key'
