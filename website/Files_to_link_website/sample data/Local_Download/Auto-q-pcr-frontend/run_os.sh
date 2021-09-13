#!/bin/bash
python3 -m venv local-env
source loca-env/bin/activate 
python3 -m pip install -r requirements.txt
python -m webbrowser http://127.0.0.1:5000/
python3 main.py
#$SHELL
