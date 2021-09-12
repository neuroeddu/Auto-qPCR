#!/bin/bash
python3 -m virtualenv venv
source venv venv/bin/activate 
python3 -m pip install -r requirements.txt
python3 -m webbrowser http://127.0.0.1:5000/
python3 main.py
#$SHELL
