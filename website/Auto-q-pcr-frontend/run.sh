#!/bin/bash
python -m virtualenv venv
source venv/Scripts/activate 
python -m pip install -r requirements.txt
python -m webbrowser http://127.0.0.1:5000/
python main.py
#$SHELL