#!/bin/bash
python3 -m venv local-env
source loca-env/bin/activate 
python3 -m pip install -r requirements.txt
python3 main.py
#$SHELL
