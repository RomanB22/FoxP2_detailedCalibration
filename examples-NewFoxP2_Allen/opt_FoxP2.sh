#!/bin/bash

nrnivmodl ./mechanisms
python opt_FoxP2.py --start
