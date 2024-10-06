#!/usr/bin/python
import sys
sys.path.insert(0,"/var/www/exosky/")
sys.path.append("/home/ubuntu/.local/lib/python3.8/site-packages")
sys.path.append("/usr/local/lib/python3.8/dist-packages")

from app import app as application
