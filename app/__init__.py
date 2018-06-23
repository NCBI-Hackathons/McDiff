from flask import Flask
from config import Config
import os

app = Flask(__name__)
app.config.from_object(Config)

from app import routes

os.path.abspath(os.path.dirname(__file__))
if not os.path.exists("{0}/static/img".format(os.path.abspath(os.path.dirname(__file__)))):
    os.makedirs("{0}/static/img".format(os.path.abspath(os.path.dirname(__file__))))
if not os.path.exists("{0}/static/uploads".format(os.path.abspath(os.path.dirname(__file__)))):
    os.makedirs("{0}/static/uploads".format(os.path.abspath(os.path.dirname(__file__))))


