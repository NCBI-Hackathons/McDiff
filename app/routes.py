from flask import render_template, redirect, url_for
from app import app
from app.forms import FaddForm
import parp_simulator
from werkzeug.utils import secure_filename
import os

@app.route('/')
def index():
    form = FaddForm()
    return(render_template('index.html', title='Home', form=form))

@app.route('/results/', methods=['GET', 'POST'])
def results():
    form = FaddForm()
    print(form.errors)
    if form.validate_on_submit():
        # CSV
        f = form.csv.data
        filename = secure_filename(f.filename)
        csv_file = os.path.join(app.root_path, 'static', 'uploads', filename)
        f.save(csv_file)
        # ROI
        f = form.roi_file.data
        filename = secure_filename(f.filename)
        roi_file = os.path.join(app.root_path, 'static', 'uploads', filename)
        f.save(roi_file)
        # MASK
        f = form.mask_file.data
        filename = secure_filename(f.filename)
        mask_file = os.path.join(app.root_path, 'static', 'uploads', filename)
        f.save(mask_file)
    else:
        return(redirect(url_for('index')))


@app.route('/acknowledgements/')
def acknowledgements():
    return(render_template('acknowledgements.html', title='Acknowledgements'))
