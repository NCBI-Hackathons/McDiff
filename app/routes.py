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
        # This is just for saving files and getting our name. Gosh
        # CSV
        f = form.csv.data
        filename = secure_filename(f.filename)
        csv_file = os.path.join(app.root_path, 'static', 'uploads', filename)
        f.save(csv_file)
        print(os.path.basename(csv_file))
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
        # Now we can actally call wrapper with our form data and files
        fit_plot_name, resid_plot_name, fit_data_name, d_final, f_final, ret_error = wrapper(csv_file, roi_file, mask_file, form.bound_d.data, form.exp_time.data, form.percent_bleached.data, form.sigma_d.data, form.sigma_f.data, form.mcmc_temp.data, form.offset.data, form.steps.data)
    else:
        return(redirect(url_for('index')))


@app.route('/acknowledgements/')
def acknowledgements():
    return(render_template('acknowledgements.html', title='Acknowledgements'))
