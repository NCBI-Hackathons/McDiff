from flask import render_template, redirect, url_for
from app import app
from app.forms import FaddForm

@app.route('/')
def index():
    form = FaddForm()
    return(render_template('index.html', title='Home', form=form))

@app.route('/results/', methods=['GET', 'POST'])
def results():
    form = FaddForm()
    if form.validate_on_submit():
        pass
    else:
        return(redirect(url_for('index')))


@app.route('/acknowledgements/')
def acknowledgements():
    return(render_template('acknowledgements.html', title='Acknowledgements'))
