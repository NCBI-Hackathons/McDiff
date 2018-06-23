from flask_wtf import FlaskForm
from wtforms import SubmitField, FloatField, IntegerField
from wtforms.validators import DataRequired, NumberRange
from flask_wtf.file import FileField, FileRequired, FileAllowed
from werkzeug import secure_filename

class FaddForm(FlaskForm):

    # CSV File
    csv = FileField(validators=[FileRequired(), FileAllowed(['csv'])])

    # File 2 
    roi_file = FileField(validators=[FileRequired(), FileAllowed(['txt'])])

    # File 3 
    mask_file = FileField(validators=[FileRequired(), FileAllowed(['txt'])])

    # bound_d
    bound_d = IntegerField('bound_d', validators=[DataRequired(message="This field is required and must be an int between 0 and 100"), NumberRange(min=1, max=100)], render_kw={"value": "20"})
    
    # exp_time
    exp_time = FloatField('exp_time', validators=[DataRequired(message="This field is required and must be an int between 0 and 60"), NumberRange(min=0, max=15)])


    # sigmaD
    sigma_d = IntegerField('sigma_d', validators=[DataRequired(message="This field is required and must be an int between 0 and 100"), NumberRange(min=1, max=100)], render_kw={"value": "2"})


    # sigmaF
    sigma_f = FloatField('sigma_f', validators=[DataRequired(message="This field is required and must be an int between 0 and 60"), NumberRange(min=0, max=15)], render_kw={"value": ".05"})

    
    # mcm_temp
    mcmc_temp = IntegerField('mcmc_temp', validators=[DataRequired(message="This field is required and must be an int between 0 and 100"), NumberRange(min=1, max=100)], render_kw={"value": "1"})


    # Offset
    offset = FloatField('Offset', validators=[DataRequired(message="This field is required and must be an int between 0 and 60"), NumberRange(min=0, max=15)], render_kw={"value": "10.5"})

    # Steps
    steps = IntegerField('Steps', validators=[DataRequired(message="This field is required and must be an int between 0 and 600"), NumberRange(min=1, max=600)], render_kw={"value": "200"})

    # Submit button
    submit = SubmitField('Run!', render_kw={"class": "btn btn-secondary"})

