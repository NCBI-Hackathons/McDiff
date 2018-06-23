from flask_wtf import FlaskForm
from wtforms import SubmitField, FloatField, IntegerField
from wtforms.validators import DataRequired, NumberRange
from flask_wtf.file import FileField, FileRequired, FileAllowed
from werkzeug import secure_filename

class FaddForm(FlaskForm):

    csv = FileField(validators=[FileRequired(), FileAllowed(['csv'])])
    offset = FloatField('Offset', validators=[DataRequired(message="This field is required and must be an int between 0 and 60"), NumberRange(min=0, max=15)], render_kw={"value": "10.5"})

    steps = IntegerField('Steps', validators=[DataRequired(message="This field is required and must be an int between 0 and 600"), NumberRange(min=1, max=600)], render_kw={"value": "60"})

    # Submit button
    submit = SubmitField('Run!', render_kw={"class": "btn btn-secondary"})

