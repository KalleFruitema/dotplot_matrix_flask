from markupsafe import Markup
from math import floor
from pprint import pprint

from flask import Flask, render_template, url_for, redirect, session, request

from packages.alignment_tool import ScoreMatrix, AlignMatrix, Score


app = Flask(__name__)
app.secret_key = b'\xf0?a\x9a\\\xff\xd4;\x0c\xcbHi'


@app.route('/', methods=['GET', 'POST'])
def home():
    session.clear()
    if request.method == "POST":
        session["seq1"] = request.form['seq1']
        session["seq2"] = request.form['seq2']
        session["matrix_type"] = request.form['matrix_type']
        session["seqtype"] = request.form['seqtype']
        return redirect(url_for('result'))
    return render_template("home.html")


def get_scorematrix():
    if session["seqtype"] == "ntseq":
        scorematrix_path = 'static/matrixes/kimura_matrix.csv'
    elif session["seqtype"] == "aaseq":
        scorematrix_path = 'static/matrixes/BLOSUM62.csv'
    scorematrix = ScoreMatrix(scorematrix_path)
    return scorematrix


def normalise(values):
    num_values = [val for val in values if type(val) in [int, float]]
    old_min = min(num_values)
    old_range = max(num_values) - old_min

    new_min = -255
    new_range = 0 + 0.99999999 - new_min

    norm = [abs(floor((n - old_min) / old_range * new_range + new_min))
            if type(n) in [int, float] else n for n in values]
    return norm


def make_color_matrix(norm_values, col_length, row_length):
    color_matrix = []
    for i in range(0, (col_length + 1) * row_length, row_length):
        color_matrix.append([Markup("<td class=\"dotplot_cell\""
                                    f"style=\"background: rgb({val}, {val},"
                                    f" {val});\"></td>")
                                   for val in norm_values[i:i + row_length]
                                   if type(val) == int])
    return color_matrix


def pairwise():
    scorematrix = get_scorematrix()
    alignmatrix = AlignMatrix(session["seq1"], session["seq2"],
                              scorematrix, aligntype='N-W')
    all_values = [((int(val.__repr__()) + 0.00000000001) / 
                   (val.get_steps() + 1))
                  if type(val) == Score else val
                  for row in alignmatrix.matrix 
                  for val in row]
    new_values = normalise(all_values)
    col_length, row_length = \
        len(alignmatrix.matrix), len(alignmatrix.matrix[0])
    color_matrix = make_color_matrix(new_values, col_length, row_length)
    return color_matrix


def direct_check():
    scorematrix = get_scorematrix()
    all_values = [int(scorematrix.get_score(nuc1, nuc2))
                  for nuc1 in session["seq2"]
                  for nuc2 in session["seq1"]]

    new_values = normalise(all_values)
    col_length, row_length = \
        len(session["seq2"]), len(session["seq1"])
    color_matrix = make_color_matrix(new_values, col_length, row_length)
    return color_matrix 


@app.route('/result')
def result():
    if session["matrix_type"] == 'pairwise':
        color_matrix = pairwise()
    elif session["matrix_type"] == 'direct_check':
        color_matrix = direct_check()
    return render_template('result.html',
                           data=color_matrix)


if __name__ == '__main__':
    app.run(debug=True)
