from markupsafe import Markup
from math import floor
from pprint import pprint

from flask import Flask, render_template, url_for, redirect, \
    session, request

from packages.alignment_tool import ScoreMatrix, AlignMatrix, Score


app = Flask(__name__)
app.secret_key = b'\xf0?a\x9a\\\xff\xd4;\x0c\xcbHi'


@app.route('/', methods=['GET', 'POST'])
def home():
    session.clear()
    if request.method == "POST":
        session["seq1"] = request.form['seq1']
        session["seq2"] = request.form['seq2']
        return redirect(url_for('result'))
    return render_template("home.html")


def normalise(values):
    num_values = [val for val in values if type(val) == int]
    old_min = min(num_values)
    old_range = max(num_values) - old_min

    new_min = -255
    new_range = 0 + 0.99999999999 - new_min

    norm = [abs(floor((n - old_min) / old_range * new_range + new_min))
            if type(n) == int else n for n in values]
    return norm


@app.route('/result')
def result():
    scorematrix = ScoreMatrix('static/matrixes/BLOSUM62.csv')
    alignmatrix = AlignMatrix(session["seq1"],
                              session["seq2"],
                              scorematrix)
    all_values = [int(val.__repr__()) if type(val) == Score else val
                  for row in alignmatrix.matrix 
                  for val in row]
    new_values = normalise(all_values)
    color_matrix = []
    col_length = len(alignmatrix.matrix)
    row_length = len(alignmatrix.matrix[0]) 
    for i in range(0, (col_length + 1) * row_length, row_length):
        color_matrix.append(Markup("<td class=\"dotplot_cell\""
                                   f"style=\"background: rgb({val}, {val}, "
                                   f"{val});\"></td>") 
                                   if type(val) == int else 
                                   Markup("<td class=\"dotplot_cell\">"
                                          f"</td>") 
                                   for val in new_values[i:i + row_length])
    return render_template('result.html',
                           data=color_matrix)


if __name__ == '__main__':
    app.run(debug=True)
