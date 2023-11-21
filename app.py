from flask import Flask, render_template, request, redirect, url_for, session
from packages.alignment_tool import ScoreMatrix, AlignMatrix


app = Flask(__name__)
app.secret_key = b'\xf0?a\x9a\\\xff\xd4;\x0c\xcbHi'


@app.route('/', methods=['GET', 'POST'])
def home():
    if request.method == "POST":
        session["seq1"] = request.form['seq1']
        session["seq2"] = request.form['seq2']
        return redirect(url_for('result'))
    return render_template("home.html")


@app.route('/result')
def result():
    scorematrix = ScoreMatrix('static/matrixes/BLOSUM62.csv')
    alignmatrix = AlignMatrix(session["seq1"],
                              session["seq2"],
                              scorematrix)
    return render_template('result.html',
                           data=alignmatrix.matrix)


if __name__ == '__main__':
    app.run(debug=True)
