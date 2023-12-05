from markupsafe import Markup
from math import floor
from pprint import pprint
import json
from typing import Any

from flask import Flask, render_template, url_for, redirect, session, request
from PIL import Image
from uuid import uuid4
from werkzeug.wrappers.response import Response

from packages.alignment_tool import ScoreMatrix, AlignMatrix, Score


app = Flask(__name__)
app.secret_key = b'\xf0?a\x9a\\\xff\xd4;\x0c\xcbHi'


def try_get_request(varname) -> bool:
    try:
        if request.form[varname] == "":
            return False
        session[varname] = request.form[varname]
        return True
    except Exception:
        return False


@app.route('/', methods=['GET', 'POST'])
def home() -> str | Response:
    session.clear()
    error = ""
    if request.method == "POST":
        for var in ["seq1", "seq2", "matrix_type", "seqtype"]:
            succes = try_get_request(var)
            if not succes:
                error = f"\"{var}\" entered incorrectly!"
                return render_template("home.html",
                                       error=error)
        if session["matrix_type"] == "direct_check":
            session["pairwise_type"] = None
        else:
            try:
                session["pairwise_type"] = request.form["pairwise_type"]
            except Exception:
                session["pairwise_type"] = "N-W"

        return redirect(url_for('result'))
    return render_template("home.html",
                           error="")


def get_scorematrix() -> ScoreMatrix:
    if session["seqtype"] == "ntseq":
        scorematrix_path = 'static/matrixes/kimura_matrix.csv'
    elif session["seqtype"] == "aaseq":
        scorematrix_path = 'static/matrixes/BLOSUM62.csv'
    scorematrix = ScoreMatrix(scorematrix_path)
    return scorematrix


def normalise(values) -> list[Any]:
    num_values = [val for val in values if type(val) in [int, float]]
    old_min = min(num_values)
    old_range = max(num_values) - old_min

    new_min = -255
    new_range = 0 + 0.99999999 - new_min

    norm = [abs(floor((n - old_min) / old_range * new_range + new_min))
            if type(n) in [int, float] else n for n in values]
    return norm


# OUTDATED
def make_color_matrix_html(norm_values, col_length, row_length) -> list[list[Markup]]:
    return [[Markup("<td class=\"dotplot_cell\" "
                    f"style=\"background: rgb({val}, {val},"
                    f" {val});\"></td>")
                    for val in norm_values[i:i + row_length]
                    if type(val) == int]
                    for i in range(0, (col_length + 1) * 
                                    row_length, row_length)]


def add_image_to_json(img_id) -> None:
    with open("static/img_db/img_db_lookup.json", "r") as json_file_read:
        json_db = json.load(json_file_read)

    json_db['alignment_images'].append({
        "image": f"{img_id}.png",
        "matrix_type": session["matrix_type"],
        "pairwise_type": session["pairwise_type"],
        "sequences": [session["seq1"], session["seq2"]],
    })
    with open("static/img_db/img_db_lookup.json", "w") as json_file_write:
        json.dump(json_db, json_file_write, indent=4)


def make_color_matrix_image(norm_values, col_length, row_length) -> str:
    all_values = [(val, val, val) for val in norm_values
                  if type(val) == int]
    ratio = row_length / col_length

    im = Image.new("RGB", (row_length, col_length))
    im.putdata(all_values, scale=ratio)
    img_id = uuid4()
    img_url = f"static/img_db/img/{img_id}.png"
    im.save(img_url, "PNG")
    add_image_to_json(img_id)
    return img_url


def pairwise() -> str:
    scorematrix = get_scorematrix()
    alignmatrix = AlignMatrix(session["seq1"], session["seq2"],
                              scorematrix, aligntype=session["pairwise_type"])

    all_values = [((int(val.__repr__()) + 0.00000000001) / 
                   (val.get_steps() + 1))
                  if type(val) == Score else val
                  for row in alignmatrix.matrix 
                  for val in row]
    new_values = normalise(all_values)
    col_length, row_length = \
        len(alignmatrix.matrix), len(alignmatrix.matrix[0])
    # color_matrix = make_color_matrix_html(new_values, col_length, row_length)
    color_matrix = make_color_matrix_image(new_values, col_length-1, row_length-1)
    return color_matrix


def direct_check() -> str:
    scorematrix = get_scorematrix()
    all_values = [int(scorematrix.get_score(nuc1, nuc2))
                  for nuc1 in session["seq2"]
                  for nuc2 in session["seq1"]]

    new_values = normalise(all_values)
    col_length, row_length = \
        len(session["seq2"]), len(session["seq1"])
    # color_matrix = make_color_matrix_html(new_values, col_length, row_length)
    color_matrix = make_color_matrix_image(new_values, col_length, row_length)
    return color_matrix 


def check_json_db() -> str | None:
    with open("static/img_db/img_db_lookup.json", "r") as json_file:
        json_db = json.load(json_file)
    for img in json_db["alignment_images"]:
        if img["matrix_type"] == session["matrix_type"] and\
        img["sequences"][0] == session["seq1"] and\
            img["sequences"][1] == session["seq2"] and\
                img["pairwise_type"] == session["pairwise_type"]:
            print("Query in database already!")
            return f'static/img_db/img/{img["image"]}'
    return None
            


@app.route('/result')
def result() -> str:
    color_matrix = check_json_db()
    if color_matrix is None:
        if session["matrix_type"] == 'pairwise':
            color_matrix = pairwise()
        elif session["matrix_type"] == 'direct_check':
            color_matrix = direct_check()
    print(color_matrix)
    return render_template('result.html',
                           image_path=Markup(f'<img src=\"../{color_matrix}\">'))


if __name__ == '__main__':
    app.run(debug=True)
 