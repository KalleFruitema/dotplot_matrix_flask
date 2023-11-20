from flask import Flask, render_template
from packages.alignment_tool import ScoreMatrix, AlignMatrix


app = Flask(__name__)


@app.route('/')
def home():
    
    return render_template("home.html")


if __name__ == '__main__':
    app.run()
