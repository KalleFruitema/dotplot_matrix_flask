from PIL import Image
import sys
from random import randint


def make_color_matrix_image(norm_values, col_length, row_length):
    all_values = [(val, val, val) for val in norm_values]
    ratio = row_length / col_length
    #TODO maak images en zet al gemaakte images in database??
    im = Image.new("RGB", (col_length, row_length))
    im.putdata(all_values)
    im.save("img_db/image.png", "PNG")
    #TODO maak


def make_grayscale_random(length):
    return [randint(0, 256) for _ in range(length)]


def main():
    length = 500*300
    vals = make_grayscale_random(length)
    make_color_matrix_image(vals, 500, 300)


if __name__ == "__main__":
    main()
