"""
Author: Mikhail Schee

Creates a gif of a set of png images in the stated directory.

Usage:
    create_gif.py NAME FILENAME FRAMES_PATH

Arguments:
    NAME            name of experiment
    FILENAME        name for gif file
    FRAMES_PATH     path to frames of animation

Script created: 2019/03/19, Mikhail Schee
Last updated: 2020/01/07, Mikhail Schee
"""

"""
Using imageio
https://imageio.readthedocs.io/en/latest/installation.html
install using:
    conda install -c conda-forge imageio
OR
    pip install imageio
Skeleton script from
https://stackoverflow.com/questions/41228209/making-gif-from-images-using-imageio-in-python
"""

import os
import imageio
# For adding arguments when running
from docopt import docopt

gif_file = 'uh_oh.gif'

if __name__ == '__main__':
    arguments = docopt(__doc__)
    NAME = arguments.get('NAME')
    gif_file = arguments.get('FILENAME')
    png_dir = arguments.get('FRAMES_PATH')
    print('Gif saving to', gif_file)

import switchboard as sbp

# To set the frames per second
frame_duration = 1.0 / sbp.fps

images = []
# need to sort because os.listdir returns a list of arbitrary order
for file_name in sorted(os.listdir(png_dir)):
    if file_name.endswith('.png'):
        file_path = os.path.join(png_dir, file_name)
        images.append(imageio.imread(file_path))
imageio.mimsave(gif_file, images, duration=frame_duration)
