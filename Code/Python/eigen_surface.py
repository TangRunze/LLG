#!/usr/bin/env python

# Copyright 2016 NeuroData (http://neurodata.io)
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

# qa_fibers.py
# Created by Vikram Chandrashekhar.
# Email: Greg Kiar @ gkiar@jhu.edu

import numpy as np
import nibabel as nib
import os
import vtk
import csv
import math

from dipy.viz import window, actor
from argparse import ArgumentParser


def visualize(atlasfile, outdir, intensityfile):
    """
    Takes fiber streamlines and visualizes them using DiPy
    Required Arguments:
        - atlasfile: Path to atlas file
        - outdir: Path to output directory
        - opacity: Opacity of overlayed brain
    Optional Arguments:
        - fname: name of output file. default is None (fname based on input
          fibfile name)
    """
    intensities, signs = parse_csv(intensityfile)
    # load atlas file
    atlas_volume = load_atlas(atlasfile, intensities, signs)

    faces = [( 1, 0, 0),
             (-1, 0, 0),
             ( 0, 1, 0),
             ( 0,-1, 0),
             ( 0, 1, 0),
             ( 0, 0, 1),
             ( 0, 0,-1),
             ( 1, 1, 0),
             ( 1, 0, 1),
             ( 0, 1, 1),
             (-1,-1, 0),
             (-1, 0,-1),
             ( 0,-1,-1),
             ]

    for i in range(len(faces)):
        # Initialize renderer
        renderer = window.Renderer()
        
        renderer.background((1, 1, 1))
        
        # Set camera orientation properties
        # TODO: allow this as an argument
        renderer.set_camera(position=faces[i])  # args are: position=(), focal_point=(), view_up=()
        
        # Add streamlines to viz session
        renderer.add(atlas_volume)

        # Saves file, if you're into that sort of thing...
        fname = os.path.split(atlasfile)[1].split('.')[0] + str(i) + '.png'
        window.record(renderer, out_path=outdir + fname, size=(600, 600))
    print('done')


def load_atlas(path, intensities, signs):
    '''
    path: path to atlas file
    opacity: opacity of overlayed atlas brain
    '''
    nifti_reader = vtk.vtkNIFTIImageReader()
    nifti_reader.SetFileName(path)
    nifti_reader.Update()


    # intensities = intensities/math.sqrt(np.mean(intensities**2))

    # The following class is used to store transparencyv-values for later
    # retrival. In our case, we want the value 0 to be completly opaque
    alphaChannelFunc = vtk.vtkPiecewiseFunction()
    alphaChannelFunc.AddPoint(0, 0.0)
    for i in range(len(intensities)):
        alphaChannelFunc.AddPoint(i+1, intensities[i])

    # This class stores color data and can create color tables from a few color
    # points. For this demo, we want the three cubes to be of the colors red
    # green and blue.
    colorFunc = vtk.vtkColorTransferFunction()
    colorFunc.AddRGBPoint(0, 0.0, 0.0, 0.0)
    for i in range(len(signs)):
        if signs[i] < 0: colorFunc.AddRGBPoint(i+1, 0.0, 0.0, intensities[i])
        elif signs[i] > 0: colorFunc.AddRGBPoint(i+1, intensities[i], 0.0, 0.0)
        else: colorFunc.AddRGBPoint(i+1, 1.0, 1.0, 1.0)

    # The previous two classes stored properties. Because we want to apply
    # these properties to the volume we want to render, we have to store them
    # in a class that stores volume prpoperties.
    volumeProperty = vtk.vtkVolumeProperty()
    volumeProperty.SetColor(colorFunc)
    volumeProperty.SetScalarOpacity(alphaChannelFunc)
    volumeProperty.ShadeOn()

    # We can finally create our volume. We also have to specify the data for
    # it, as well as how the data will be rendered.
    volumeMapper = vtk.vtkSmartVolumeMapper()
    volumeMapper.SetInputDataObject(nifti_reader.GetOutput())

    # The class vtkVolume is used to pair the preaviusly declared volume as
    # well as the properties to be used when rendering that volume.
    volume = vtk.vtkVolume()
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperty)

    return volume

def parse_csv(csvfile):
    with open(csvfile, "rb") as f:
        reader = csv.reader(f)
        intensities = [float(col) for row in reader for col in row]
    signs = np.sign(intensities) 
    return normalize_intensities(intensities), signs

def normalize_intensities(intensities):
    intensities = np.absolute(intensities)
    minval = intensities.min()
    maxval = intensities.max()
    intensities = (intensities - minval)/maxval
    print intensities
    return intensities

def main():
    """
    Argument parser. Takes organization and atlas
    information and produces a dictionary of file lists based on datasets
    of interest and then passes it off for processing.
    Required parameters:
        atlasfile:
            - Path to regdti file
        opacity:
            - Opacity of overlayed atlas file
        outdir:
            - Path to fa png save location
    """
    parser = ArgumentParser(description="Generates a fiber png based on fibers"
                            " and atlas")
    parser.add_argument("atlasfile", action="store", help="base directory loc")
    parser.add_argument("outdir", action="store", help="base directory loc")
    parser.add_argument("intensities", action="store", help="array of intensities")
    result = parser.parse_args()

    visualize(result.atlasfile, result.outdir, result.intensities)

if __name__ == "__main__":
    main()
