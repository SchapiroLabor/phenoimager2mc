#!/usr/bin/env python
########################################################################################################
# AUTHOR: Chiara Schiller <chiara.schiller@uni-heidelberg.de>
# DESCRIPTION: Convert PhenoImager output .tif tiles to a single .tif file with correct OME-XML metadata
########################################################################################################
import argparse
import warnings
import copy
import platform
import os

import numpy as np
import tifffile as tiff
from PIL import Image
import xml.etree.ElementTree as ET
import ome_types
from ome_types.model import OME, Pixels, TiffData, Channel, Plane
from ome_types.model.pixels import DimensionOrder, PixelType
from ome_types.model.simple_types import ChannelID, ImageID, PixelsID, Color
from uuid import uuid4
from ome_types import to_xml
from skimage import util

from ._version import __version__


def getOptions(myopts=None):
    """ Function to pull in arguments """
    description = """ PhenoImager2MC """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    standard = parser.add_argument_group(
        title='Standard Inputs',
        description='Standard input for staging module.')
    standard.add_argument(
        "-i",
        "--indir",
        nargs='+',
        required=True,
        help='Folder name or list of files to process')
    standard.add_argument(
        "-m",
        "--num_markers",
        dest="num_markers",
        action="store",
        required=True,
        type=int,
        help="Provide number of markers in the image.")
    standard.add_argument(
        "-n",
        "--normalization",
        dest="normalization",
        action="store",
        required=False,
        default="99th",
        choices=["99th", "max"],
        help="Provide method to normalize marker intensities per channel. Options = [99th, max]")
    standard.add_argument("--version", action="version", version=f"v{__version__}")

    output = parser.add_argument_group(title='Required output')
    output.add_argument(
        "-o",
        "--output",
        dest="output",
        action='store',
        required=True,
        help="Output file, existing or will be newly created.")

    args = parser.parse_args(myopts)
    return args


def concatenate_tiles(input_dir, output_file):
    """
    Stack all .tif tiles per cycle into a single .tif file

    :Arguments:
        :type input_dir: list
        :param input_dir: Folder path (as single-element list) or list of .tif file paths

        :type output_file: str
        :param output_file: Path for the stacked .tif output file

    :Returns:
        :rtype stacked_image: np.ndarray
        :return stacked_image: Stacked image saved as numpy array
    """
    print(input_dir)
    if len(input_dir) == 1:
        folder = input_dir[0]
        input_files = [
            os.path.join(folder, f) for f in os.listdir(folder)
            if f.endswith('.tif')
        ]
    else:
        input_files = input_dir

    input_files = sorted(input_files)

    first_img = np.asarray(tiff.imread(input_files[0]))
    num_tiles = len(input_files)
    num_channels = first_img.shape[0]
    # Initialize an empty numpy array to hold all images
    # Shape: (Tiles, Channels, Y, X)
    stack_shape = (num_tiles, num_channels, first_img.shape[1],
                   first_img.shape[2])
    stacked_image = np.empty(stack_shape, dtype=first_img.dtype)

    for idx, file in enumerate(input_files):
        img = np.asarray(tiff.imread(file))
        stacked_image[idx] = img

    tiff.imwrite(output_file,
                 stacked_image)
    return stacked_image

def compute_normalization_factors(input_files, num_channels, normalization):
    """Stream tiles to compute global per-channel normalization factors."""
    # Accumulate per-channel values for percentile computation
    channel_values = [[] for _ in range(num_channels)]
    
    for file in input_files:
        img = np.asarray(tiff.imread(file))  # shape: (num_channels, Y, X)
        for ch in range(num_channels):
            channel_values[ch].append(img[ch].ravel())
    
    factors = {}
    for ch in range(num_channels):
        all_values = np.concatenate(channel_values[ch])
        if normalization == "99th":
            factors[ch] = np.percentile(all_values, 99)
        elif normalization == "max":
            factors[ch] = np.max(all_values)
        print(f"Channel {ch} normalization factor: {factors[ch]}")
    
    return factors

def normalize_and_write(input_files, output_file, factors, num_channels):
    """Normalize each tile using precomputed factors and write to output."""
    first_img = np.asarray(tiff.imread(input_files[0]))
    num_tiles = len(input_files)
    stack_shape = (num_tiles, num_channels, first_img.shape[1], first_img.shape[2])
    
    # Pre-allocate output array — only one tile in memory at a time during processing
    stacked = np.empty(stack_shape, dtype=np.uint16)
    
    for idx, file in enumerate(input_files):
        img = np.asarray(tiff.imread(file), dtype=np.float32)
        print(img.shape, file)
        for ch in range(num_channels):
            print(img[ch].shape)
            normalized = img[ch] / factors[ch]
            normalized = np.clip(normalized, 0, 1)
            stacked[idx, ch] = util.img_as_uint(normalized)
    tiff.imwrite(output_file, stacked)
    return stacked

def extract_metadata(tif_file_path):
    """
    Extract metadata from PhenoImager corrupted .tif files

    :Arguments:
        :type tif_file_path: str
        :param tif_file_path: File path per tile in one cycle

    :Returns:
        :rtype metadata_list: list
        :return metadata_list: List of extracted metadata per tile across channels.
    """
    metadata_list = []

    with tiff.TiffFile(tif_file_path) as tif:
        for i, page in enumerate(tif.pages):
            page_metadata = {
                "page_index": i,
                "image_width": page.imagewidth,
                "image_length": page.imagelength,
                "bits_per_sample": page.bitspersample,
                "x_position": None,
                "y_position": None,
                "z_position": i,
                "image_description": None,
                "marker": None,
                "exposure_time": None,
                "signal_unit": None,
                "device_name": None,
                "no_of_channels": None,
                "pixel_unit": "PixelsPerCentimeter",
                "dpi": None,
                "x_physical_size": None,
                "y_physical_size": None
            }

            img = Image.open(tif_file_path)
            info = img.info
            page_metadata["dpi"] = info['dpi'][0]

            pixels_per_micrometer = page_metadata["dpi"] / 25400
            micrometers_per_pixel = 1 / pixels_per_micrometer

            x_position = page.tags.get('XPosition')
            y_position = page.tags.get('YPosition')
            if x_position and y_position:
                page_metadata["x_position"] = (x_position.value[0] / x_position.value[1]) * 10000
                page_metadata["y_position"] = (y_position.value[0] / y_position.value[1]) * 10000
                page_metadata["x_physical_size"] = micrometers_per_pixel
                page_metadata["y_physical_size"] = micrometers_per_pixel

            xml_string = page.tags.get('ImageDescription').value
            root = ET.fromstring(xml_string)

            name_element = root.find('./Name')
            if name_element is not None:
                page_metadata["marker"] = name_element.text

            exposure_time_element = root.find('./ExposureTime')
            if exposure_time_element is not None:
                page_metadata["exposure_time"] = exposure_time_element.text

            signal_units_element = root.find('./SignalUnits')
            if signal_units_element is not None:
                page_metadata["signal_unit"] = signal_units_element.text

            scan_profile = root.find('.//ScanProfile')
            if scan_profile is not None:
                pixel_units_element = scan_profile.find('.//root')
                if pixel_units_element is not None:
                    pixel_units_element = scan_profile.find('.//PixelUnits')
                    if pixel_units_element is not None:
                        page_metadata["pixel_unit"] = pixel_units_element.text

            metadata_list.append(page_metadata)

    return metadata_list


def create_ome(all_tile_metadata, output_file, num_markers):
    """
    Creates and overwrites correct omexml metadata for the OME-TIFF files

    :Arguments:
        :type all_tile_metadata: list
        :param all_tile_metadata: List of extracted metadata per tile across channels from extract_metadata().

        :type output_file: str
        :param output_file: Path to the stacked .tif file to annotate with OME-XML.

        :type num_markers: int
        :param num_markers: Number of markers in the image.
    """
    markers = [all_tile_metadata[0][i]['marker'] for i in range(num_markers)]

    no_of_channels = num_markers
    bits_per_sample = all_tile_metadata[0][0]['bits_per_sample']
    no_of_tiles = len(all_tile_metadata)
    tiff.tiffcomment(output_file, '')

    tiff_block = [
        TiffData(first_c=ch, ifd=ch, plane_count=1)
        for ch in range(no_of_channels)
    ]

    plane_block = [
        Plane(the_c=ch, the_t=0, the_z=0, position_x=0, position_y=0, position_z=0, exposure_time=0)
        for ch in range(no_of_channels)
    ]

    chann_block = [
        Channel(
            id=ChannelID('Channel:{x}'.format(x=ch)),
            name=chann_name,
            color=Color((255, 255, 255)),
            emission_wavelength=1,
            excitation_wavelength=1,
        )
        for ch, chann_name in enumerate(markers)
    ]

    pix_block = []
    ifd_counter = 0
    for t in range(no_of_tiles):
        template_plane_block = copy.deepcopy(plane_block)
        template_chann_block = copy.deepcopy(chann_block)
        template_tiffdata_block = copy.deepcopy(tiff_block)
        for ch, mark in enumerate(markers):
            template_plane_block[ch].position_x = all_tile_metadata[t][0]['x_position']
            template_plane_block[ch].position_y = all_tile_metadata[t][0]['y_position']
            template_plane_block[ch].exposure_time = all_tile_metadata[t][ch]['exposure_time']
            template_chann_block[ch].id = 'Channel:{y}:{x}'.format(x=ch, y=100 + t)
            template_chann_block[ch].name = mark
            template_tiffdata_block[ch].ifd = ifd_counter
            ifd_counter += 1
        pix_block.append(Pixels(
            id=PixelsID('Pixels:{x}'.format(x=t)),
            dimension_order=DimensionOrder('XYCZT'),
            size_c=no_of_channels,
            size_t=1,
            size_x=all_tile_metadata[t][0]['image_width'],
            size_y=all_tile_metadata[t][0]['image_length'],
            size_z=1,
            type=PixelType('uint32'),
            big_endian=False,
            channels=template_chann_block,
            interleaved=False,
            physical_size_x=all_tile_metadata[t][0]['x_physical_size'],
            physical_size_y=all_tile_metadata[t][0]['y_physical_size'],
            physical_size_z=1.0,
            planes=template_plane_block,
            bits_per_sample=bits_per_sample,
            tiff_data_blocks=template_tiffdata_block))

    img_block = [
        ome_types.model.Image(
            id=ImageID('Image:{x}'.format(x=t)),
            pixels=pix_block[t])
        for t in range(no_of_tiles)
    ]

    ome_custom = OME()
    ome_custom.creator = " ".join([
        ome_types.__name__, ome_types.__version__,
        '/ python version-', platform.python_version()
    ])
    ome_custom.images = img_block
    ome_custom.uuid = uuid4().urn
    ome_xml = to_xml(ome_custom)
    tiff.tiffcomment(output_file, ome_xml)


def main(args):
    """Core pipeline logic."""
    if len(args.indir) == 1:
        input_dir = args.indir[0]
        input_files = [
            os.path.join(input_dir, f) for f in os.listdir(input_dir)
            if f.endswith('.tif')
        ]
    else:
        input_files = args.indir

    input_files = sorted(input_files)

    factors = compute_normalization_factors(input_files, args.num_markers, args.normalization)
    normalize_and_write(input_files, args.output, factors, args.num_markers)

    all_tile_metadata = []
    for tif_file_path in input_files:
        tile_metadata = extract_metadata(tif_file_path)
        all_tile_metadata.append(tile_metadata)

    create_ome(all_tile_metadata, args.output, args.num_markers)

def cli():
    """CLI entry point."""
    args = getOptions()
    warnings.filterwarnings("ignore", category=DeprecationWarning)
    main(args)

if __name__ == '__main__':
    cli()
