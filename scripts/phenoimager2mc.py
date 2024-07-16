#!/usr/bin/env python
########################################################################################################
# AUTHOR: Chiara Schiller <chiara.schiller@uni-heidelberg.de>
# DESCRIPTION: Convert PhenoImager output .tif tiles to a single .tif file with correct OME-XML metadata
########################################################################################################

import os
import numpy as np

import argparse
import warnings
from argparse import RawDescriptionHelpFormatter

import tifffile as tiff
from PIL import Image
import xml.etree.ElementTree as ET
import ome_types
from ome_types.model import OME,Instrument,Pixels,TiffData,Channel,Plane,Pixels_DimensionOrder
from ome_types.model.simple_types import PixelsID
from ome_types.model.pixels import DimensionOrder
import copy
import platform
from uuid import uuid4
from ome_types import to_xml
from skimage import util


def getOptions(myopts=None):
    """ Function to pull in arguments """
    description = """ PhenoImager2MC """
    parser = argparse.ArgumentParser(
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    # Standard Input
    standard = parser.add_argument_group(
        title='Standard Inputs',
        description='Standard input for staging module.')
    # provide complete input folder containing subfolders per cycles with tif filed
    standard.add_argument(
        "-i",
        "--indir",
        dest="indir",
        action='store',
        required=True,
        help="Input folder with .tif files from PhenoImager.")
    standard.add_argument("-m",
                          "--num_markers",
                          dest="num_markers",
                          action="store",
                          required=True,
                          type=int,
                          help="Provide number of markers in the image.")

    # Tool Output
    output = parser.add_argument_group(title='Required output')
    output.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        action='store',
        required=True,
        help="Output folder, existing or will be newly created.")

    args = parser.parse_args(myopts)

    # Standardize paths
    args.indir = os.path.abspath(args.indir)
    return (args)


# concatentae tiles of all channels and tiles in one cycle
def concatenate_tiles(input_dir, output_file):
    # Stack the images with numpy
    """
    Stack all .tif tiles per cycle into a single .tif file

    :Arguments:
        :type input_dir: folder
        :param input_dir: Folder with all .tif input tiles in one cycle

        :type output_file: .tif file
        :param output_file: Stacked .tif file with all input tiles per cycle.

    :Returns:
        :rtype stacked_image: np.sarray
        :return scores_df: Stacked image saved as numpy array
    """
    # Make list of all files with a ".tif" extension in the input directory
    input_files = [
        os.path.join(input_dir, file) for file in os.listdir(input_dir)
        if file.endswith('.tif')
    ]
    # Sort the list of files alphabetically
    input_files = sorted(input_files)

    # Load the first image to determine the shape
    first_img = np.asarray(tiff.imread(input_files[0]))

    # Determine shape of the stacked image
    num_tiles = len(input_files)
    num_channels = first_img.shape[0]

    # Initialize an empty numpy array to hold all images
    # Shape: (Tiles, Channels, Z, Y, X)
    stack_shape = (num_tiles, num_channels, first_img.shape[1],
                   first_img.shape[2])
    stacked_image = np.empty(stack_shape, dtype=first_img.dtype)

    # Loop through each file and stack the images
    for idx, file in enumerate(input_files):
        img = np.asarray(tiff.imread(file))
        stacked_image[idx] = img

    # Save the stacked image as a TIFF file
    tiff.imsave(output_file, stacked_image)
    return stacked_image


# create function normalize_image
def normalize_image(output_file):
    """
    99th percentile normalization per channel

    :Arguments:
        :type output_file: np.sarray
        :param output_file: Stacked image from concatenate_tiles()

    :Returns:
        :rtype img_normalized: np.sarray
        :return img_normalized: Normalized image per channel
    """
    # Initialize an empty array for the normalized image
    # Read the float TIFF
    img = tiff.imread(output_file)

    img_normalized = np.empty_like(img, dtype=np.float32)

    # Normalize each channel separately
    for ch in range(img.shape[1]):  # Now channels are the second dimension
        percentile_99 = np.percentile(img[:, ch, :], 99)
        img_normalized[:, ch, :] = img[:, ch, :] / percentile_99
        img_normalized[:,
                       ch, :][img_normalized[:,
                                             ch, :] > 1] = 1  # Cap values at 1
    img_normalized = util.img_as_uint(img_normalized)
    # print the max value of the normalized image
    print("Max value of the normalized image: ", np.max(img_normalized))
    # Save the converted image
    tiff.imwrite(output_file, img_normalized)

    return img_normalized

# create function extract_metadata

def extract_metadata(tif_file_path):
    """
    Extract metadata from PhenoImager corrupted .tif files

    :Arguments:
        :type tif_file_path: file path
        :param tif_file_path: File path per tile in one cycle

    :Returns:
        :rtype metadata_list: list
        :return metadata_list: List of extracted metadata per tile across channels.
    """
    # List to store metadata for each page
    metadata_list = []

    # Open the TIFF file
    with tiff.TiffFile(tif_file_path) as tif:
        # Iterate over all pages (IFDs) in the TIFF file
        for i, page in enumerate(tif.pages):
            #if i >= 6:
            #    break
            # Dictionary to hold metadata for the current page
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
                #"pixel_unit": None,
                "device_name": None,
                "no_of_channels": None,
                "pixel_size" : 20046.2,
                #"physical_size_y" : 20046.2,
                "pixel_unit" : "PixelsPerCentimeter",
                "dpi" : None,
                "x_physical_size" : None,
                "y_physical_size" : None
            }
            #print(page.tags)

            #extract dpi information
            # Open the image file
            img = Image.open(tif_file_path)
            # Get detailed image info
            info = img.info
            page_metadata["dpi"] = info['dpi'][0]

            #calculate pixel size
            pixels_per_micrometer = page_metadata["dpi"] / 25400
            # Calculate the physical size in micrometers per pixel
            micrometers_per_pixel = 1 / pixels_per_micrometer
            #micrometers_per_pixel = inches_per_pixel * 25400

            # Extract position information
            x_position = page.tags.get('XPosition')
            y_position = page.tags.get('YPosition')
            if x_position and y_position:
                page_metadata["x_position"] = (x_position.value[0] / x_position.value[1])*10000
                page_metadata["y_position"] = (y_position.value[0] / y_position.value[1])*10000
                page_metadata["x_physical_size"] = micrometers_per_pixel
                page_metadata["y_physical_size"] = micrometers_per_pixel

            # Extract more information stored in xml format
            xml_string = page.tags.get('ImageDescription').value
            root = ET.fromstring(xml_string)

            # Extract <Name> tag
            name_element = root.find('./Name')
            if name_element is not None:
                page_metadata["marker"] = name_element.text

            # Extract <ExposureTime> tag
            exposure_time_element = root.find('./ExposureTime')
            if exposure_time_element is not None:
                page_metadata["exposure_time"] = exposure_time_element.text

            # Extract <SignalUnits> tag
            signal_units_element = root.find('./SignalUnits')
            if signal_units_element is not None:
                page_metadata["signal_unit"] = signal_units_element.text

            # Extract <PixelUnits> tag within <ScanProfile>
            scan_profile = root.find('.//ScanProfile')
            if scan_profile is not None:
                pixel_units_element = scan_profile.find('.//root')
                if pixel_units_element is not None:
                    pixel_units_element = scan_profile.find('.//PixelUnits')
                    if pixel_units_element is not None:
                        page_metadata["pixel_unit"] = pixel_units_element.text

            # Append the page metadata to the list
            metadata_list.append(page_metadata)

    return metadata_list


def create_ome(all_tile_metadata, input_dir, num_markers):
    """
    Creates and overwrites correct omexml metadata for the OME-TIFF files

    :Arguments:
        :type all_tile_metadata: list
        :param all_tile_metadata: List of extracted metadata per tile across channels from extract_metadata().

        :type input_dir: directory
        :param input_dir: input_dir of folder containg all tiles in a cycle.

        :type num_markers: integer
        :param num_markers: Number of markers in the image.

    """
    # define markers
    markers = []
    for i in range(0, num_markers):  #TODO do not hardcode here
        marker = all_tile_metadata[0][i]['marker']
        markers.append(marker)
    img_name = "test"
    #no_of_channels = len(all_tile_metadata[0])
    no_of_channels = args.num_markers
    bits_per_sample = all_tile_metadata[0][0]['bits_per_sample']
    no_of_tiles = len(all_tile_metadata)
    tiff.tiffcomment(input_dir, '')
    #--Generate tiff_data_blocks--#
    tiff_block = []
    #uuid_obj=UUID(file_name=img_name,value=uuid4().urn)
    for ch in range(0, no_of_channels):
        tiff_block.append(
            TiffData(
                first_c=ch,
                ifd=ch,
                plane_count=1  #,
                #uuid=uuid_obj
            ))
    #--Generate planes block (contains the information of each tile)--#
    plane_block = []
    #length_units=ome_types.model.simple_types.UnitsLength('µm')
    for ch in range(0, no_of_channels):
        plane_block.append(
            Plane(
                the_c=ch,
                the_t=0,
                the_z=0,
                position_x=0,  #x=0 is just a place holder
                position_y=0,  #y=0 is just a place holder
                position_z=0,
                exposure_time=0
                #position_x_unit=pixel_units,
                #position_y_unit=pixel_units
            ))
    #--Generate channels block--#
    chann_block = []
    #for ch in range(0,no_of_channels):
    for ch, chann_name in enumerate(markers):
        chann_block.append(
            Channel(
                id=ome_types.model.simple_types.ChannelID(
                    'Channel:{x}'.format(x=ch)),
                name=chann_name,
                color=ome_types.model.simple_types.Color((255, 255, 255)),
                emission_wavelength=1,  #place holder
                excitation_wavelength=1,  #place holder
            ))
    #--Generate pixels block--#
    pix_block = []
    ifd_counter = 0
    for t in range(0, no_of_tiles):
        template_plane_block = copy.deepcopy(plane_block)
        template_chann_block = copy.deepcopy(chann_block)
        template_tiffdata_block = copy.deepcopy(tiff_block)
        for ch, mark in enumerate(markers):
            #template_plane_block[ch].position_x=int(round(all_tile_metadata[t][0]['x_position']*10000))
            #template_plane_block[ch].position_y=int(round(all_tile_metadata[t][0]['y_position']*10000))
            #dont round them and try if ASHLAR can deal with it
            template_plane_block[ch].position_x = all_tile_metadata[t][0][
                'x_position']  #*10000
            template_plane_block[ch].position_y = all_tile_metadata[t][0][
                'y_position']  #*10000
            template_plane_block[ch].exposure_time = all_tile_metadata[t][ch][
                'exposure_time']
            template_chann_block[ch].id = 'Channel:{y}:{x}'.format(x=ch,
                                                                   y=100 +
                                                                   t)  ### why?
            template_chann_block[ch].name = mark
            template_tiffdata_block[ch].ifd = ifd_counter
            ifd_counter += 1
        pix_block.append(
            Pixels(
                id=ome_types.model.simple_types.PixelsID(
                    'Pixels:{x}'.format(x=t)),
                dimension_order=ome_types.model.pixels.DimensionOrder(
                    'XYCZT'),  ### check if the order is correct!!
                size_c=no_of_channels,
                size_t=1,
                size_x=all_tile_metadata[t][0]['image_width'],
                size_y=all_tile_metadata[t][0]['image_length'],
                size_z=1,
                type=ome_types.model.pixels.PixelType('uint32'),
                #type=bit_depth,
                big_endian=False,
                channels=template_chann_block,
                interleaved=False,
                physical_size_x=all_tile_metadata[t][0]['x_physical_size'],
                physical_size_y=all_tile_metadata[t][0]['y_physical_size'],
                #physical_size_x=pixel_size,
                #physical_size_x_unit=pixel_units,
                #physical_size_y=pixel_size,
                #physical_size_y_unit=pixel_units,
                physical_size_z=1.0,
                planes=template_plane_block,
                bits_per_sample=bits_per_sample,
                tiff_data_blocks=template_tiffdata_block))
    #--Generate image block--#
    img_block = []
    for t in range(0, no_of_tiles):
        img_block.append(
            ome_types.model.Image(id=ome_types.model.simple_types.ImageID(
                'Image:{x}'.format(x=t)),
                                  pixels=pix_block[t]))
    #--Create the OME object with all previously defined blocks--#
    ome_custom = OME()
    ome_custom.creator = " ".join([
        ome_types.__name__, ome_types.__version__, '/ python version-',
        platform.python_version()
    ])
    ome_custom.images = img_block
    ome_custom.uuid = uuid4().urn
    ome_xml = to_xml(ome_custom)
    tiff.tiffcomment(args.outdir, ome_xml)


# create ome tifs per cycle
def main(args):
    """
    Stack, normalize and add ome-metadata for each cycle to have one file per cycle with OME-XML metadata

    :Arguments:
        :type args.indir: directory
        :param args.indir: Input directory containing folders per cycle with .tif files.

        :type args.outdir: directory
        :param args.outdir: directory to save the output .tif files per cycle with OME-XML metadata. Will be created if not existent.
        
    """
    # Concatenate the tiles into a single stacked image
    concatenate_tiles(args.indir, args.outdir)
    # Normalize the image
    normalize_image(args.outdir)
    # Extract metadata from the TIFF files
    input_files = [
        os.path.join(args.indir, file) for file in os.listdir(args.indir)
        if file.endswith('.tif')
    ]
    # Sort the list of files alphabetically
    input_files = sorted(input_files)

    # Process all TIFF files in the input directory
    all_tile_metadata = []
    for tif_file_path in input_files:
        tile_metadata = extract_metadata(tif_file_path)
        all_tile_metadata.append(tile_metadata)

    # Create OME-XML metadata
    #if not os.path.exists(args.outdir):
    #    os.makedirs(args.outdir)

    #index = 1
    #base_name = os.path.basename(args.indir)
    #output_filename = args.outdir + base_name + "_cycle" + str(index) + ".tif"
    create_ome(all_tile_metadata, args.outdir, args.num_markers)

#  function to go one folder structure higher
# def main(args):
#     """
#     Loops ome_per_cycle function over all folders in the input directory

#     :Arguments:
#         :type args.indir: directory
#         :param args.indir: Input directory containing folders per cycle with .tif files.

#         :type args.outdir: directory
#         :param args.outdir: directory to save the output .tif files per cycle with OME-XML metadata. Will be created if not existent.

#     """
#loop the ome_per_cycle function per folder in the input path
# input_base = args.indir
# output_base = args.outdir

# # Create ome_output directory within the specified output path if it doesn't exist
# if not os.path.exists(output_base):
#     os.makedirs(output_base)

# index = 0
# for folder in sorted(os.listdir(input_base)):
#     subfolder_path = os.path.join(input_base, folder)

#     if os.path.isdir(subfolder_path):  # Check if it is a directory
#         # name output like foldername with _cycle1.tif as extension
#         index = index + 1
#         args.indir = subfolder_path
#         output_filename = folder + "_cycle" + str(index) + ".tif"
#         output_filename = output_filename.replace(' ',
#                                                   '_').replace('\t', '_')
#         args.outdir = os.path.join(output_base, output_filename)
#         print(args.outdir)
#         ome_per_cycle(args)


if __name__ == '__main__':
    """Tool is called on the command-line"""

    args = getOptions()
    warnings.filterwarnings("ignore", category=DeprecationWarning)

    main(args)
