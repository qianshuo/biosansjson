
import sys
import os
import time
import json
from mantid.simpleapi import mtd, logger
from drtsans.mono import biosans as sans
from drtsans.api import subtract_background
from drtsans.plots import plot_detector, plot_IQmod
from drtsans.sensitivity import load_sensitivity_workspace
from drtsans.instruments import extract_run_number
import drtsans.mask_utils as masking
from drtsans.path import registered_workspace
from drtsans.mono.biosans import convert_to_q
from drtsans.iq import determine_1d_linear_bins,determine_1d_log_bins
from drtsans.save_ascii import save_ascii_binned_1D, save_ascii_binned_2D
from drtsans.transmission import calculate_transmission
from biosans_common import *
start_time = time.time()
#json dictionary

if os.path.isfile(sys.argv[1]):
    print(sys.argv[1])
    with open(sys.argv[1], 'r') as fd:
        reduction_input = json.load(fd)
else:
    print("Usage: python biosansjson.py biosans.json")
    raise


#default mask
if reduction_input["configuration"]["useDefaultMask"]:
    default_mask=[{'Pixel':'1-18,239-256'}, {'Bank':'18-24,42-48'}]
else:
    default_mask = []

# chekcing if output directory exists, if it doesn't, creates the folder
output_dir = reduction_input["configuration"]["outputDir"]
for subfolder in ['1D','2D']:
    output_folder=os.path.join(output_dir,subfolder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)



# load all required files at the beginning, and transform them to histograms
prefix = ""
instrument_name = reduction_input["instrumentName"]
ipts = reduction_input["iptsNumber"]
sample = reduction_input["runNumber"]
sample_trans = reduction_input["transmission"]["runNumber"]
bkgd = reduction_input["background"]["runNumber"]
bkgd_trans = reduction_input["background"]["transmission"]["runNumber"]
empty = reduction_input["emptyTrans"]["runNumber"]
center = reduction_input["beamCenter"]["runNumber"]

if reduction_input["configuration"]["useBlockedBeam"]:
    blocked_beam = reduction_input["configuration"]["BlockBeamFileName"]
else:
    blocked_beam = None

load_params = {}
wavelength = reduction_input["configuration"]["wavelength"]
wavelengthSpread = reduction_input["configuration"]["wavelengthSpread"]
if wavelength and wavelengthSpread:
    load_params["wavelengthinstrumentName"] = wavelength
    load_params["wavelengthSpread"] = wavelengthSpread

for run_number in [center, sample, bkgd, empty, sample_trans, bkgd_trans, blocked_beam]:
    if run_number:
        ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
        if not registered_workspace(ws_name):
            filename = f"/HFIR/{instrument_name}/IPTS-{ipts}/nexus/{instrument_name}_{run_number}.nxs.h5"
            print(f"Loading filename {filename}")
            sans.load_events_and_histogram(filename, output_workspace=ws_name, **load_params)
            for btp_params in default_mask:
                sans.api.apply_mask(ws_name, **btp_params)

raw_sample_ws = mtd[f'{prefix}_{instrument_name}_{sample}_raw_histo']
raw_bkgd_ws = mtd[f'{prefix}_{instrument_name}_{bkgd}_raw_histo'] if bkgd else None
raw_blocked_ws = mtd[f'{prefix}_{instrument_name}_{blocked_beam}_raw_histo'] if blocked_beam else None

# do the same for dark current if exists
dark_current_main = None
dark_current_wing = None
if reduction_input["configuration"]["useDarkFileName"]:
    dark_current_file_main = reduction_input["configuration"]["darkMainFileName"]
    dark_current_file_wing = reduction_input["configuration"]["darkWingFileName"]
    if dark_current_file_main and dark_current_file_wing:
        run_number = extract_run_number(dark_current_file_main)
        ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
        if not registered_workspace(ws_name):
            print(f"Loading filename {dark_current_file_main}")
            dark_current_main = sans.load_events_and_histogram(dark_current_file_main, output_workspace=ws_name,
                                                               **load_params)
            for btp_params in default_mask:
                sans.api.apply_mask(ws_name, **btp_params)
        run_number = extract_run_number(dark_current_file_wing)
        ws_name = f'{prefix}_{instrument_name}_{run_number}_raw_histo'
        if not registered_workspace(ws_name):
            print(f"Loading filename {dark_current_file_wing}")
            dark_current_wing = sans.load_events_and_histogram(dark_current_file_wing, output_workspace=ws_name,
                                                               **load_params)
            for btp_params in default_mask:
                sans.api.apply_mask(ws_name, **btp_params)

# load required processed_files
sensitivity_main_ws_name = None
sensitivity_wing_ws_name = None
if reduction_input["configuration"]["useSensitivityFileName"]:
    flood_file_main = reduction_input["configuration"]["sensitivityMainFileName"]
    flood_file_wing = reduction_input["configuration"]["sensitivityWingFileName"]
    if flood_file_main and flood_file_wing:
        sensitivity_main_ws_name = f'{prefix}_main_sensitivity'
        sensitivity_wing_ws_name = f'{prefix}_wing_sensitivity'
        if not registered_workspace(sensitivity_main_ws_name):
            print(f"Loading filename {flood_file_main}")
            load_sensitivity_workspace(flood_file_main, output_workspace=sensitivity_main_ws_name)
        if not registered_workspace(sensitivity_wing_ws_name):
            print(f"Loading filename {flood_file_wing}")
            load_sensitivity_workspace(flood_file_wing, output_workspace=sensitivity_wing_ws_name)

mask_ws = None
if reduction_input["configuration"]["useMaskFileName"]:
    custom_mask_file = reduction_input["configuration"]["maskFileName"]
    if custom_mask_file:
        mask_ws_name = f'{prefix}_mask'
        if not registered_workspace(mask_ws_name):
            print(f"Loading filename {custom_mask_file}")
            mask_ws = masking.load_mask(custom_mask_file, output_workspace=mask_ws_name)
print('Done loading')

#setting up for different absolute scaling method
flux_method = reduction_input["configuration"]["normalization"]
if len(reduction_input["configuration"]["mmRadiusForTransmission"]) == 0:
    transmission_radius = None
else:
    transmission_radius = reduction_input["configuration"]["mmRadiusForTransmission"]
solid_angle = reduction_input["configuration"]["useSolidAngleCorrection"]
sample_trans_value = reduction_input["transmission"]["value"]
bkg_trans_value = reduction_input["background"]["transmission"]["value"]
theta_deppendent_transmission = reduction_input["configuration"]["useThetaDepTransCorrection"]
mask_panel = None #reduction_input["configuration"]["useMaskBackTubes"]
output_suffix = ''
try:
    thickness = float(reduction_input['thickness'])
except ValueError:
    thickness=1.
absolute_scale_method=reduction_input["configuration"]["absoluteScaleMethod"]
try:
    beam_radius=float(reduction_input["configuration"]["DBScalingBeamRadius"])
except ValueError:
    beam_radius=None
try:
    absolute_scale = float(reduction_input["configuration"]["StandardAbsoluteScale"])
except ValueError:
    absolute_scale = 1.

# find the center
ws_center = mtd[f'{prefix}_{instrument_name}_{center}_raw_histo']
xc, yc, yw = sans.find_beam_center(ws_center)
print("Center  =", xc, yc, yw)

#manual beam center
xc = -0.01320
yc = -0.01440
yw = -0.00036

####plot_detector(ws_center,None,backend='mpl')




############################after defining a few functions in biosans common, start more reduction prep
#PREPRARE for Transmission
# empty beam transmission workspace
if empty:
    raw_ws_name = f'{prefix}_{instrument_name}_{empty}_raw_histo'
    empty_trans_ws_name = f'{prefix}_empty'
    empty_trans_ws = prepare_data_workspaces(raw_ws_name,
                                             flux_method=flux_method,
                                             center_x=xc,
                                             center_y=yc,
                                             center_y_wing=yw,
                                             solid_angle=False,
                                             sensitivity_ws=sensitivity_main_ws_name,
                                             output_workspace=empty_trans_ws_name)
else:
    empty_trans_ws = None

# background transmission
if bkgd_trans:
    raw_ws_name = f'{prefix}_{instrument_name}_{bkgd_trans}_raw_histo'
    bkgd_trans_ws_name = f'{prefix}_bkgd_trans'
    bkgd_trans_ws_processed = prepare_data_workspaces(raw_ws_name,
                                                      flux_method=flux_method,
                                                      center_x=xc,
                                                      center_y=yc,
                                                      center_y_wing=yw,
                                                      solid_angle=False,
                                                      sensitivity_ws=sensitivity_main_ws_name,
                                                      output_workspace=bkgd_trans_ws_name)
    bkgd_trans_ws = calculate_transmission(bkgd_trans_ws_processed, empty_trans_ws,
                                           radius=transmission_radius,radius_unit="mm")
    print ('Background transmission =',bkgd_trans_ws.extractY()[0,0])
else:
    bkgd_trans_ws = None

# sample transmission
if sample_trans:
    raw_ws_name = f'{prefix}_{instrument_name}_{sample_trans}_raw_histo'
    sample_trans_ws_name = f'{prefix}_sample_trans'
    sample_trans_ws_processed = prepare_data_workspaces(raw_ws_name,
                                                        flux_method=flux_method,
                                                        center_x=xc,
                                                        center_y=yc,
                                                        center_y_wing=yw,
                                                        solid_angle=False,
                                                        sensitivity_ws=sensitivity_main_ws_name,
                                                        output_workspace=sample_trans_ws_name)
    sample_trans_ws = calculate_transmission(sample_trans_ws_processed, empty_trans_ws,
                                             radius=transmission_radius,radius_unit="mm")
    print ('Sample transmission =',sample_trans_ws.extractY()[0,0])
else:
    sample_trans_ws = None
# Real reduction starts here
################################################
processed_data_main = process_single_configuration(raw_sample_ws,
                                                   sample_trans_ws=sample_trans_ws,
                                                   sample_trans_value=sample_trans_value,
                                                   bkg_ws_raw=raw_bkgd_ws,
                                                   bkg_trans_ws=bkgd_trans_ws,
                                                   bkg_trans_value=bkg_trans_value,
                                                   blocked_ws_raw=raw_blocked_ws,
                                                   theta_deppendent_transmission=theta_deppendent_transmission,
                                                   center_x=xc, center_y=yc, center_y_wing=yw,
                                                   dark_current=dark_current_main,
                                                   flux_method=flux_method,    # normalization (time/monitor)
                                                   mask_detector='wing_detector',  # main or wing
                                                   mask_ws=mask_ws,     # apply a custom mask from workspace
                                                   mask_panel=mask_panel,     # mask back or front panel
                                                   mask_btp=None, #{'Bank':'23,24,47,48'},     # mask bank/tube/pixel
                                                   solid_angle=True,
                                                   sensitivity_workspace=mtd[sensitivity_main_ws_name],
                                                   output_workspace='processed_data_main',
                                                   output_suffix=output_suffix,
                                                   thickness=thickness,
                                                   absolute_scale_method=absolute_scale_method,
                                                   empty_beam_ws=empty_trans_ws,
                                                   beam_radius=beam_radius,
                                                   absolute_scale=absolute_scale,
                                                   keep_processed_workspaces=False)
#plot_detector(processed_data_main,None,'mpl')


processed_data_wing = process_single_configuration(raw_sample_ws,
                                                   sample_trans_ws=sample_trans_ws,
                                                   sample_trans_value=sample_trans_value,
                                                   bkg_ws_raw=raw_bkgd_ws,
                                                   bkg_trans_ws=bkgd_trans_ws,
                                                   bkg_trans_value=bkg_trans_value,
                                                   blocked_ws_raw=raw_blocked_ws,
                                                   theta_deppendent_transmission=theta_deppendent_transmission,
                                                   center_x=xc, center_y=yc, center_y_wing=yw,
                                                   dark_current=dark_current_wing,
                                                   flux_method=flux_method,    # normalization (time/monitor)
                                                   mask_detector='detector1',  # main or wing
                                                   mask_ws=mask_ws,     # apply a custom mask from workspace
                                                   mask_panel=mask_panel,     # mask back or front panel
                                                   mask_btp=None,     # mask bank/tube/pixel
                                                   solid_angle=True,
                                                   sensitivity_workspace=mtd[sensitivity_wing_ws_name],
                                                   output_workspace='processed_data_wing',
                                                   output_suffix=output_suffix,
                                                   thickness=thickness,
                                                   absolute_scale_method=absolute_scale_method,
                                                   empty_beam_ws=empty_trans_ws,
                                                   beam_radius=beam_radius,
                                                   absolute_scale=absolute_scale,
                                                   keep_processed_workspaces=False)
#plot_detector(processed_data_wing,None,'mpl')

#Binning setup #################################################
nxbins_main = reduction_input["configuration"]["numMainQxQyBins"]
nybins_main = nxbins_main
nxbins_wing = reduction_input["configuration"]["numWingQxQyBins"]
nybins_wing = nxbins_wing
bin1d_type = reduction_input["configuration"]["1DQbinType"]
linear_binning = reduction_input["configuration"]["QbinType"] != 'log'
even_decades = reduction_input["configuration"]["EvenDecades"]
nbins_main = reduction_input["configuration"]["numMainQBins"]
nbins_wing = reduction_input["configuration"]["numWingQBins"]
outputFilename = reduction_input["outputFilename"]

import drtsans
import matplotlib.pyplot as plt
import numpy as np
from mantid.simpleapi import MaskBTP
from drtsans.iq import BinningMethod, determine_1d_linear_bins, determine_1d_log_bins

if reduction_input["configuration"]["useErrorWeighting"]:
    bin_method = BinningMethod.WEIGHTED
else:
    bin_method = BinningMethod.NOWEIGHT

# Q 2D main
q_data = sans.convert_to_q(processed_data_main, mode='azimuthal')
qx_min = np.min(q_data.qx)
qx_max = np.max(q_data.qx)
binning_x = determine_1d_linear_bins(qx_min, qx_max, nxbins_main)
qy_min = np.min(q_data.qy)
qy_max = np.max(q_data.qy)
binning_y = determine_1d_linear_bins(qy_min, qy_max, nybins_main)

iq_output = drtsans.iq.bin_intensity_into_q2d(q_data,
                                       binning_x,
                                       binning_y,
                                       method=bin_method)

filename = os.path.join(output_dir,'2D',f'{outputFilename}_2D_main.txt')
save_ascii_binned_2D(filename, "I(Qx,Qy)", iq_output)
fig, ax = plt.subplots()
pcm = ax.pcolormesh(iq_output.qx, iq_output.qy, iq_output.intensity.T,
                    cmap='jet')
ax.set_title("Main")
ax.set_xlabel("$Q_x (\AA^{-1})$")
ax.set_ylabel("$Q_y (\AA^{-1})$")
fig.colorbar(pcm, ax=ax)

filename = os.path.join(output_dir,'2D',f'{outputFilename}_2D_main.png')
fig.savefig(filename)

#2D wing
q_data = sans.convert_to_q(processed_data_wing, mode='azimuthal')
qx_min = np.min(q_data.qx)
qx_max = np.max(q_data.qx)
binning_x = determine_1d_linear_bins(qx_min, qx_max, nxbins_wing)
qy_min = np.min(q_data.qy)
qy_max = np.max(q_data.qy)
binning_y = determine_1d_linear_bins(qy_min, qy_max, nybins_wing)

iq_output = drtsans.iq.bin_intensity_into_q2d(q_data,
                                       binning_x,
                                       binning_y,
                                       method=bin_method)

filename = os.path.join(output_dir,'2D',f'{outputFilename}_2D_wing.txt')
save_ascii_binned_2D(filename, "I(Qx,Qy)", iq_output)
fig, ax = plt.subplots()
pcm = ax.pcolormesh(iq_output.qx, iq_output.qy, iq_output.intensity.T,
                    cmap='jet')
ax.set_title("Wing")
ax.set_xlabel("$Q_x (\AA^{-1})$")
ax.set_ylabel("$Q_y (\AA^{-1})$")
fig.colorbar(pcm, ax=ax)

filename = os.path.join(output_dir,'2D',f'{outputFilename}_2D_wing.png')
fig.savefig(filename)


# 1D|Q|
from mantid.simpleapi import MaskAngle

MaskAngle(processed_data_wing, MinAngle=57)
MaskAngle(processed_data_main, MaxAngle=0.165)

q_data_main = sans.convert_to_q(processed_data_main, mode='scalar')
q_min_main = q_data_main.mod_q.min()
q_max_main = q_data_main.mod_q.max()
q_min_main = 0.003
q_max_main = 0.045
print('Qmin1d',q_min_main)
q_data_wing = sans.convert_to_q(processed_data_wing, mode='scalar')
q_min_wing = q_data_wing.mod_q.min()
q_max_wing = q_data_wing.mod_q.max()
q_min_wing = 0.035
q_max_wing = 0.850

if bin1d_type == 'scalar':
    if linear_binning:
        linear_bins_main = determine_1d_linear_bins(q_min_main, q_max_main, nbins_main)
        iq_output_main = drtsans.iq.bin_intensity_into_q1d(q_data_main, linear_bins_main, bin_method=bin_method)
        linear_bins_wing = determine_1d_linear_bins(q_min_wing, q_max_wing, nbins_wing)
        iq_output_wing = drtsans.iq.bin_intensity_into_q1d(q_data_wing, linear_bins_wing, bin_method=bin_method)
    else:
        log_bins_main = determine_1d_log_bins(q_min_main, q_max_main, nbins_main, even_decade=even_decades)
        iq_output_main = drtsans.iq.bin_intensity_into_q1d(q_data_main, log_bins_main, bin_method=bin_method)
        log_bins_wing = determine_1d_log_bins(q_min_wing, q_max_wing, nbins_wing, even_decade=even_decades)
        iq_output_wing = drtsans.iq.bin_intensity_into_q1d(q_data_wing, log_bins_wing, bin_method=bin_method)
else:
    raise NotImplementedError

ascii_1D_filename = os.path.join(output_dir, '1D', f'{outputFilename}_1D_main.txt')
save_ascii_binned_1D(ascii_1D_filename, "I(Q)", iq_output_main)
ascii_1D_filename = os.path.join(output_dir, '1D', f'{outputFilename}_1D_wing.txt')
save_ascii_binned_1D(ascii_1D_filename, "I(Q)", iq_output_wing)

OLT_Qmin = 0.03 #q_min_wing
OLT_Qmax = 0.05 #q_max_main
iq_output_both = sans.stitch_profiles(profiles=[iq_output_main, iq_output_wing],
                                      overlaps=[OLT_Qmin, OLT_Qmax],
                                      target_profile_index=0)

ascii_1D_filename = os.path.join(output_dir, '1D', f'{outputFilename}_1D_both.txt')
save_ascii_binned_1D(ascii_1D_filename, "I(Q)", iq_output_both)

png_file = os.path.join(output_dir, '1D', f'{outputFilename}_1D.png')
plot_IQmod([iq_output_main, iq_output_wing, iq_output_both], png_file, loglog=True, backend='mpl')


#plt.close(fig)


#Shuo: print some info
end_time = time.time()
print(end_time-start_time)
localtime = time.asctime( time.localtime(time.time()) )
print("Local current time :", localtime)
print("beam center", xc, yc, yw)
print ('Sample transmission =',sample_trans_ws.extractY()[0,0])
print ('Background transmission =',bkgd_trans_ws.extractY()[0,0])
print("DONEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE")
