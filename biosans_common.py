import os
from mantid.simpleapi import mtd, MaskDetectors

from drtsans.solid_angle import solid_angle_correction
from drtsans.mono.biosans.beam_finder import center_detector
from drtsans.mask_utils import apply_mask
from drtsans.mono.load import load_events, transform_to_wavelength
from drtsans.mono.normalization import normalize_by_monitor, normalize_by_time
from drtsans.mono.dark_current import subtract_dark_current
from drtsans.sensitivity import apply_sensitivity_correction, load_sensitivity_workspace
from drtsans.transmission import apply_transmission_correction
from drtsans.thickness_normalization import normalize_by_thickness
from drtsans.mono.absolute_units import empty_beam_scaling
from drtsans.mono.gpsans.attenuation import attenuation_factor
from drtsans.path import registered_workspace
from drtsans import subtract_background


def prepare_data_workspaces(data,
                            center_x=None, center_y=None, center_y_wing=None,
                            dark_current=None,
                            flux_method=None,  # normalization (time/monitor)
                            mask_ws=None,  # apply a custom mask from workspace
                            mask_detector=None,  # main or wing
                            mask_panel=None,  # mask back or front panel
                            mask_btp=None,  # mask bank/tube/pixel
                            solid_angle=True,
                            sensitivity_workspace=None,
                            output_workspace=None,
                            output_suffix='', **kwargs):
    r"""
    Given a " raw"data workspace, this function provides the following:

        - centers the detector
        - subtracts dark current
        - normalize by time or monitor
        - applies masks
        - corrects for solid angle
        - corrects for sensitivity

    All steps are optional. data, mask_ws, dark_current are either None
    or histogram workspaces. This function does not load any file.

    Parameters
    ----------
    data: ~mantid.dataobjects.Workspace2D
        raw workspace (histogram)
    center_x: float
        Move the center of the detector to this X-coordinate. If :py:obj:`None`, the
        detector will be moved such that the X-coordinate of the intersection
        point between the neutron beam and the detector array will have ``x=0``.
    center_y: float
        Move the center of the detector to this Y-coordinate. If :py:obj:`None`, the
        detector will be moved such that the Y-coordinate of the intersection
        point between the neutron beam and the detector array will have ``y=0``.
    dark_current: ~mantid.dataobjects.Workspace2D
        histogram workspace containing the dark current measurement
    flux_method: str
        Method for flux normalization. Either 'monitor', or 'time'.
    mask_ws: ~mantid.dataobjects.Workspace2D
        Mask workspace
    mask_panel: str
        Either 'front' or 'back' to mask whole front or back panel.
    mask_btp: dict
        Additional properties to Mantid's MaskBTP algorithm
    solid_angle: bool
        Apply the solid angle correction
    sensitivity_workspace: str, ~mantid.api.MatrixWorkspace
        workspace containing previously calculated sensitivity correction. This
        overrides the sensitivity_filename if both are provided.
    output_workspace: str
        The output workspace name. If None will create data.name()+output_suffix
    output_suffix: str
        replace '_raw_histo' in the output workspace name.
        If empty, the default is '_processed_histo'

    Returns
    -------
    ~mantid.dataobjects.Workspace2D
        Reference to the processed workspace
    """
    if not output_workspace:
        output_workspace = str(data)
        output_workspace.replace('_raw_histo', '') + '_processed_histo'

    mtd[str(data)].clone(OutputWorkspace=output_workspace)  # name gets into workspace

    if center_x is not None and center_y is not None and center_y_wing is not None:
        center_detector(output_workspace, center_x=center_x, center_y=center_y, center_y_wing=center_y_wing)

    # Dark current
    if dark_current is not None:
        subtract_dark_current(output_workspace, dark_current)

    # Normalization
    if str(flux_method).lower() == 'monitor':
        normalize_by_monitor(output_workspace)
    elif str(flux_method).lower() == 'time':
        normalize_by_time(output_workspace)

    # Mask either detector
    if mask_detector is not None:
        MaskDetectors(output_workspace, ComponentList=mask_detector)

    # Additional masks
    if mask_btp is None:
        mask_btp = dict()
    apply_mask(output_workspace, panel=mask_panel, mask=mask_ws, **mask_btp)

    # Solid angle
    if solid_angle:
        solid_angle_correction(output_workspace)

    # Sensitivity
    if sensitivity_workspace is not None:
        apply_sensitivity_correction(output_workspace,
                                     sensitivity_workspace=sensitivity_workspace)

    return mtd[output_workspace]


def process_single_configuration(sample_ws_raw,
                                 sample_trans_ws=None,
                                 sample_trans_value=None,
                                 bkg_ws_raw=None,
                                 bkg_trans_ws=None,
                                 bkg_trans_value=None,
                                 blocked_ws_raw=None,
                                 theta_deppendent_transmission=True,
                                 center_x=None,
                                 center_y=None,
                                 dark_current=None,
                                 flux_method=None,  # normalization (time/monitor)
                                 mask_ws=None,  # apply a custom mask from workspace
                                 mask_detector=None,
                                 mask_panel=None,  # mask back or front panel
                                 mask_btp=None,  # mask bank/tube/pixel
                                 solid_angle=True,
                                 sensitivity_workspace=None,
                                 output_workspace=None,
                                 output_suffix='',
                                 thickness=1.,
                                 absolute_scale_method='standard',
                                 empty_beam_ws=None,
                                 beam_radius=None,
                                 absolute_scale=1.,
                                 keep_processed_workspaces=True,
                                 **kwargs):
    r"""
    This function provides full data processing for a single experimental configuration,
    starting from workspaces (no data loading is happening inside this function)

    Parameters
    ----------
    sample_ws_raw: ~mantid.dataobjects.Workspace2D
        raw data histogram workspace
    sample_trans_ws: ~mantid.dataobjects.Workspace2D
        optional histogram workspace for sample transmission
    sample_trans_value: float
        optional value for sample transmission
    bkg_ws_raw: ~mantid.dataobjects.Workspace2D
        optional raw histogram workspace for background
    bkg_trans_ws: ~mantid.dataobjects.Workspace2D
        optional histogram workspace for background transmission
    bkg_trans_value: float
        optional value for background transmission
    blocked_ws_raw: ~mantid.dataobjects.Workspace2D
        optional histogram workspace for blocked beam
    theta_deppendent_transmission: bool
        flag to apply angle dependent transmission
    center_x: float
        x center for the beam
    center_y: float
        y center for the beam
    dark_current: ~mantid.dataobjects.Workspace2D
        dark current workspace
    flux_method: str
        normalization by time or monitor
    mask_ws: ~mantid.dataobjects.Workspace2D
        user defined mask
    mask_panel: str
        mask fron or back panel
    mask_btp: dict
        optional bank, tube, pixel to mask
    solid_angle: bool
        flag to apply solid angle
    sensitivity_workspace: ~mantid.dataobjects.Workspace2D
        workspace containing sensitivity
    output_workspace: str
        output workspace name
    output_suffix:str
        suffix for output workspace
    thickness: float
        sample thickness (cm)
    absolute_scale_method: str
        method to do absolute scaling (standard or direct_beam)
    empty_beam_ws: ~mantid.dataobjects.Workspace2D
        empty beam workspace for absolute scaling
    beam_radius: float
        beam radius for absolute scaling
    absolute_scale: float
        absolute scaling value for standard method
    keep_processed_workspaces: bool
        flag to keep the processed blocked beam and background workspaces

    Returns
    -------
    ~mantid.dataobjects.Workspace2D
        Reference to the processed workspace
    """
    if not output_workspace:
        output_workspace = output_suffix + '_sample'

    # create a common configuration for prepare data
    prepare_data_conf = {'center_x': center_x,
                         'center_y': center_y,
                         'dark_current': dark_current,
                         'flux_method': flux_method,
                         'mask_ws': mask_ws,
                         'mask_detector': mask_detector,
                         'mask_panel': mask_panel,
                         'mask_btp': mask_btp,
                         'solid_angle': solid_angle,
                         'sensitivity_workspace': sensitivity_workspace}

    # process blocked
    if blocked_ws_raw:
        blocked_ws_name = output_suffix + '_blocked'
        if not registered_workspace(blocked_ws_name):
            blocked_ws = prepare_data_workspaces(blocked_ws_raw,
                                                 output_workspace=blocked_ws_name,
                                                 **prepare_data_conf)
        else:
            blocked_ws = mtd[blocked_ws_name]

    # process sample
    sample_ws = prepare_data_workspaces(sample_ws_raw,
                                        output_workspace=output_workspace,
                                        **prepare_data_conf)
    if blocked_ws_raw:
        sample_ws = subtract_background(sample_ws, blocked_ws)
    # apply transmission to the sample
    if sample_trans_ws or sample_trans_value:
        sample_ws = apply_transmission_correction(sample_ws,
                                                  trans_workspace=sample_trans_ws,
                                                  trans_value=sample_trans_value,
                                                  theta_dependent=theta_deppendent_transmission,
                                                  output_workspace=output_workspace)

    # process background, if not already processed
    if bkg_ws_raw:
        bkgd_ws_name = output_suffix + '_background'
        if not registered_workspace(bkgd_ws_name):
            bkgd_ws = prepare_data_workspaces(bkg_ws_raw,
                                              output_workspace=bkgd_ws_name,
                                              **prepare_data_conf)
            if blocked_ws_raw:
                bkgd_ws = subtract_background(bkgd_ws, blocked_ws)
            # apply transmission to bkgd
            if bkg_trans_ws or bkg_trans_value:
                bkgd_ws = apply_transmission_correction(bkgd_ws,
                                                        trans_workspace=bkg_trans_ws,
                                                        trans_value=bkg_trans_value,
                                                        theta_dependent=theta_deppendent_transmission,
                                                        output_workspace=bkgd_ws_name)
        else:
            bkgd_ws = mtd[bkgd_ws_name]
        # subtract background
        sample_ws = subtract_background(sample_ws, bkgd_ws, output_workspace=output_workspace)

        if not keep_processed_workspaces:
            bkgd_ws.delete()

    if blocked_ws_raw and not keep_processed_workspaces:
        blocked_ws.delete()

    # finalize with absolute scale and thickness
    sample_ws = normalize_by_thickness(sample_ws, thickness)

    # standard method assumes absolute scale from outside
    if absolute_scale_method == 'direct_beam':
        try:
            empty = mtd[str(empty_beam_ws)]
        except KeyError:
            raise ValueError(f"Could not find empty beam {str(empty_beam_ws)}")

        ac, ace = attenuation_factor(empty)
        sample_ws = empty_beam_scaling(sample_ws,
                                       empty,
                                       beam_radius=beam_radius,
                                       unit='mm',
                                       attenuator_coefficient=ac,
                                       attenuator_error=ace,
                                       output_workspace=output_workspace)
    else:
        sample_ws *= absolute_scale

    return mtd[str(sample_ws)]
