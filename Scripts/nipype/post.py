from nipype.interfaces.base import (InputMultiPath, OutputMultiPath, BaseInterface,
                                    BaseInterfaceInputSpec, traits, File)
from nipype.utils.filemanip import split_filename
from ..util import postproc as pp

import nibabel as nib
import numpy as np
import os

class PVProcessInputSpec(BaseInterfaceInputSpec):
    in_files = InputMultiPath( File(exists=True), copyfile=False, desc='image, or multi-channel set of images,' \
                                                                       'used for segmentation', mandatory=True )
    in_maps  = InputMultiPath( File(exists=True), copyfile=False, desc='Resulting tissue probability maps', mandatory=True )
    parameters = File( exists=True, desc='CSV file containing parameters for all the classes', mandatory=True)
    pure_tissues = traits.ListInt( value=[0,2,4], minlen=1, desc='Identifiers of the pure tissue classes' )
    dist = traits.String( value="euclidean", desc="Distance definition to be used" )
    reorder = traits.Bool( value=True, desc='Reorder classes if the classification is not ordered by first contrast means' )

class PVProcessOutputSpec(TraitedSpec):
    out_files = OutputMultiPath( File( desc='filenames, one for each pure tissue according to prefix_msegXX.nii.gz' ) )

class PVProcess(BaseInterface):
    input_spec = PVProcessInputSpec
    output_spec = PVProcessOutputSpec

    def _run_interface(self, runtime):
        params = pp.loadParameters(self.inputs.parameters)
        out_names = pp.fusePV( self.inputs.in_files,
                               self.inputs.in_maps,
                               params,:
                               self.inputs.distance,
                               self.inputs.reorder,
                               os.path.abspath( os.path.join( os.getcwd() ) )
        fname = self.inputs.volume
        img = nb.load(fname)
        data = np.array(img.get_data())

        active_map = data > self.inputs.threshold

        thresholded_map = np.zeros(data.shape)
        thresholded_map[active_map] = data[active_map]

        new_img = nb.Nifti1Image(thresholded_map, img.get_affine(), img.get_header())
        _, base, _ = split_filename(fname)
        nb.save(new_img, base + '_thresholded.nii')

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        fname = self.inputs.volume
        _, base, _ = split_filename(fname)
        outputs["thresholded_volume"] = os.path.abspath(base + '_thresholded.nii')
        return outputs
