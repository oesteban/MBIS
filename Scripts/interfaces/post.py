from nipype.interfaces.base import (InputMultiPath, OutputMultiPath, BaseInterface,
                                    BaseInterfaceInputSpec, TraitedSpec, traits, File)
from nipype.utils.filemanip import split_filename
import Scripts.util.postproc as pp

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
    _fnames = []

    def _run_interface(self, runtime):
        params = pp.loadParameters(self.inputs.parameters)

        _, base, _ = split_filename( self.inputs.in_maps[0] )

        self._fnames = pp.fusePV( self.inputs.in_files,
                               self.inputs.in_maps,
                               params,
                               self.inputs.distance,
                               self.inputs.reorder,
                               os.path.abspath( base ) )

        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        _, base, _ = split_filename( self.inputs.in_maps[0] )
        outputs["out_files"] = [ '%s_mseg%02d.nii.gz' % ( base, i ) for i in range(0,len(self.inputs.pure_tissues)) ]
        return outputs
