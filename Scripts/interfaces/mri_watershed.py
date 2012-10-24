# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
""" MBIS nipype interface definition
    Based on FAST interface definition

    Change directory to provide relative paths for doctests
    >>> import os
    >>> filepath = os.path.dirname( os.path.realpath( __file__ ) )
    >>> datadir = os.path.realpath(os.path.join(filepath, '../../testing/data'))
    >>> os.chdir(datadir)
"""

import os, os.path as op
import warnings
import numpy as np

from nipype.interfaces.base import (TraitedSpec, File, InputMultiPath,
                                    OutputMultiPath, Undefined, traits,
                                    isdefined, OutputMultiPath, 
                                    CommandLineInputSpec, CommandLine )
from nipype.utils.filemanip import split_filename,fname_presuffix

from nibabel import load


warn = warnings.warn
warnings.filterwarnings('always', category=UserWarning)



class WatershedInputSpec( CommandLineInputSpec ):
    in_file = File(exists=True, desc='T1 image to be segmented', argstr='%s', position=-2, mandatory=True )
    out_file = File(exists=True, desc='Segmented brain', argstr='%s', position=-1,genfile=True )

class WatershedOutputSpec( TraitedSpec ):
    out_file = File(desc='Segmented brain',mandatory=True )


class Watershed(CommandLine):
    _cmd = 'mri_watershed'
    input_spec = WatershedInputSpec
    output_spec = WatershedOutputSpec

    def _gen_outfilename(self):
        out_file = self.inputs.out_file
        if not isdefined(out_file) and isdefined(self.inputs.in_file):
            out_file = self._gen_fname(self.inputs.in_file,
                                       suffix='_brain')
        return os.path.abspath(out_file)

    def _list_outputs(self):
         outputs = self.output_spec().get()
         outputs['out_file'] = self._gen_outfilename()
         return outputs
 

    def _gen_fname(self, prefix, suffix=None, ext='.nii.gz' ):
        """Generate a filename based on the given parameters.

        The filename will take the form: preffix<suffix><ext>.

        Parameters
        ----------
        prefix : str
            Filename to base the new filename on.

        suffix : str
            Suffix to add to the `basename`.  (defaults is '' )

        ext : str
            Desired extension (default is nii.gz)

        Returns
        -------
        fname : str
            New filename based on given parameters.

        """
        if ((prefix == '') or (prefix is None) ):
            prefix = './'
        if suffix is None:
            suffix = ''
        suffix = ''.join((suffix,ext))
        fname = fname_presuffix(prefix, suffix=suffix, use_ext=False)
        return fname

    def _gen_filename(self, name):
        if name == 'out_file':
            return self._gen_outfilename()
        return None

