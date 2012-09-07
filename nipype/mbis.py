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



class MBISInputSpec( CommandLineInputSpec ):
    in_files = InputMultiPath( File(exists=True), copyfile=False,
                            desc='image, or multi-channel set of images, ' \
                                'to be segmented',
                            argstr='-C %s', position=-1, mandatory=True )

    mask = File(exists=True, desc='binary mask file', argstr='-x %s' )

    mask_auto = traits.Bool( desc='channels are implicitly masked', argstr='-M' )

    out_prefix = File(desc='base name of output files',
                        argstr='-o %s')  # uses in_file name as basename if none given

    number_classes = traits.Range(low=2, high=10, argstr='-n %d',
                                  desc='number of tissue-type classes', value=3)

    output_steps = traits.Bool(desc='output intermediate steps',
                                   argstr='--output-steps')


    output_biasfield = traits.Bool(desc='output estimated bias field',
                                   argstr='--bias-output')

    output_biascorrected = traits.Bool(desc='output restored image ' \
                                           '(bias-corrected image)',
                                       argstr='--bias-corrected-output')

    probability_maps = traits.Bool(desc='outputs a separate binary image for each ' \
                               'tissue type',
                           argstr='-g', value=True )

    priors = InputMultiPath(File(exist=True), desc='initialize with prior images',
                               argstr='-P %s', minlen=3, maxlen=10)

    no_bias = traits.Bool(desc='do not remove bias field',
                         argstr='--bias-skip', value=True )
    em_iters = traits.Range(low=1, high=50, value=3,
                                 desc='number of EM iterations',
                                 argstr='--em-iterations %d')
    mrf_iters = traits.Range(low=1, high=10,
                                 desc='number of MRF iterations',
                                 argstr='--mrf-iterations %d')

    mrf_lambda = traits.Range(low=0.01, high=1.0,
                         desc='MRF lambda parameter (segmentation spatial smoothness)',
                         argstr='-l %.3f')

    manual_init = File(exists=True, desc='Filename containing intensities',
                     argstr='-f %s')
    manual_init_means = File(exists=True, desc='Filename containing intensities',
                     argstr='--init-means %s')

#    no_pve = traits.Bool(desc='turn off PVE (partial volume estimation)',
#                        argstr='--nopve')
 
#    use_priors = traits.Bool(desc='use priors throughout',
#                             argstr='-P')   # must also set -a!,
#                                              # mutually inclusive??
#                                              # No, conditional
#                                              # mandatory... need to
#                                              # figure out how to
#                                              # handle with traits.
#   verbose = traits.Bool(desc='switch on diagnostic messages',
#                         argstr='-v')



class MBISOutputSpec( TraitedSpec ):
    out_classified = File(exists=True,
                            desc='path/name of binary segmented volume file' \
                            ' one val for each class  _mrf')
 
    bias_field = OutputMultiPath(File(desc='Estimated bias field _bias'))
    probability_maps = OutputMultiPath(File(desc='filenames, one for each class, for each ' \
                                'input, mrf_x'))
    restored_image = OutputMultiPath(File(desc='restored images (one for each input image) ' \
                              'named according to the input images _corrected_chXX'))


class MBIS(CommandLine):
    """ Use MBIS for segmentation and bias correction.

    Examples
    --------
    >>> from nipype.interfaces import mbis
    >>> from nipype.testing import example_data

    Assign options through the ``inputs`` attribute:

    >>> mbisr = mbis.MBIS()
    >>> mbisr.inputs.in_files = example_data('structural.nii')
    >>> out = mbisr.run() #doctest: +SKIP

    """
    _cmd = 'brain_seg'
    input_spec = MBISInputSpec
    output_spec = MBISOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()
        if not isdefined(self.inputs.number_classes):
            nclasses = 3
        else:
            nclasses = self.inputs.number_classes
        # when using multichannel, results basename is based on last
        # input filename
        if isdefined(self.inputs.out_prefix):
            basefile = self.inputs.out_prefix
        else:
            basefile = ''

        outputs['out_classified'] = self._gen_fname(basefile,
                                                      suffix='_mrf')
        if self.inputs.probability_maps:
            outputs['probability_maps'] = []
            for  i in range(nclasses):
                outputs['probability_maps'].append(
                        self._gen_fname(basefile, suffix='_mrf_seg%02d' % i))

        if isdefined(self.inputs.output_biascorrected):
            outputs['restored_image'] = []
            for val, f in enumerate(self.inputs.in_files):
                 # image numbering is 1-based
                 outputs['restored_image'].append(
                         self._gen_fname(basefile, suffix='_corrected_ch%02d' % (val + 1)))


        if self.inputs.output_biasfield:
            outputs['bias_field'] = []
            for val, f in enumerate(self.inputs.in_files):
                # image numbering is 1-based
                outputs['bias_field'].append(
                        self._gen_fname(basefile, suffix='_bias_ch%02d' % (val + 1)))
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
