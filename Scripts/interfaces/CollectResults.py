
from nipype.interfaces.base import (TraitedSpec, BaseInterface, traits,
                    BaseInterfaceInputSpec, File, isdefined)

class CollectResultsInputSpec(TraitedSpec):
    in_volume  = traits.ListFloat(desc='Input volume', mandatory=True)
    subject_id = traits.String(desc='Subject ID', mandatory=True)
    out_file = traits.File(desc='Output file', mandatory=True)

class CollectResultsOutputSpec(TraitedSpec):
    out_file = traits.File()

class CollectResults(BaseInterface):
    input_spec = CollectResultsInputSpec
    output_spec = CollectResultsOutputSpec

    def _run_interface(self, runtime):
        with open( self.inputs.out_file, 'a' ) as csvfile:
            writer = csv.writer(csvfile, delimiter=',',quoting=csv.QUOTE_NONE)
            writer.writerow( [ self.inputs.subject_id] + [ '%.3f' % el for el in self.inputs.in_volume] )
        csvfile.close()
        runtime.returncode=0
        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_file'] = self.inputs.out_file
        return outputs

