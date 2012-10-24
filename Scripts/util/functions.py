#!/usr/bin/python

def outlierFilter( in_file, quantile=99.998 ):
    import nibabel as nib
    import os
    from scipy.stats import scoreatpercentile
    img = nib.load( in_file )
    data = img.get_data()
    sample = data.reshape(-1)
    threshold = scoreatpercentile( sample, quantile )
    data[ data > threshold ] = 0.0
    nii = nib.Nifti1Image( data, img.get_affine(), img.get_header() )
    
    outfilename, ext = os.path.splitext( os.path.basename( in_file ) )
    if ext == '.gz':
        outfilename,ext = os.path.splitext( outfilename )
        ext+= '.gz'
    
    outfilename = os.path.abspath( os.path.join( os.getcwd(), '%s_filtered%s' % ( outfilename, ext ) ) )
    nib.save( nii, outfilename )
    return outfilename

def plot_slice(in_file, z_idx=5):
    import nibabel as nib
    import numpy as np
    import matplotlib.pyplot as plt
    # Load the image and collect the data
    # and orientation information
    img = nib.load(in_file)
    data = img.get_data()
    aff = img.get_affine()
    # Find the center of the brain matrix
    ctr = np.dot(np.linalg.inv(aff), [0, 0, 0, 1])[:3]
    # Plot the data
    #vmin, vmax = (0, 1) if data.dtype == np.int16 else (30, 150)
    vmin = np.min( data.reshape(-1) )
    vmax = np.max( data.reshape(-1) )
    plt.imshow(data[:, :, ctr[2] + z_idx], 
           cmap="gray", vmin=vmin, vmax=vmax)
    plt.xticks(()); plt.yticks(())

def measureVolume( in_files, icv_mask='' ):
    import nibabel as nib
    import os
    import numpy as np
    icv = 1.0
    
    if os.path.exists(icv_mask):
        mask = nib.load( icv_mask )
        icv = float( np.sum( mask.get_data().reshape(-1) ) )
    
    result = []
    for in_file in in_files:
        img = nib.load( in_file )
        volume = float( np.sum( img.get_data().reshape(-1) ) )
        result.append( (100.0 * volume) / icv )
    return result


