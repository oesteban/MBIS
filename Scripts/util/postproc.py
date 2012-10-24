
def loadParameters( fname ):
    import csv
    import numpy as np

    data = []
    with open( fname ) as csvfile:
        for row in csvfile:
            row = row.replace("[","").replace("]","").replace( " ", "")
            data.append( [ float(val) for val in row.split(',') ] )
    
    nparam = len( data[0] )
    nclasses = int( 0.5 * (nparam*4+1)**0.5 - 0.5 )
    # Load parameters
    mu = []
    sigma = []

    for row in data:
        mu.append( np.matrix( [ float(val) for val in row[0:nclasses] ] ) )
        sigma.append( np.matrix( np.reshape( [ float(val) for val in row[nclasses:] ], (nclasses,nclasses) ) ) )
    
    # Return arrays
    return ( mu, sigma )

def distancePV ( sample, mask, params_tissue1, params_tissue2, distance='euclidean' ):
    from scipy.spatial.distance import mahalanobis,euclidean
    import numpy as np

    # Direction vector between pure tissues
    d_vect = np.ravel(params_tissue2[0] - params_tissue1[0]).T
    mu1 = np.ravel(params_tissue1[0])
    mu2 = np.ravel(params_tissue2[0])
    SI1 = params_tissue1[1].getI()
    SI2 = params_tissue2[1].getI()

    if distance=='mahalanobis':
        norm = np.array( [ 1/(1+ mahalanobis(pix,mu2,SI2)/ mahalanobis(pix,mu1,SI1)) for pix in sample[mask==1] ] )
    elif distance=='dummy':
        norm = mask*0.5
    else:
        norm = np.array( [ 1/(1+ euclidean(pix,mu2)/ euclidean(pix,mu1)) for pix in sample[mask==1] ] )
    result = np.zeros( np.shape( mask ) )
    result[mask==1] = norm
    return result

def fusePV( in_files, in_maps, parameters, pt_list=[ 0, 2, 4 ], distance='mahalanobis', reorder=True, prefix='./' ):
    import nibabel as nib
    import numpy as np
    import os
    
    nmaps = len( in_maps )
    nclasses = np.shape( parameters )[1]
    assert nmaps == nclasses
    
    corder = range( 0, nclasses )
    npts = len( pt_list )
    pt_list = np.sort( pt_list )

    # Load probability maps
    initmaps = [ nib.load(f) for f in in_maps ]

    # If reorder is True, find unordered tissue signatures        
    if reorder:
        means = parameters[0]
        firstmeans = np.array( [ np.ravel(val)[0] for val in means ] )
        m = np.argsort(firstmeans)
        corder = np.take(corder, m)   
    new_idx = np.take( corder, pt_list )

    # Load images
    channels = [ nib.load(c) for c in in_files ]

    # Prepare sample
    data_shape = np.shape( channels[0] )
    data_sample = [ channel.get_data().reshape(-1) for channel in channels ]
    data_sample = np.swapaxes( data_sample, 0, 1)


    # Split between pv (partial volume) and pt (pure tissue) maps
    pt_niis = []
    pv_niis = []
    pt_param = []

    for tmap,i in zip(initmaps,range(0,nclasses)):
        idx = np.where( new_idx==i )[0]
        if len(idx) == 1:
            pt_niis.append( tmap )
            pt_param.append( [ parameters[0][i], parameters[1][i] ] )
        else:
            pv_niis.append( tmap )

    # Compute the steps required 
    steps = [ val-pt_list[i-1]-1 for val,i in zip( pt_list[1:], range(1,len(pt_list[1:])+1 ) ) ]

    # Extract data and initialize normalizer
    pt_maps = [ m.get_data().reshape(-1) for m in pt_niis ]
    pv_maps = [ m.get_data().reshape(-1) for m in pv_niis ]
    normalizer = np.zeros( np.shape( pt_maps[0] ) )

    # Process
    for pt_map,i in zip( pt_maps[:-1],range(0,len(pt_maps[:-1])) ):
        curr_steps = steps[i]
        for step in range(1,curr_steps+1):
            pv_idx = (step-1)+i
            pv_map = pv_maps[pv_idx]
            mask = np.zeros( np.shape( pt_map ) )
            mask[pv_map>0.001] = 1

            if not step == curr_steps:  # Direct addition of the pv map to the last pure tissue map
                pt_map+= pv_map
            else:                       # Split pv fraction proportionally to the distance to a contiguous pure tissue
                dist = distancePV( data_sample, mask, pt_param[i], pt_param[i+1], distance )
                pt_map+= dist * pv_map
                pt_maps[i+1]+= (1-dist) * pv_map

    # Compute normalizer
    for pt_map in pt_maps:
        normalizer[pt_map>0.001]+= pt_map[pt_map>0.001]
    
    # Generate output names
    fnames = [ '%s_mseg%02d.nii.gz' % ( prefix, i ) for i in range(0,len(pt_niis)) ]

    # Normalize and save
    for i in range(0,len(pt_niis)):
        pt_data = pt_maps[i]
        pt_data[normalizer>0]/= normalizer[normalizer>0]
        nii =  nib.Nifti1Image( np.reshape( pt_data, data_shape) , pt_niis[0].get_affine(), pt_niis[0].get_header() )
        nib.save( nii, fnames[i] )
    return fnames


#testpath= os.path.abspath( './' )
#in_files = [ '%s/%s.nii.gz' % (testpath, val) for val in ['T1', 'T2', 'PD' ] ]
#in_mask  = testpath + '/mask.nii.gz' 
#in_maps  = [ '%s/mbis_mrf_seg%02d.nii.gz' % ( testpath, val ) for val in range(1,6) ]
    

