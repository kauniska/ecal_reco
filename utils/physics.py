def dr(hit_ref, hits):
    ## eval array of vectors connecting all hits with ref hit
    dr = hits - hit_ref

    # dr_u =  np.abs(dr[:,1])
    # dr_d =  np.abs(dr[:,0])
    ## projections
    dr_u =  dr[:,1]
    dr_d =  dr[:,0]
    ## eval module
    dr_mod =  np.sqrt(dr_u*dr_u + dr_d*dr_d)
    return dr_mod