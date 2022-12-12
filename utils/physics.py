def dr(hit_ref, hits):
    ## eval array of vectors connecting all hits with ref hit
    dr = hits - hit_ref

    # dr_u =  np.abs(dr[:,1])
    # dr_d =  np.abs(dr[:,0])
    ## projections
    dr_u =  dr[:,1]
    dr_d =  dr[:,0]
    ## eval module
    dr_mod =  (dr_u*dr_u + dr_d*dr_d)**0.5
    return dr_mod

def overlap_length(int1,int2):
    """
    arguments : 
        - intervals encoded as list of 2 boundaries
    returns   : 
        - length of the overlap between the 2 intervals
    """
    if max(int1)<=min(int2) or max(int2)<=min(int1):
        return 0
    else:
        return min(max(int1)-min(int2),max(int2)-min(int1))

# Returns the distance between 
def dist_line_pt(t,x0,point):
    """
    arguments : 
        - t     :   tangent of the line
        - x0    :   intercept of the line at z = 0
        - point :   either a list of 2 coordinates or a numpy.array of 2 components (x,z)
    returns   : 
        - d     :   distance between the line and the point
    """
    d = abs(point[1]*t-point[0]+x0)/(t**2+1)**0.5
    return d

def dist_line_rect(t,x0,center,height,width):
    """
    arguments : 
        - t      :  tangent of the line
        - x0     :  intercept of the line at z = 0
        - center :  either a list of 2 coordinates or a numpy.array of 2 components (x,z) corresponding to the center of the rectangle
        - height :  height of the rectangle (z axis)
        - width  :  width of the rectangle (x axis)
    returns   : 
        - d     :   distance between the line and the rectancle
    """
    interval_x_rect = [center[0]-width/2,center[0]+width/2]
    interval_x_line = [t*(center[1]-height/2)+x0,t*(center[1]+height/2)+x0]
    if overlap_length(interval_x_rect,interval_x_line) > 0:
        return 0
    elif t == 0:
        return min(abs(center[0]-width/2-x0),abs(center[0]+width/2-x0))
    elif t > 0:
        return min(dist_line_pt(t,x0,[center[0]+width/2,center[1]-height/2]),dist_line_pt(t,x0,[center[0]-width/2,center[1]+height/2]))
    else:
        return min(dist_line_pt(t,x0,[center[0]+width/2,center[1]+height/2]),dist_line_pt(t,x0,[center[0]-width/2,center[1]-height/2]))

