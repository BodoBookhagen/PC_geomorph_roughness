def evans_topo(z, spacing):

    def EVANS_FX(window):
        s = int(np.sqrt(window.shape[0]))
        window = window.reshape((s,s))
        #Hengl and Reuter (2009)
        Z1, Z2, Z3 = window[0,:]
        Z4, Z5, Z6 = window[1,:]
        Z7, Z8, Z9 = window[2,:]
        
        fx = (Z3 + Z6 + Z9 - Z1 - Z4 - Z7) / (6*spacing)
        #fx = ((Z3 + Z6 + Z9)/3 - (Z1 + Z4 + Z7)/3) / (2*spacing)
        #fy = (Z1 + Z2 + Z3 - Z7 - Z8 - Z9) / (6*spacing)
        
        #dx = (window[1,2] + window[0,2] + window[2,2])/3. - (window[1,0] + window[0,0] + window[2,0])/3.
        #a = (Z3 + Z6 + Z9) / 3
        #b = (Z2 + Z5 + Z8) / 3
        #c = (Z1 + Z4 + Z7) / 3
        #d = spacing
        
        #fx = (Z6 - Z4) / 2*d
        #fxx = (((a - b) / d) - ((b - c) / d)) / d
        
        #a = (Z1 + Z2 + Z3) / 3
        #b = (Z4 + Z5 + Z6) / 3
        #c = (Z7 + Z8 + Z9) / 3
        #d = spacing
        
        #fy = (Z2 - Z8) / 2*d
        #fyy = (((a - b) / d) - ((b - c) / d)) / d
        
        #fxy = (-Z1 + Z3 + Z7 - Z9) / (4*(spacing**2))
        
        ##From Schmidt et al. 2003
        #mean_curv = -1*(((1 + fy**2)*fxx-2*fxy*fx*fy+(1+fx**2)*fyy)/(2*(fx**2 + fy**2 + 1)**1.5))
        #gauss_curv = (fxx*fyy-fxy**2)/((fx**2 + fy**2 + 1)**2)
        #total_curv = fxx**2 + 2*fxy**2 + fyy**2
        #profile_curv = -1*((fxx*fx**2+2*fxy*fx*fy+fyy*fy**2)/((fx**2+fy**2)*(fx**2+fy**2+1)**1.5))
        #contour_curv = -1*((fxx*fx**2+2*fxy*fx*fy+fyy*fy**2)/((fx**2+fy**2+1)**1.5))
        
        #mag = np.sqrt(fx**2 + fy**2)
        #return np.arctan(mag) * 180 / np.pi
        return fx
    
    def EVANS_FY(window):
        s = int(np.sqrt(window.shape[0]))
        window = window.reshape((s,s))
        #FLORINSKY 1998 (p, q = dx, dy)
        Z1, Z2, Z3 = window[0,:]
        Z4, Z5, Z6 = window[1,:]
        Z7, Z8, Z9 = window[2,:] 
        
        #window[1,2] - window[1,0]
        
        #a = (Z3 + Z6 + Z9) / 3
        #b = (Z2 + Z5 + Z8) / 3
        #c = (Z1 + Z4 + Z7) / 3
        #d = spacing
        
        #fx = (Z6 - Z4) / 2*d
        #fxx = (((a - b) / d) - ((b - c) / d)) / d
        
        #fx = (Z3 + Z6 + Z9 - Z1 - Z4 - Z7) / (6*spacing)
        fy = (Z1 + Z2 + Z3 - Z7 - Z8 - Z9) / (6*spacing)
        #fy = ((Z1 + Z2 + Z3)/3 - (Z7 + Z8 + Z9)/3) / (2*spacing)
        #dy = window[2,1] - window[0,1]
        #a = (Z1 + Z2 + Z3) / 3
        #b = (Z4 + Z5 + Z6) / 3
        #c = (Z7 + Z8 + Z9) / 3
        #d = spacing
        
        #fy = (Z2 - Z8) / 2*d
        #fyy = (((a - b) / d) - ((b - c) / d)) / d
        
        #fxy = (-Z1 + Z3 + Z7 - Z9) / (4*(spacing**2))
        
        ##From Schmidt et al. 2003
        #mean_curv = -1*(((1 + fy**2)*fxx-2*fxy*fx*fy+(1+fx**2)*fyy)/(2*(fx**2 + fy**2 + 1)**1.5))
        #gauss_curv = (fxx*fyy-fxy**2)/((fx**2 + fy**2 + 1)**2)
        #total_curv = fxx**2 + 2*fxy**2 + fyy**2
        #profile_curv = -1*((fxx*fx**2+2*fxy*fx*fy+fyy*fy**2)/((fx**2+fy**2)*(fx**2+fy**2+1)**1.5))
        #contour_curv = -1*((fxx*fx**2+2*fxy*fx*fy+fyy*fy**2)/((fx**2+fy**2+1)**1.5))
        
        #return (np.arctan2(-1*fy, fx) * 180 / np.pi - 180)*-1
        #return (np.arctan2(fy, fx) * 180 / np.pi - 180)*-1
        return fy
    
    dx = np.empty(z.shape)
    generic_filter(z, EVANS_FX, size=(3,3), footprint=None, output=dx, mode='constant', cval=np.nan, origin=0)
    
    dy = np.empty(z.shape)
    generic_filter(z, EVANS_FY, size=(3,3), footprint=None, output=dy, mode='constant', cval=np.nan, origin=0)

    asp = (np.arctan2(-1*dy, dx) * 180 / np.pi - 180)*-1
    mag = np.sqrt(dx**2 + dy**2)
    slp = np.arctan(mag) * 180 / np.pi
    
    return slp, asp