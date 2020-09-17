import numpy as np
from uncertainties import ufloat
# Foley et al LOSS sample likelihoods
def morphology_likelihood(hubbletype):
    # hubbletype is list of the useful classify fractions from table to determine a morphology likelihood
    [funclassifiable,fps,fsphere,fdisk,firr,fds,fspiral,fbar,ftidalarms,fDiskDom,fBulgeDom] = hubbletype
    if all(i <= fps for i in hubbletype) or all(i <= funclassifiable for i in hubbletype):
        # point source or unclassifiable ~ no sn type likelihood available 
        return  -99.0
    else:
        # elliptical
        lIaE = 0.141
        dlIaE = [0.021,0.018] # pm
        lCCE = 0.002
        dCCE = [0.003,0.001]
        lIaE_upper = ufloat(lIaE,dlIaE[0])
        lIaE_lower = ufloat(lIaE,dlIaE[1])
        lCCE_upper = ufloat(lCCE,dCCE[0])
        lCCE_lower = ufloat(lCCE,dCCE[1])
        # irregular
        lIaIrr = 0.003
        dlIaIrr = [0.004,0.002] # pm
        lCCIrr = 0.015
        dCCIrr = [0.006,0.004]
        lIaIrr_upper = ufloat(lIaIrr,dlIaIrr[0])
        lIaIrr_lower = ufloat(lIaIrr,dlIaIrr[1])
        lCCIrr_upper = ufloat(lCCIrr,dCCIrr[0])
        lCCIrr_lower = ufloat(lCCIrr,dCCIrr[1])
        # disk+sph ~ lenticular
        lIaS0 = 0.217
        dlIaS0 = [0.0026,0.0023] # pm
        lCCS0 = 0.017
        dCCS0 = [0.007,0.005]
        lIaS0_upper = ufloat(lIaS0,dlIaS0[0])
        lIaS0_lower = ufloat(lIaS0,dlIaS0[1])
        lCCS0_upper = ufloat(lCCS0,dCCS0[0])
        lCCS0_lower = ufloat(lCCS0,dCCS0[1])
        # Here begin the spirals, a->d tight arms big nucleus->looser more arms
        # I make an average spiral type likelihood
        # Sa
        lIaSa = 0.149
        dlIaSa = [.022,0.019] # pm
        lCCSa = .142
        dCCSa = [0.017,0.015]
        # Sb
        lIaSb = 0.177
        dlIaSb = [0.023,0.021] # pm
        lCCSb = 0.188
        dCCSb = [0.020,0.018]
        # Sbc
        lIaSbc = 0.117
        dlIaSbc = [0.019,0.017] # pm
        lCCSbc = 0.231
        dCCSbc = [0.022,0.020]
        # Sc
        lIaSc = 0.120
        dlIaSc = [0.019,0.017] # pm
        lCCSc = 0.218
        dCCSc = [0.021,0.019]
        # Scd
        lIaScd = 0.076
        dlIaScd = [0.016,0.013] # pm
        lCCScd = 0.188
        dCCScd = [0.020,0.018]

        """
        # this doesn't make sense for the uncertainty on spiral, the np average of upper and lower returns single uncertainty mixing them up
        # I will do w ufloat again
        # grouping the spirals together to make aforementioned spiral average likelihood
        # may return to to see if simple to recalc likelihoods w this grouping
        
        lIaS = np.average([lIaSa,lIaSb,lIaSbc,lIaSc,lIaScd])
        dlIaS = np.average([dlIaSa,dlIaSb,dlIaSbc,dlIaSc,dlIaScd])
        lCCS = np.average([lCCSa,lCCSb,lCCSbc,lCCSc,lCCScd])
        dCCS = np.average([dCCSa,dCCSb,dCCSbc,dCCSc,dCCScd])

        # use all the classification tools to get weighted likelihood
        # E,S0,Irr,S~(Sa,Sb,Sbc,Sc,Scd)...fsph,fdisk+sph,firr,fspiral
        k = fsphere + fds + firr + fspiral # weights normalization
        lIa = (fsphere*lIaE+fds*lIaS0+firr*lIaIrr+fspiral*lIaS)/k 
        lCC = (fsphere*lCCE+fds*lCCS0+firr*lCCIrr+fspiral*lCCS)/k
        dlIa = (fsphere*dlIaE+fds*dlIaS0+firr*dlIaIrr+fspiral*dlIaS)/k
        dCC = (fsphere*dCCE+fds*dCCS0+firr*dCCIrr+fspiral*dCCS)/k
        """

        # calcs w error analysis built in for combining the spiral likelihoods...& then combining the different morphologies
        # upper + ~ 0 idx;lower - ~ 1 idx... dlIas, dlCCs are lists [+delta,-delta]
        lIaS = [lIaSa,lIaSb,lIaSbc,lIaSc,lIaScd]
        dlIaS = [dlIaSa,dlIaSb,dlIaSbc,dlIaSc,dlIaScd]
        lCCS = [lCCSa,lCCSb,lCCSbc,lCCSc,lCCScd]
        dlCCS = [dCCSa,dCCSb,dCCSbc,dCCSc,dCCScd]

        upper_lIaS_ufloats = []
        upper_lCCS_ufloats = []
        lower_lIaS_ufloats = []
        lower_lCCS_ufloats = []
        for i in range(len(lIaS)):
            upper_lIaS_ufloats.append(ufloat(lIaS[i],dlIaS[i][0]))
            upper_lCCS_ufloats.append(ufloat(lCCS[i],dlCCS[i][0]))
            lower_lIaS_ufloats.append(ufloat(lIaS[i],dlIaS[i][1]))
            lower_lCCS_ufloats.append(ufloat(lCCS[i],dlCCS[i][1]))
        # take the average of the ufloats for upper and lower spiral likelihoods
        # for some unknown reason np average no good here, do it ye olde fashioned way
        lIaS_upper = sum(upper_lIaS_ufloats)/len(upper_lIaS_ufloats)
        lCCS_upper = sum(upper_lCCS_ufloats)/len(upper_lCCS_ufloats)
        lIaS_lower = sum(lower_lIaS_ufloats)/len(lower_lIaS_ufloats) 
        lCCS_lower = sum(lower_lCCS_ufloats)/len(lower_lCCS_ufloats)
        # get the nominal values and std devs for spirals averaged  
        lIaS = lIaS_upper.nominal_value
        dlIaS = [lIaS_upper.std_dev,lIaS_lower.std_dev]
        lCCS = lCCS_upper.nominal_value
        dCCS = [lCCS_upper.std_dev,lCCS_lower.std_dev]

        # now use all the classification tools (different morphology fractions) to get weighted likelihood
        # E,S0,Irr,S~(Sa,Sb,Sbc,Sc,Scd)...fsph,fdisk+sph,firr,fspiral
        try:
            # [0] on fracs gives value so they don't come as column objects...ie with their names, if thats what you desire .name
            k = fsphere[0] + fds[0] + firr[0] + fspiral[0] # weights normalization
            lIa_upper = (fsphere[0]*lIaE_upper+fds[0]*lIaS0_upper+firr[0]*lIaIrr_upper+fspiral[0]*lIaS_upper)/k 
            lCC_upper = (fsphere[0]*lCCE_upper+fds[0]*lCCS0_upper+firr[0]*lCCIrr_upper+fspiral[0]*lCCS_upper)/k
            lIa_lower =  (fsphere[0]*lIaE_lower+fds[0]*lIaS0_lower+firr[0]*lIaIrr_lower+fspiral[0]*lIaS_lower)/k 
            lCC_lower =  (fsphere[0]*lIaE_lower+fds[0]*lIaS0_lower+firr[0]*lIaIrr_lower+fspiral[0]*lIaS_lower)/k 
        except:
            # here in case provided fractions by hand (ie not column object no need to index)
            k = fsphere + fds + firr + fspiral # weights normalization
            lIa_upper = (fsphere*lIaE_upper+fds*lIaS0_upper+firr*lIaIrr_upper+fspiral*lIaS_upper)/k 
            lCC_upper = (fsphere*lCCE_upper+fds*lCCS0_upper+firr*lCCIrr_upper+fspiral*lCCS_upper)/k
            lIa_lower =  (fsphere*lIaE_lower+fds*lIaS0_lower+firr*lIaIrr_lower+fspiral*lIaS_lower)/k 
            lCC_lower =  (fsphere*lIaE_lower+fds*lIaS0_lower+firr*lIaIrr_lower+fspiral*lIaS_lower)/k 

        lIa = lIa_upper.nominal_value
        lCC = lCCS_upper.nominal_value
        dlIa = [lIa_upper.std_dev,lIa_lower.std_dev]
        dCC = [lCC_upper.std_dev,lCC_lower.std_dev]
        
        #dlIa = (fsphere*dlIaE+fds*dlIaS0+firr*dlIaIrr+fspiral*dlIaS)/k
        #dCC = (fsphere*dCCE+fds*dCCS0+firr*dCCIrr+fspiral*dCCS)/k


        # other useful info to gather
        elliptical = (fsphere,lIaE,lCCE,dlIaE,dCCE)
        lenticular = (fds,lIaS0,lCCS0,dlIaS0,dCCS0)
        irregular = (firr,lIaIrr,lCCIrr,dlIaIrr,dCCIrr)
        spiral = (fspiral,lIaS,lCCS,dlIaS,dCCS)
        # all together now!
        # this is necessary for plot params classify if want to use shapes to show morphology 
        hubblefractions = [elliptical,lenticular,irregular,spiral]
        return [lIa,lCC,dlIa,dCC,hubblefractions]
    
def color_likelihood(color):
    if color < 0 or color > 6.25:
        # not in one of the bins avail
        return -99.0
    else: # find the bin
        if color >= 0 and color < 1.75:
            lIa = 0.026
            dlIa = [0.010,0.007] # pm
            lCC = 0.044
            dCC = [0.010,0.008]
        elif color >= 1.75 and color < 2.25:
            lIa = .023
            dlIa = [.010,.007]
            lCC = .075
            dCC = [.013,.011]
        elif color >= 2.25 and color < 2.5:
            lIa = .037
            dlIa = [.012,.009]
            lCC = .069
            dCC = [.013,.011]
        elif color >= 2.5 and color < 2.75:
            lIa = .043
            dlIa = [.013,.010]
            lCC = .115
            dCC = [.016,.014]
        elif color >= 2.75 and color < 3.0:
            lIa = .111
            dlIa = [.019,.016]
            lCC = .181
            dCC = [.020,.018]
        elif color >= 3.0 and color < 3.25:
            lIa = .131
            dlIa = [.021,.018]
            lCC = .185
            dCC = [.020,.018]
        elif color >= 3.25 and color < 3.5:
            lIa = .154
            dlIa = [.022,.020]
            lCC = .137
            dCC = [.018,.016]
        elif color >= 3.5 and color < 3.75:
            lIa = .125
            dlIa = [.020,.018]
            lCC = .115
            dCC = [.016,.014]
        elif color >= 3.75 and color < 4.0:
            lIa = .168
            dlIa = [.023,.021]
            lCC = .048
            dCC = [.011,.009]
        elif color >= 4.0 and color < 4.25:
            lIa = .103
            dlIa = [.019,.016]
            lCC = .022
            dCC = [.008,.006]
        elif color >= 4.25 and color < 6.25:
            lIa = .08
            dlIa = [.017,.014]
            lCC = .01
            dCC = [.006,.004]
        return [lIa,lCC,dlIa,dCC]
        
def MK_likelihood(MK):
    if MK > -17.1 or MK < -26.5:
        # value not in LOSS bin
        return -99.0
    else: # find bin
        if MK >= -21.5 and MK < -17.1:
            lIa = 0.020
            dlIa = [0.009,0.006] # pm
            lCC = .055
            dCC = [0.011,0.009]
        elif MK >= -22.3 and MK < -21.5:
            lIa = .037
            dlIa = [.012,.009]
            lCC = .060
            dCC = [.012,.010]
        elif MK >= -22.7 and MK < -22.3:
            lIa = .031
            dlIa = [.011,.008]
            lCC = .058
            dCC = [.012,.010]
        elif MK >= -23.1 and MK < -22.7:
            lIa = .048
            dlIa = [.013,.010]
            lCC = .115
            dCC = [.016,.014]
        elif MK >= -23.5 and MK < -23.1:
            lIa = .040
            dlIa = [.012,.009]
            lCC = .133
            dCC = [.017,.015]
        elif MK >= -23.9 and MK < -23.5:
            lIa = .130
            dlIa = [.021,.018]
            lCC = .136
            dCC = [.017,.015]
        elif MK >= -24.3 and MK < -23.9:
            lIa = .153
            dlIa = [.022,.019]
            lCC = .125
            dCC = [.017,.015]
        elif MK >= -24.7 and MK < -24.3:
            lIa = .190
            dlIa = [.025,.022]
            lCC = .146
            dCC = [.018,.016]
        elif MK >= -25.1 and MK < -24.7:
            lIa = .164
            dlIa = [.023,.020]
            lCC = .125
            dCC = [.017,.015]
        elif MK >= -25.5 and MK < -25.1:
            lIa = .108
            dlIa = [.019,.016]
            lCC = .029
            dCC = [.009,.007]
        elif MK >= -26.5 and MK < -25.5:
            lIa = .079
            dlIa = [.016,.014]
            lCC = .018
            dCC = [.007,.005]
        return [lIa,lCC,dlIa,dCC]

def eo_likelihood(eo):
    if eo < 0 or eo > 5.25:
        return -99.0
    else:
        if eo >= 0 and eo < .05:
            lIa = 0.043
            dlIa = [0.012,0.010] # pm
            lCC = 0.032
            dCC = [0.009,0.007]
        elif eo >= .05 and eo < 0.1:
            lIa = .098
            dlIa = [.018,.015]
            lCC = .086
            dCC = [.014,.012]
        elif eo >= 0.1 and eo < 0.15:
            lIa = .130
            dlIa = [.020,.018]
            lCC = .091
            dCC = [.014,.012]
        elif eo >= .15 and eo < 0.2:
            lIa = .09
            dlIa = [.017,.014]
            lCC = .091
            dCC = [.014,.012]
        elif eo >= 0.2 and eo < 0.25:
            lIa = .071
            dlIa = [.015,.013]
            lCC = .114
            dCC = [.016,.014]
        elif eo >= 0.25 and eo < 0.3:
            lIa = .057
            dlIa = [.014,.011]
            lCC = .084
            dCC = [.013,.012]
        elif eo >= 0.3 and eo < 0.35:
            lIa = .068
            dlIa = [.015,.012]
            lCC = .058
            dCC = [.011,.009]
        elif eo >= .35 and eo < 0.4:
            lIa = .060
            dlIa = [.014,.011]
            lCC = .054
            dCC = [.011,.009]
        elif eo >= 0.4 and eo < 0.45:
            lIa = .038
            dlIa = [.012,.009]
            lCC = .067
            dCC = [.012,.010]
        elif eo >= 0.45 and eo < 0.5:
            lIa = .052
            dlIa = [.013,.011]
            lCC = .063
            dCC = [.012,.01]
        elif eo >= 0.5 and eo < 0.6:
            lIa = .060
            dlIa = [.014,.011]
            lCC = .089
            dCC = [.014,.012]
        elif eo >= .6 and eo < 0.75:
            lIa = .062
            dlIa = [.014,.012]
            lCC = .048
            dCC = [.01,.009]
        elif eo >= .75 and eo < 1:
            lIa = .073
            dlIa = [.016,.013]
            lCC = .045
            dCC = [.01,.008]    
        elif eo >= 1 and eo < 1.4:
            lIa = .052
            dlIa = [.013,.011]
            lCC = .041
            dCC = [.01,.008]    
        elif eo >= 1.4 and eo < 5.25:
            lIa = .046
            dlIa = [.013,.010]
            lCC = .037
            dCC = [.009,.007]
        return [lIa,lCC,dlIa,dCC]
        