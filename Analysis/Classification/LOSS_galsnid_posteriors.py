from LOSS_likelihoods import color_likelihood,MK_likelihood,eo_likelihood,morphology_likelihood
from LOSS_hubbletype import hubbletype
import numpy as np
from uncertainties import ufloat
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70, Om0=0.3, Ob0=0.05)

def galsnid_params(eo,z,RFcolors,Morphology):
    # parameters used from SN host to assign probability of Ia
    # effective offset, K absolute mag, B-K color, morphology

    # eo used to determine likely host, float val already available
    # RFcolors should be the 3dhst table (the row for this obj)
    # Morphology should be the visual classification table (obj row)
    
    # rest frame K & B band flux densities
    fnuK = RFcolors['L163'][0] 
    fnuB = RFcolors['L136'][0] 
    m_K = -2.5*np.log10(fnuK) + 25
    m_B = -2.5*np.log10(fnuB) + 25
    # flcdm distance modulus ~ m - M
    mu = cosmo.distmod(z)
    # absolute mags
    MK = m_K - mu.value[0]
    MB = m_B - mu.value[0]
    Color = MB-MK # B-K

    # this hubbletype is the classification fractions 
    HubbleType = hubbletype(Morphology)
    return [eo,MK,Color,HubbleType]

def plot_params_colorbarclassification(SNlist,x_param='MK',y_param='Color',pIalower=0.4,pIaupper=0.9,prior=0.5):
    #want to show the classification probabilities with the parameters 
    #ie a scatter marker location will show the parameter values and a color the classification Ia prob 
    
    # the SNlist should have [SN,LOSS_params,LOSS_galsnids]
    # ie need to loop through galsnid for each candels SN to make & pickle list
    # then load it in to feed here for the plot
    # at some point may want to include a loop w galsnid inside this function
    # that way could save a step if want to play w prior 
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.cm import coolwarm

    fig,ax = plt.subplots()
    cax = fig.add_axes([0.27,1,0.5,0.05]) # initialize the colorbar axis, location

    # LOSS_params will have list like [eo,MK,Color,HubbleType]
    if x_param=='eo':
        x_idx =  0
    elif x_param=='MK':
        x_idx = 1
    elif x_param=='Color':
        x_idx = 2
    elif x_param=='HubbleType':
        x_idx = 3
    else:
        print('Your x param wasnt given as eo,MK,Color,HubbleType')
        return 
    if y_param=='eo':
        y_idx =  0
    elif y_param=='MK':
        y_idx = 1
    elif y_param=='Color':
        y_idx = 2
    elif y_param=='HubbleType':
        y_idx = 3
    else:
        print('Your y param wasnt given as eo,MK,Color,HubbleType')
        return 

    # list of host values in parameter space and the classification prob for colorbar
    x_vals = []
    y_vals = []
    vs = []
    for i in SNlist:
        SN = i[0]
        LOSS_params = i[1]
        LOSS_galsnid = i[2]
        
        # the host in parameter space 
        x_vals.append(LOSS_params[x_idx])
        y_vals.append(LOSS_params[y_idx])

        # classification for the color bar setting
        vIa = LOSS_galsnid[0][0] # pIa ufloat; is tuple for an upper and lower error
        # [1] would be the CC ufloat also as tuple for upper and lower err
        # access prob use attribute .nominal_value, or error using .std_dev
        vs.append(vIa.nominal_value)

    s = ax.scatter(x_vals,y_vals,c=vs,cmap=coolwarm)
    fig.colorbar(s,cax=cax,orientation='horizontal')
    print('done did it')
    return 

def galsnid(params,prior=0.5):
    # params should be a list [eo,MK,color,hubbletype]
    eo,MK,color,HubbleType=params
    # can't id if don't have any of the parameters
    if eo == -99.0 and MK == -99.0 and color == -99.0 and HubbleType == -99.0:
        return -99.0
    else:
        avail_params = []
        for i in range(len(params)):
            if params[i] == -99.0:
                # we don't have this parameter can't use for id
                pass
            else:
                avail_params.append([i,params[i]])
        
        used_params = [] # if its avail it still needs to be in loss bin to be used
        used_paramflags = [] # 0,1,2,3 ~ eo,MK,col,hubtype...if want to know which used 
        for i in range(len(avail_params)):
            if avail_params[i][0] == 0:
                # eo
                if eo_likelihood(eo) != -99.0: # falls into LOSS bin
                    used_params.append(eo_likelihood(eo))
                    used_paramflags.append(avail_params[i][0])
            elif avail_params[i][0] == 1:
                # MK
                if MK_likelihood(MK) != -99.0: #in LOSS bin
                    used_params.append(MK_likelihood(MK)) 
                    used_paramflags.append(avail_params[i][0])
            elif avail_params[i][0] == 2:
                # color
                if color_likelihood(color) != -99.0: #in LOSS bin
                    used_params.append(color_likelihood(color))
                    used_paramflags.append(avail_params[i][0])
            elif avail_params[i][0] == 3:
                # hubble type
                if morphology_likelihood(HubbleType) != -99.0: # in LOSS bin
                    used_params.append(morphology_likelihood(HubbleType))
                    used_paramflags.append(avail_params[i][0])
        lIas,lCCs,dlIas,dlCCs = [],[],[],[]
        for i in range(len(used_params)):
            lIas.append(used_params[i][0])
            lCCs.append(used_params[i][1])
            dlIas.append(used_params[i][2])
            dlCCs.append(used_params[i][3])

        """
        no longer need this block
        I did the calc this way before including error analysis
        the nominal values are consistent w results from ufloat
        # the normalization
        k = prior*np.prod(lIas) + (1-prior)*np.prod(lCCs)
        # the posterior, p(Ia|D)
        pIa = prior*np.prod(lIas)/k
        # p(CC|D), calculate this as sanity check for normalization... pIa + pCC = 1
        pCC = (1-prior)*np.prod(lCCs)/k
        """

        # calcs w error analysis built in
        # upper + ~ 0 idx;lower - ~ 1 idx... dlIas, dlCCs are lists [+delta,-delta]
        upper_lIa_ufloats = []
        upper_lCC_ufloats = []
        for i in range(len(used_params)):
            upper_lIa_ufloats.append(ufloat(lIas[i],dlIas[i][0]))
            upper_lCC_ufloats.append(ufloat(lCCs[i],dlCCs[i][0]))
        upper_k_ufloat = prior*np.prod(upper_lIa_ufloats) + (1-prior)*np.prod(upper_lCC_ufloats)

        lower_lIa_ufloats = []
        lower_lCC_ufloats = []
        for i in range(len(used_params)):
            lower_lIa_ufloats.append(ufloat(lIas[i],dlIas[i][1]))
            lower_lCC_ufloats.append(ufloat(lCCs[i],dlCCs[i][1]))
        lower_k_ufloat = prior*np.prod(lower_lIa_ufloats) + (1-prior)*np.prod(lower_lCC_ufloats)
        
        """
        some useful ufloat properties, ufloat.nominal_value,std_dev,error_components
        """

        pIa_upper_ufloat = prior*np.prod(upper_lIa_ufloats)/upper_k_ufloat
        pIa_lower_ufloat = prior*np.prod(lower_lIa_ufloats)/lower_k_ufloat
        pCC_upper_ufloat = (1-prior)*np.prod(upper_lCC_ufloats)/upper_k_ufloat
        pCC_lower_ufloat = (1-prior)*np.prod(lower_lCC_ufloats)/lower_k_ufloat

        return [(pIa_upper_ufloat,pIa_lower_ufloat),(pCC_upper_ufloat,pCC_lower_ufloat),used_paramflags]
