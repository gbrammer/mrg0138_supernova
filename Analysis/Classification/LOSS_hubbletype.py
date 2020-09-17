def hubbletype(morphology,weighted=True):
    # morphology is the visual classification table 
    # many independent classifications are made on the galaxies
    # they are allowed to select various parameters
    # the table reports fraction of classifiers who select a given parameter
    
    # this function takes the useful classification fractions for selecting Hubble Type 
    if morphology == -99.0: # don't have morphology ~ visual classification table
        return -99.0
    else:
        # base classifications  
        funclassifiable = morphology['f(Unclassifiable)']
        fps = morphology['f(PointSource)']
        fsphere = morphology['f(Spheroid)']
        fdisk = morphology['f(Disk)']
        firr = morphology['f(Irregular)']
        # other potentially useful?
        fds = morphology['f(DS)'] # disk and sphere
        fspiral = morphology['f(Spiral)']
        fbar = morphology['f(Bar)']
        ftidalarms = morphology['f(TidalArms)'] # more of an interaction/merger feature
        fDiskDom = morphology['f(DiskDom)']
        fBulgeDom = morphology['f(BulgeDom)']

        classification_fractions = [funclassifiable,fps,fsphere,fdisk,firr,fds,fspiral,fbar,ftidalarms,fDiskDom,fBulgeDom]
        return classification_fractions

    