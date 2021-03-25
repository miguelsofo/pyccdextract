import extract

config = {
    #"in_folder":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/DM2-16Mar2020_toti-BRC_2",
    #"out_file":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/DM2-16Mar2020_toti-BRC_2/events.hdf5",
    #"gain":[596.56,536.86,487.86,611.41],
    #"noise":[114.25,125.26,120.35,123.23]
    
    #"in_folder":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/22Ago2020_toti_run1",
    #"out_file":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/22Ago2020_toti_run1/events.hdf5",
    #"gain":[596.56,536.86,487.86,611.41],
    #"noise":[114.25,125.26,120.35,123.23]
    
    #"in_folder":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/phys/2020-09-05-bin1x5-rolling",
    #"out_file":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/phys/2020-09-05-bin1x5-rolling/events.hdf5",
    #"gain":[1.0,1.0,1.0,1.0],
    #"noise":[0.1,0.1,0.1,0.1]
    
    #"in_folder":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/2020-09-05-bin1x5-rolling",
    #"out_file":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/2020-09-05-bin1x5-rolling/events.hdf5",
    #"gain":[603.03,553.57,505.71,620.70], #from XB
    #"noise":[108.9054, 121.4884, 121.3152, 124.218 ]
    
    #"in_folder":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/testimgs",
    #"out_file":"/Users/kiwi/Documents/PROYECTOS/DM2/datos/testimgs/events.hdf5",
    #"gain":[603.03,553.57,505.71,620.70], #from XB
    #"noise":[108.9054, 121.4884, 121.3152, 124.218 ]

    # ===================================================================================
    # Regular-CCD (or 1-sample images)    
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/xray-IW-scan/IW100/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/xray-IW-scan/IW100/proc/events_IW100.hdf5",
    #"noise":[14.0,15.0,15.0,15.0] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/xray-IW-scan/IW200/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/xray-IW-scan/IW200/proc/events_IW200.hdf5",
    #"noise":[13.0,15.0,16.0,18.0] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/xray-IW-scan/IW300/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/xray-IW-scan/IW300/proc/events_IW300.hdf5",
    #"noise":[12.0,20.0,20.0,20.0] # HDU noise in ADU
}
extract.runExtract(config)
