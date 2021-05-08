import extract_nico

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
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW30/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW30/proc/events.hdf5",
    #"noise":[170,170,170,170] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW40/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW40/proc/events.hdf5",
    #"noise":[220,220,220,220] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW50/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW50/proc/events.hdf5",
    #"noise":[250,250,250,250] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW60/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW60/proc/events.hdf5",
    #"noise":[290,290,290,290] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW100/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW100/proc/events.hdf5",
    #"noise":[430,430,430,430] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW120/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW120/proc/events.hdf5",
    #"noise":[475,475,475,475] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW200/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_buff_xrays/IW200/proc/events.hdf5",
    #"noise":[720,720,720,720] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_acds_xrays/iw200/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_acds_xrays/iw200/proc/events.hdf5",
    #"noise":[160,160,160,160] # HDU noise in ADU
    
    #"in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_acds_xrays_skp/iw30/proc",
    #"out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_acds_xrays_skp/iw30/proc/events.hdf5",
    #"noise":[320,320,320,320] # HDU noise in ADU
    
    "in_folder":"/Volumes/ExtremeSSD/datos-ACDS/imgs_rnd/proc",
    "out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_rnd/proc/events.hdf5",
    "noise":[2050,2050,2050,2050] # HDU noise in ADU
}
extract_nico.runExtract(config)
