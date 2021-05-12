import extract

configA = {
    # Example:
    #   - Regular CCD
    #   - Input folder
    "input":"/Volumes/ExtremeSSD/datos-ACDS/test/proc",
    "out_file":"/Volumes/ExtremeSSD/datos-ACDS/test/event.hdf5",
    "noise":[24,24,24,24] # HDU noise in ADU
}
    
configB = {
    # Example:
    #   - Regular CCD
    #   - Input csv
    "input":"/Volumes/ExtremeSSD/datos-ACDS/imgs_acds_skp_loop/proc/ovs_report.csv",
    "out_file":"/Volumes/ExtremeSSD/datos-ACDS/imgs_acds_skp_loop/proc/events.hdf5",
}

configTest = {
    # Example:
    #   - Regular CCD
    #   - Input folder
    "input":"/Volumes/ExtremeSSD/datos-ACDS/test",
    "out_file":"/Volumes/ExtremeSSD/datos-ACDS/test/event.hdf5",
    "T1":[234,567,123,100] # HDU noise in ADU
}

extract.runExtract(configTest)
