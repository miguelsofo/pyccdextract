import extract
    
config = {
    "input":"/media/kiwi/ExtremeSSD/datos-ACDS/imgs_acds_xrays_skp_loop2/proc/ovs_report.csv",
    "out_file":"/media/kiwi/ExtremeSSD/datos-ACDS/imgs_acds_xrays_skp_loop2/proc/events.hdf5",
}

extract.runExtract(config)
