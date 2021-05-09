import extract

config = {
    # Example:
    #   - Regular CCD
    #   - Input folder
    "input":"/Volumes/ExtremeSSD/datos-ACDS/test/proc",
    "out_file":"/Volumes/ExtremeSSD/datos-ACDS/test/event.hdf5",
    "noise":[24,24,24,24] # HDU noise in ADU
}
extract.runExtract(config)
