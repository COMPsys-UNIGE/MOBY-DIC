function folder = getvhdlpath()

folder = which('MOBYDICVHDLFolder');
folder = folder(1:end-19);