function folder = getcpppath()

folder = which('MOBYDICC++Folder');
folder = folder(1:end-18);