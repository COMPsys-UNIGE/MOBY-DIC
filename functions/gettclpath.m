function folder = gettclpath()

folder = which('MOBYDICTCLFolder');
folder = folder(1:end-18);