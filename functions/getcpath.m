function folder = getcpath()

folder = which('MOBYDICCFolder');
folder = folder(1:end-16);