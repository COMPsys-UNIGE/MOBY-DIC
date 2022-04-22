function folder = getsimulinkpath()

folder = which('MOBYDICSimulinkFolder');
folder = folder(1:end-23);