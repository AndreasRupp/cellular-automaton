mex('mex/CAM/cellular_automaton.cxx', '-Iinclude', 'COMPFLAGS=$COMPFLAGS -O3')

cellular_automaton(5,0.94,4,1,zeros(100,6))