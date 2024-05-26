import animation_core.visualizer as ac
import numpy as np
import pandas as pd

df = pd.read_csv('_out/simulation_results.txt')

if 'type' in df.columns:
    if df['type'].iloc[0] == 'attitude':
        qlist = []
        rlist = []
        vlist = []
        for i in range(len(df)):
            qlist.append(np.array([ df['q1'].iloc[i], df['q2'].iloc[i], 
                                    df['q3'].iloc[i], df['q4'].iloc[i] ]))
            
            rlist.append(np.array([ df['x'].iloc[i], df['y'].iloc[i], df['z'].iloc[i] ]))
            vlist.append(np.array([ df['vx'].iloc[i], df['vy'].iloc[i], df['vz'].iloc[i] ]))

        ac.runVisualizer(qlist, rlist, vlist, 0.01, buff_amt=100)

    elif df['type'].iloc[1] == 'orbit':
        pass

    else:
        print('unexpected simulation type, idk what to do with this??')

else:
    print('unexpected simulation type, idk what to do with this??')