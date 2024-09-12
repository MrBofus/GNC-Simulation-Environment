import numpy as np
import math
import matplotlib.pyplot as plt

import numba as nb
import warnings
warnings.filterwarnings('ignore')


## assumes no error :)) 
# (see top comment)
# https://stackoverflow.com/questions/37908104/lyapunov-exponent-python-implementation

# also i had to skip repeats which occur every int(physicsHz/softwareHz) iterations
# ( which equals counter below i did the math :) ) 

@nb.jit(nopython=True)
def _return_counter(dataset):
    first_counter = 0
    prev_step = dataset[0]
    for i in range(len(dataset)):
        if prev_step == dataset[i]:
            first_counter += 1
        else:
            break
    
    counter = 0
    prev_step = dataset[first_counter]
    for i in range(first_counter, len(dataset)):
        if prev_step == dataset[i]:
            counter += 1
        else:
            break
    
    return counter


@nb.jit(nopython=True)
def _drop_zeros(dataset):
    
    counter = 0
    for i in range(len(dataset)):
        if not dataset[i] == 0.0:
            break
        counter += 1

    return dataset[counter:]


@nb.jit(nopython=True)
def _remove_repeats(dataset):
    
    counter = _return_counter(dataset)

    clean_dataset = []
    for i in range(0, len(dataset), counter):
        clean_dataset.append(dataset[i])

    return clean_dataset


# @nb.jit(nopython=True)
def _compute_lyapunovs(dataset, eps, skip_step):

    clean_set = _remove_repeats(_drop_zeros(dataset))

    ll = []
    ts = []
    tcounter = 0

    N = len(clean_set)
 
    for i in range(0, N, skip_step):
        
        if i % 1000 == 0: print("\r" + str(int(100*i/N)), end="")
        lyapunovs = []
        for j in range(i + 1, N, skip_step):

            if N - j  < skip_step: continue

            # print(j)
            if np.abs(clean_set[i] - clean_set[j]) < eps and not np.abs(clean_set[i] - clean_set[j]) == 0:
                
                for k in range(1, min(N - i, N - j), skip_step):
                    d0 = np.abs(clean_set[i] - clean_set[j])
                    dn = np.abs(clean_set[i + k] - clean_set[j + k])
                    # print(d0, dn)
                    lyapunovs.append(math.log(dn) - math.log(d0))
                
        # print(lyapunovs)
        if len (lyapunovs) > 0:
            ll.append(np.mean(np.array(lyapunovs)))
            ts.append(tcounter)
            tcounter += 1

    return ts, ll

# ts, ll = _compute_lyapunovs(q1e, 0.01, 10)
# ts, ll = _compute_lyapunovs(q1e, 0.1, 200)


def _return_slices(dataset_list, n_dimensions, spread):
    
    qlist = []
    for dataset_ in dataset_list:

        dataset = _drop_zeros(dataset_)
        counter = _return_counter(dataset)

        q_s = []
        for i in range(n_dimensions):
            q_s.append(dataset[i*spread*counter:(-(n_dimensions-i)*spread*counter)])
        
        qlist.append(q_s)

    return qlist

def plot_parameter_space(dataset, n_dimensions, spread):

    qlist = _return_slices(dataset, n_dimensions, spread)

    if n_dimensions == 2:
        plt.figure(1)
        for q_s in qlist:
            plt.plot(q_s[0], q_s[1])
            plt.plot([0]*len(q_s[0]), q_s[1], 'k')
            plt.plot(q_s[0], [0]*len(q_s[1]), 'k')
        plt.grid()
        plt.show()

    elif n_dimensions == 3:
        figure = plt.figure()
        axes = figure.add_subplot(projection='3d')
        
        for q_s in qlist:
            axes.plot(q_s[0], q_s[1], q_s[2])
            axes.plot(q_s[0], [0]*len(q_s[1]), [0]*len(q_s[2]), 'k')
            axes.plot([0]*len(q_s[0]), q_s[1], [0]*len(q_s[2]), 'k')
            axes.plot([0]*len(q_s[0]), [0]*len(q_s[1]), q_s[2], 'k')
            

        plt.show()


if __name__ == "__main__":

    import pandas as pd

    df_spec = pd.read_csv('_out/simulation_results_spec.txt')
    df_jagsat = pd.read_csv('_out/simulation_results.txt')

    # [14:2420]

    # https://www.jovanodavic.com/publication/application-of-largest-lyapunov-exponent-analysis-on-the-studies-of-dynamics-under-external-forces/application-of-largest-lyapunov-exponent-analysis-on-the-studies-of-dynamics-under-external-forces.pdf
    # this might explain what i'm seeing

    q1e_jagsat = _remove_repeats(_drop_zeros(np.array(df_jagsat['q1e'])))[14:200]
    q2e_jagsat = _remove_repeats(_drop_zeros(np.array(df_jagsat['q2e'])))[14:200]
    q3e_jagsat = _remove_repeats(_drop_zeros(np.array(df_jagsat['q3e'])))[14:200]

    q1e_spec = _remove_repeats(_drop_zeros(np.array(df_spec['q1e'])))
    q2e_spec = _remove_repeats(_drop_zeros(np.array(df_spec['q2e'])))
    q3e_spec = _remove_repeats(_drop_zeros(np.array(df_spec['q3e'])))


    tt, ll = _compute_lyapunovs(q1e_spec, 0.005, 10)
    plt.plot(tt, ll)
    plt.grid()
    plt.show()

    tt, ll = _compute_lyapunovs(q1e_jagsat, 0.005, 1)
    plt.plot(tt, ll)
    plt.grid()
    plt.show()

    plot_parameter_space([q1e_jagsat], 3, 1)
    plot_parameter_space([q1e_spec], 3, 100)

    import nolds

    print('\n')
    print('jagsat controls: ', nolds.lyap_r(q1e_jagsat))
    print('spec controls: ', nolds.lyap_r(q1e_spec))
    print('\n')