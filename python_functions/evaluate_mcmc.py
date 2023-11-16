import matplotlib.pyplot as plt
import numpy as np


#root_path = "final_results_max_accuracy_max_iteration/"
#secnearios = ["illite_medium_goethite_fine", "illite_medium_goethite_coarse", "illite_fine_goethite_fine", "illite_fine_goethite_coarse"]
secnearios = ["result_mcmcm"]
#root_path = "result_mcmc_illite_medium_goethite_fine/"
#root_path = "result_mcmc_illite_medium_goethite_coarse/"
#root_path = "result_mcmc_illite_fine_goethite_fine/"
#root_path = "result_mcmc_illite_fine_goethite_coarse/"

for secnario in secnearios:

    root_path = "result_mcmc_3D_nx_50_v2/" 
    fig,ax = plt.subplots()

    dist = [ 0.25,0.5,0.75]
    jp = [0,5,10,15,20]
    nof_configs = len(dist) * len(jp)
    parameter_calc = np.zeros((nof_configs,2))
    parameter_east = np.zeros((nof_configs,2))
    labels = []
    cluster = ['^','.']
    #print(parameter_calc)
    index = 0;
    for d in dist:
        for j in jp:
            parameter_calc[index][0] = d
            parameter_calc[index][1] = j
            index = index +1
    #print(parameter_calc)
    index = 0;
    for d in dist:
        for j in jp:
            #filepath = "result_mcmcm_distribution_" + str(d) +'_jump_parameter_' + str(j)+ '_n_steps_10/'

            #filepath = "illite_medium_goethite_fine_distribution_" + str(d) +'_jump_parameter_' + str(j)+ '_n_steps_10/'
            #filepath = "illite_medium_goethite_coarse_distribution_" + str(d) +'_jump_parameter_' + str(j)+ '_n_steps_10/'
            #filepath = "illite_fine_goethite_fine_distribution_" + str(d) +'_jump_parameter_' + str(j)+ '_n_steps_10/'
            filepath = secnario + "_distribution_" + str(d) +'_jump_parameter_' + str(j)+ '_n_steps_10/'


            filepath =  root_path + "/" + filepath + "chainfile.txt"
            with open(filepath) as f:
                for line in f:
                    pass
                last_line = line.split(' ')
                parameter_east[index][0] = float(last_line[1])
                parameter_east[index][1] = float(last_line[0])
            labels.append( "calc [" + str(round(d,2)) + "," + str(round(j,2)) + "] east " + "[" + str(round(parameter_east[index][0],2)) + "," + str(round(parameter_east[index][1],2)) + "]")
            x = [parameter_east[index][0], parameter_calc[index][0]]
            y = [parameter_east[index][1], parameter_calc[index][1]]

            plt.plot(x,y, linestyle='dotted', label = labels[-1])
            for xp, yp, m in zip(x, y, cluster):
                plt.scatter([xp],[yp], marker=m)
            
            #print(last_line)
            #print(float(last_line[0]), float(last_line[1]))
            index = index +1
    #ax.legend(labels)
    plt.suptitle(root_path, fontsize = 12)
    ax.set_xlabel('distribution')
    ax.set_ylabel('jump parameter')
    minor_ticks_x = np.arange(-0.2, 1.2, 0.05)
    ax.set_xticks(minor_ticks_x, minor=True)

    minor_ticks_y = np.arange(-20, 50, 5)
    ax.set_yticks(minor_ticks_y, minor=True)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))



    plt.axis([-0.2, 1.2, -20, 50])
    ax.grid(which='minor', alpha=0.5)
    plt.savefig(root_path + "result_secnarios", dpi = 400)
    plt.show()
    #   print(labels)

