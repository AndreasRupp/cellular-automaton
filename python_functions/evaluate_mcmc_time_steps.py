import matplotlib.pyplot as plt
import numpy as np
import math

fig,ax = plt.subplots()#figsize=(12, 8)
plt.rcParams["figure.figsize"] = [7.50, 3.50]
plt.rcParams["figure.autolayout"] = True
root_path = "result_mcmc_3D_nx_50_v2/"
n_steps = []
for i in range(10):
    n_steps.append(i+1)
# for i in range(55):
#     ii = i +1
#     n_steps.append(10*ii)
nof_configs = len(n_steps)
parameter_calc = np.zeros((nof_configs,2))
parameter_east = np.zeros((nof_configs,2))
sqrt_error = np.zeros(nof_configs)
labels = []
cluster = ['^','.']
#print(parameter_calc)
index = 0;
for index in range(len(n_steps)):
        parameter_calc[index][0] = 0.5
        parameter_calc[index][1] = 20
#print(parameter_calc)
index = 0;
for n in n_steps:
        filepath = "result_mcmcm_illite_medium_goethite_coarse_n_steps_"+str(n)+'/'
        filepath =  root_path + filepath + "chainfile.txt"
        with open(filepath) as f:
            for line in f:
                pass
            last_line = line.split(' ')
            parameter_east[index][0] = float(last_line[1])
            parameter_east[index][1] = float(last_line[0])
        
        x = [parameter_east[index][0], parameter_calc[index][0]]
        y = [parameter_east[index][1], parameter_calc[index][1]]
        sqrt_error[index] = math.sqrt((parameter_east[index][0]- parameter_calc[index][0])**2 + (parameter_east[index][1]- parameter_calc[index][1])**2)
        
        labels.append( "n_steps: "+ str(n)+ " east " + "[" + str(round(parameter_east[index][0],2)) + "," + str(round(parameter_east[index][1],2)) + "] error " + str(round(sqrt_error[index], 2)))
        
        #plt.plot(x,y, linestyle='dotted', label = labels[-1]) #v1
        #plt.plot(n,sqrt_error[index], linestyle='dotted', label = labels[-1])

        plt.scatter(n,sqrt_error[index], label = labels[-1]) #v2
        ax.annotate(str(n), (n, sqrt_error[index]), fontsize=6) #v2

        # for xp, yp, m in zip(x, y, cluster): #v1
        #     plt.scatter([xp],[yp], marker=m) #v1
        #     ax.annotate(str(n), (xp, yp), fontsize=6) #v1
        
        #print(last_line)
        #print(float(last_line[0]), float(last_line[1]))
        index = index +1
#----plot1-------
# plt.suptitle(root_path, fontsize = 12)
# ax.set_xlabel('distribution')
# ax.set_ylabel('jump parameter')
# minor_ticks_x = np.arange(0.2, 0.7, 0.05)
# ax.set_xticks(minor_ticks_x, minor=True)

# minor_ticks_y = np.arange(-15, 40, 5)
# ax.set_yticks(minor_ticks_y, minor=True)
# box = ax.get_position()
# ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# # Put a legend to the right of the current axis
# ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': 8})



# plt.axis([0.2, 0.7, 15, 40])
# ax.grid(which='minor', alpha=0.5)

# plt.savefig(root_path + "timesteps", dpi = 800)
# plt.show()


#----plot2-------
plt.suptitle(root_path, fontsize = 12)
ax.set_xlabel('number of time steps')
ax.set_ylabel('sqrt error')

minor_ticks_x = np.arange(-1, 15,1 )#600
ax.set_xticks(minor_ticks_x, minor=True)

minor_ticks_y = np.arange(-5, 40, 1)
ax.set_yticks(minor_ticks_y, minor=True)

box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),prop={'size': 8})
plt.xticks(np.arange(0, 15, 1))#600
plt.yticks(np.arange(-5, 40, 1))
ax.grid(which='minor', alpha=0.5)
plt.savefig(root_path + "timesteps_error", dpi = 800)
plt.show()


