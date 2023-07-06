from __future__ import print_function
from datetime import datetime

import os, sys
import numpy as np
import itertools


def cam_particles(n_steps, jump_parameter, distribution,porosity, type_g, charges_g,type_i, charges_i, filename = 'output/cam.png', debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)

  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + "import")
    import CAM

  try:
    from plot import plot, plot_to_file, plot_to_vtk
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + 
      "python_functions")
    from plot import plot, plot_to_file, plot_to_vtk

  
  const                 = CAM.config()
  const.nx              = [500, 500]
  const.debug_mode      = debug_mode
  PyCAM = CAM.include(const)
  Domain = PyCAM(jump_parameter)

  n_fields = np.prod(const.nx)
  n_solids = round((1 - porosity) * n_fields)
  # print('number of fields',n_fields, 'number of solids', n_solids)



  placement_i = n_solids * distribution
  placement_g = n_solids * (1- distribution)
  print('Number cells of goethite', placement_g,'Number cells of illite', placement_i)
  print(type_g)
  print(type_i)
  size_g = np.prod(type_g)
  size_i = np.prod(type_i)
  n_g = round(placement_g/size_g)
  n_i = round(placement_i/size_i)
  # print('size g', size_g, 'size_i', size_i)
  print('number of goethite',n_g,'number of illite',n_i)

  def stencil_size(const_stencil, area):
    return np.floor((const_stencil/(np.sqrt(area))) + 0.5)
  #typ 1
  summe = 0
  success = 0
  jp_g = stencil_size(jump_parameter,size_g)
  while success < np.floor(n_g/2):
    success = success + Domain.place_plane(jp_g, type_g, -1, charges_g )#y #[1] * 4 
  print(success)
  summe = summe + success

  success = 0
  type_g = list(np.roll(type_g,1))
  charges_g = list(np.roll(charges_g,2))
  while success <  np.ceil(n_g/2):
    success = success + Domain.place_plane(jp_g,type_g,  -1, charges_g )#x #[1] * 4 
  print(success)
  summe = summe + success

  # placement_i = n_solids - summe * size_g
  # n_i = np.floor(placement_i/size_i)


  #type 2
  jp_i = stencil_size(jump_parameter,size_i)
  success = 0
  while success <  np.floor(n_i/2):
    success = success + Domain.place_plane(jp_i, type_i, -1,   charges_i)#y_neg #[-1] * 4
  print(success)
  summe = summe + success

  success = 0
  type_i = list(np.roll(type_i,1))
  charges_i = list(np.roll(charges_i,2))
  while success < np.ceil(n_i/2):
    success = success + Domain.place_plane(jp_i,type_i ,-1,charges_i)#x_neg #[-1] * 4
  print(success)
  summe = summe + success
  print('porosity ', (n_g * size_g + n_i * size_i)/n_fields)
  print('jump_parameter goethite', jp_g,'jump_parameter illite', jp_i)

  # Domain.print_array()

  # Domain.place_particles()

  # save_data = np.zeros( (n_steps + 1, np.prod(const.nx)) ) 
  plot_data = np.zeros( (2, np.prod(const.nx)) ) 

  # save_data[0] = Domain.fields()
  plot_data[0] = Domain.fields()
  # plot_to_vtk("output/cam", [save_data[0]], const.nx)
  for step in range(n_steps):
    Domain.do_cam()
    # save_data[step+1] = Domain.fields()
    print(step)
  plot_data[1] = Domain.fields()
  eval = Domain.eval_measures()
  # for i in range(12):
  #   print(eval[i])
  print('specific_surface_area',eval[7])
  print('contact_area_per_volume', eval[8])
  # print(Domain.average_particle_size())
  # print(Domain.bulk_distance(save_data[0],save_data[-1]))
  end_time = datetime.now() 
  print("Program ended at", end_time, "after", end_time-start_time)
  # save_data[(save_data <= n_g) & (save_data > 0)] = 1
  # save_data[save_data > n_g] = 2

  plot_data[(plot_data <= n_g) & (plot_data > 0)] = 1
  plot_data[plot_data > n_g] = 2
  # if not os.path.exists('output'):  os.makedirs('output')
  # plot_to_vtk("output/cam", [save_data[0]], const.nx)
  # if not os.path.exists('output2D'):  os.makedirs('output2D')
  # plot_to_vtk("output/cam", save_data, const.nx)
  # if not os.path.exists('output1'):  os.makedirs('output1')
  # plot_to_vtk("output1/cam", [save_data[-1]], const.nx)
  print("done")
  text = 'SSA ' + str(round(eval[7],2)) + ' (L^-1), CA/V ' + str(round(eval[8],2)) + ' (L^-1)'
  print(filename)
  plot_to_file(const.nx, plot_data[-1], filename, text)
  # plot(const.nx, plot_data, 0)
# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  n_steps =5000
  porosity = 0.9
  jump_const = 20
  type_i_addition = 0.95


  illite_fine = [2, 6]
  illite_medium = [2, 30]
  goethite_fine = [1, 17]
  goethite_coarse = [2, 34]

  uniform_positive = [1] * 4 
  uniform_negative = [-1] * 4 

  face_to_edge = [1,1, -1, -1]
  face_to_face = [-1,-1, -1, -1]
  #goethite positive charges
  #basal of illite = long edges
  #illite negative charges on basal planes and pos or negative charges on other egdes

  y = [0 ,0 , 1 ,1]
  y_neg = [0 ,0 , -1 ,-1]
  x =  [ 1,1,0,0]
  x_neg =  [ -1,-1,0,0]

  type_i_addition = 1- 0.05
  filename = 'output/goethite_illite_5.png'
  cam_particles(n_steps,jump_const,type_i_addition, porosity,  goethite_coarse,uniform_positive, illite_fine,uniform_negative ,filename,debug_mode)
  type_i_addition = 1 - 0.45
  filename = 'output/goethite_illite_45.png'
  cam_particles(n_steps,jump_const,type_i_addition, porosity,  goethite_coarse,uniform_positive, illite_fine,uniform_negative ,filename,debug_mode)
  type_i_addition = 1 - 0.55
  filename = 'output/goethite_illite_55.png'
  cam_particles(n_steps,jump_const,type_i_addition, porosity,  goethite_coarse,uniform_positive, illite_fine,uniform_negative ,filename,debug_mode)
  type_i_addition = 1 - 0.95
  filename = 'output/goethite_illite_95.png'
  cam_particles(n_steps,jump_const,type_i_addition, porosity,  goethite_coarse,uniform_positive, illite_fine,uniform_negative ,filename,debug_mode)



  type_i_addition = 0.5
  filename = 'output/illite_medium_edge_to_face.png'
  cam_particles(n_steps,jump_const,type_i_addition, porosity,  illite_medium, face_to_edge, illite_medium, face_to_edge ,filename, debug_mode)
  filename = 'output/illite_fine_edge_to_face.png'
  cam_particles(n_steps,jump_const,type_i_addition, porosity,  illite_fine, face_to_edge, illite_fine, face_to_edge ,filename, debug_mode)

  filename = 'output/illite_medium_face_to_face.png'
  cam_particles(n_steps,jump_const,type_i_addition, porosity,  illite_medium, face_to_face, illite_medium, [i * (-1) for i in face_to_face],filename,debug_mode)
  filename = 'output/illite_fine_face_to_face.png'
  cam_particles(n_steps,jump_const,type_i_addition, porosity,  illite_fine, face_to_face, illite_fine, [i * (-1) for i in face_to_face] ,filename, debug_mode)




# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
