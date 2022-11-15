from __future__ import print_function

from datetime import datetime

import os, sys


# --------------------------------------------------------------------------------------------------
# Function bilaplacian_test.
# --------------------------------------------------------------------------------------------------
def cam_test(debug_mode=False):
  start_time = datetime.now()
  print("Starting time is", start_time)

  try:
    import CAM
  except (ImportError, ModuleNotFoundError) as error:
    sys.path.append(os.path.dirname(os.path.abspath(__file__)) + os.sep  + ".." + os.sep + "import")
    import CAM
  
  const                 = CAM.config()
  const.nx              = [10, 10]
  const.debug_mode      = debug_mode

  PyCAM = CAM.include(const)
  CAM_wrapper = PyCAM( 0.7, 5 )

  helper = CAM_wrapper.move_particles()
  CAM_wrapper.print_array(helper)
  
  end_time = datetime.now()
  print("Program ended at", end_time, "after", end_time-start_time)
  

# --------------------------------------------------------------------------------------------------
# Function main.
# --------------------------------------------------------------------------------------------------
def main(debug_mode):
  cam_test(debug_mode)


# --------------------------------------------------------------------------------------------------
# Define main function.
# -------------------------------------------------------------------------------------------------- 
if __name__ == "__main__":
  main(len(sys.argv) > 1 and sys.argv[1] == "True")
