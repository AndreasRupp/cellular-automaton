import os, re

## Return path to directory of file.
def this_dir():
  path = os.path.dirname(os.path.abspath(__file__))
  path  =  re.sub(r'\\', '/', path)
  return path

## Return path to main directory of HyperHDG.
def main_dir():
  return this_dir() + "/../.."
