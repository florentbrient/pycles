__author__ = 'pressel'

from distutils.core import setup
import distutils.core
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np
import mpi4py as mpi4py



#Get platform specific compile and linking arguments. Will use distutils.core.sys.platform to determine the system
#architecture
extra_compile_args=[]
extra_link_args=[]
platform = distutils.core.sys.platform
print(platform)
if platform == 'linux2':
    extra_compile_args+=['-Wno-maybe-uninitialized','-O3','-Wno-cpp','-Wno-array-bounds','-march=native','-Wno-unused','-Wno-#warnings','-fopenmp','-fPIC']
    extra_link_args+=['-fopenmp','-Wno-array-bounds']
#elif platform =='linux2':
#    extra_compile_args+=['-O3','-openmp']
#    extra_link_args+=['-openmp']


#Now get include paths from relevant python modules
include_path = [np.get_include()]
include_path += [mpi4py.get_include()]

#Now actually configure the extension build process
extensions = [
    Extension("*", ["*.pyx"],
        include_dirs = include_path,
        libraries = ["gfortran"],
        library_dirs = ["/cluster/apps/gcc/4.8.2/lib64/"], # Current on Euler as of Feb 01 2016
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args, 
        #extra_objects=['./RRTMG/build/rrtmg_lw.o','./RRTMG/build/rrtmg_sw.o'])
        extra_objects=['./RRTMG/build/rrtmg_combined.o'])
        ]
setup(
    name = "My hello app",
    ext_modules = cythonize(extensions),
)
