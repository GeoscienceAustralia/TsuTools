from numpy.distutils.core import Extension, setup

setup(name='DC3D',
       ext_modules=[Extension(name='DC3D', sources=['TsuSource/DC3D/DC3D.f'])],
       )
