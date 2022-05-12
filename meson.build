project('gefera', 'fortran', version : '0.1')

py3_mod = import('python3')
py3 = py3_mod.find_python('python3.8')
py3_dep = dependency('python3', required : false)

photlib = py3_mod.extension_module('phot',
	['fortran/phot.f90', 'fortran/phot_nograd.f90', 'fortran/ellip.f90'],
    	dependencies : py3_dep,
)
hrchlib = py3_mod.extension_module('hrch',
	['fortran/hrch.f90', 'fortran/kep.f90'],
	dependencies : py3_dep
)
conflib = py3_mod.extension_module('conf',
	['fortran/conf.f90', 'fortran/kep.f90'],
	dependencies : py3_dep	
)
