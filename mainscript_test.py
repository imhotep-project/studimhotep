"""Script to test the package
This the main script to test the package is well installed. 
It calls functions from  my library of tools in libimhotep/libSLXclassIMHOTEP.py
"""

from libimhotep import libSLXtoolsIMHOTEP as slx

print('This the main script to test the package is well installed. If so, it should output a global map without data.')


print(slx.main())