from setuptools import setup, find_packages

setup(
    name='MarkovPatternTool',
    version='0.2',
    author='Florian Mock',
    py_modules=['cli','FractalMatrix'],
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'Click',
        'numpy',
        'pyparsing',
        'pytz',
        'six',
        'matplotlib',
        'editdistance',
    ],
    entry_points={
        'console_scripts': [
            'MPT =cli:cli'
        ]}

)


'''
> MANIFEST.in tells Distutils what files to include in the source distribution
but it does not directly affect what files are installed. For that you need to
include the appropriate files in the setup.py file, generally either as package
data or as additional files. -- stackoverflow, 3596979

https://docs.python.org/3/distutils/setupscript.html#installing-package-data
https://docs.python.org/3/distutils/sourcedist.html#manifest
'''
