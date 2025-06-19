from setuptools import setup, find_packages

setup(
    name='wavetrace-guinea',
    version='4.0.4',
    author='Justin Black',
    author_email='justinblack@college.harvard.edu',
    packages=find_packages(exclude=('tests', 'docs')),
    url='https://github.com/jblack020/wavetrace-guinea',
    license=license,
    data_files=[('', ['LICENSE.txt'])],
    description='A fork of wavetrace to produce radio signal coverage reports for Guinea instead of New Zealand',
    long_description_content_type='text/x-rst',
    install_requires=[
        'requests>=2.20.0',
        'Shapely>=1.6.4.post2',
        'click>=7.0',
    ],
    entry_points={
        'console_scripts': ['wavey=wavetrace.cli:wavey'],
    },
    python_requires='>=3.5',
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Topic :: Scientific/Engineering :: GIS',
    ],
)
