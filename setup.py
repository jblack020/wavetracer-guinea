from setuptools import setup, find_packages

setup(
    name='wavetracer-guinea',
    version='4.0.4',
    author='Justin Black',
    author_email='justinblack@college.harvard.edu',
    packages=find_packages(exclude=('tests', 'docs')),
    url='https://github.com/jblack020/wavetrace-guinea',
    description='A fork of wavetrace to produce radio signal coverage reports for Guinea instead of New Zealand',
    install_requires=[
        'requests>=2.20.0',
        'Shapely>=1.6.4.post2',
        'click>=7.0',
    ],
    entry_points={
        'console_scripts': ['wavey=wavetrace.cli:wavey'],
    },
)
