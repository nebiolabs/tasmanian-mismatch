from setuptools import setup, find_packages #find_namespace_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='tasmanian missmatch analysis tool',
      version='0.0.29',
      description='Tasmanian tool to analyze mismatches at read and position in high throughput sequencing data',
      url='https://github.com/nebiolabs/tasmanian',
      author='Ariel Erijman and Brad Langhorst',
      author_email='aerijman@neb.com',
      license='MIT',
      packages=find_packages(), #['tasmanian'],
      install_requires=[
          'numpy',
          'pandas',
          'scipy',
          'plotly'
      ],
      zip_safe=False,
      python_requires='>=3.6',
      scripts=[
           'bin/run_tasmanian',
           'bin/run_intersections'
      ])
