from setuptools import setup, find_packages #find_namespace_packages

def readme():
    with open('README.md') as f:
        return f.read()

setup(name='tasmanian-mismatch',
      version='1.0.8',
      description='Tasmanian tool to analyze mismatches at read and position in high throughput sequencing data',
      long_description=readme(),
      long_description_content_type="text/markdown",
      url='https://github.com/nebiolabs/tasmanian-mismatch',
      author='Ariel Erijman and Brad Langhorst',
      author_email='aerijman@neb.com',
      license='GNU',
      packages=find_packages(), #['tasmanian'],
      install_requires=[
          'numpy',
          'pandas',
          'scipy',
          'plotly'
      ],
      zip_safe=False,
      python_requires='>=3.10,<3.13',
      scripts=[
           'bin/run_tasmanian',
           'bin/run_intersections'
      ])
