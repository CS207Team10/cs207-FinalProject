from setuptools import setup

setup(name='chemkin_g10',
      version='0.3',
      description='A chemical kinetics library',
      url='https://github.com/CS207Team10/cs207-FinalProject',
      author='Hidenori Tanaka, Jiachen Song, Xiangru Shu',
      author_email="jsong1@g.harvard.edu",
      license="MIT",
      packages=['chemkin_g10'],
      install_requires=[
          'numpy',
          'scipy'
      ],
      setup_requires=['pytest-runner'],
      tests_require=['pytest', 'pytest-cov'],
      include_package_data=True,
      )