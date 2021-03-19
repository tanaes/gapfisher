from setuptools import setup, find_packages

setup(
    name='gapfisher',
    version='0.1.0',
    description='A project to use Adaptive Sampling to target metagenomic assembly gaps ',
    url='https://github.com/tanaes/gapfisher',
    author='Jon G Sanders',
    author_email='jonsan@gmail.com',
    license='MIT License',
    packages=find_packages(),
    install_requires=['numpy',
                      'scikit-bio',
                      'click'],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Linux',
        'Operating System :: MacOS :: MacOS X'
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
    ],
)