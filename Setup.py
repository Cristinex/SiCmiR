from setuptools import setup, find_packages

setup(
    name='SiCmiR',
    version='0.1.0',
    author='Xiao-Xuan Cai et al.',
    description='SiCmiR: Predict miRNA expression from L1000 landmark genes',
    packages=find_packages(where='src'),  
    package_dir={'': 'src'},
    install_requires=[
        'numpy',
        'pandas',
        'scipy',
        'torch>=2.0'
    ],
    entry_points={
        'console_scripts': [
            'sicmir=SiCmiR.SiCmiR:main'
        ]
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
    python_requires='>=3.8',
)
