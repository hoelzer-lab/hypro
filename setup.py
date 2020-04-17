from setuptools import setup, find_packages

requirements=['pandas==0.25.2','mygene==3.1.0'] # delete versions: 'pandas==0.25.2', 'mygene==3.1.0'

long_description="HyPro uses mmseqs2 to search for homology of sequence segments that have been annotated by Prokka as hypothetical proteins. These sequences are searched in a nucleotide database to obtain a more specific Prokka annotation. In this way, HyPro updates the Prokka output, as the program is designed to be easily integrated into custom analysis pipelines."
setup(
    name="HyPro",
    version='1.0',
    author='Maximilian Arlt and Martin Hoelzer',
    author_email='hoelzer.martin@gmail.com',
    license='GPLv3',
    description='Extension for homology searches of hypothetical proteins to enhance Prokka annotations',
    install_requires=requirements,
    scripts=['hypro-1.0/scripts/hypro.py', 'hypro-1.0/scripts/mmseqs2.sh'],
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/hoelzer-lab/hypro.git',       
    packages=find_packages(), # 're', 'os', 'sys', 'argparse', subprocess ? # change: install_requires
    classifiers=['Programming Language :: Python :: 3.7.6',  
    'Operating System :: Unix',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    ]
)

