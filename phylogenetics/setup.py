from setuptools import setup

setup(
    name='phylogenetics',
    version='0.1',
    packages=[''],
    url='',
    license='MIT',
    author='smurphy333',
    author_email='steven.t.murphy3@gmail.com',
    description='A tool to generate phylogenetic data',
    install_requires=['biopython', 'matplotlib'],
    zip_safe=False,
    entry_points={
        'make_phylogenetic_tree.main': [
            'make_tree = make_phylogenetic_tree:main',
        ]
    }
)
